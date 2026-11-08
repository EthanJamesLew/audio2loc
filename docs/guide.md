# audio2loc guide

Geolocates an audio or video recording from the Doppler curves of propeller aircraft passing
overhead, matched against historical ADS-B flight data. It's a C++20 reimplementation of the method in
colsto's *"I discovered a sound that leaks your exact location"* (youtube.com/watch?v=UytYrhYLOKI).

The core has no dependencies beyond zlib: its own FFT, WAV reader, JSON reader, PNG writer and SRTM
reader. Video and other compressed audio are decoded through `ffmpeg` if it's installed. Results come
out as GeoJSON plus a ready-to-open QGIS project (`audio2loc.qgs`). Optional builds link
`libqgis_core` to write a `.qgz` through the native QGIS API, and run the heavy stages on an NVIDIA GPU
(CUDA). It also ships as a Docker image that runs an MCP server, so an AI agent can drive it.

## Build

```sh
cmake -S . -B build && cmake --build build -j
# optional native QGIS link (needs QGIS dev headers, Qt5 or Qt6):
cmake -S . -B build -DAUDIO2LOC_WITH_QGIS=ON [-DQGIS_PREFIX=/usr] && cmake --build build -j
# optional CUDA (needs the CUDA toolkit; nvcc is found in /usr/local/cuda automatically):
cmake -S . -B build -DAUDIO2LOC_CUDA=ON [-DCMAKE_CUDA_ARCHITECTURES="75;86;89"] && cmake --build build -j
```

A CUDA build uses the GPU automatically and falls back to the CPU if there's no device or a GPU call
fails. `--no-gpu` forces the CPU path, which stays the reference implementation.

## Workflow

```sh
# 1. find the Doppler curves (zoom into roughly 40-400 Hz; they're often inaudible but visible)
audio2loc spectrogram --audio clip.mp4 --t0 0 --t1 600 --fmin 40 --fmax 400 -o spec.png

# 2. put a box around one tone of each curve and check the fit (speed, CPA distance, CPA time, f0)
audio2loc measure --audio clip.mp4 --event t0=95,t1=210,fmin=95,fmax=170,mult=1:2:3 \
                                  --event t0=1010,t1=1130,fmin=70,fmax=120

# 3. ADS-B for the day (streams adsb.lol's global daily history into ./audio2loc_db)
audio2loc update --date 2025-09-08 --bbox -125,24,-66,50

# 4. search (--auto-update fetches any missing day in the window)
audio2loc locate --audio clip.mp4 --db audio2loc_db --auto-update \
  --event t0=95,t1=210,fmin=95,fmax=170,mult=1:2:3 --event t0=1010,t1=1130,fmin=70,fmax=120 \
  --start 2025-09-08 --end 2025-09-09 --dem srtm/ --out result/
qgis result/audio2loc.qgs
```

If you know roughly when the recording started, `--rec-start ISO --rec-start-tol SEC` narrows the
search a lot. With a tight time you can locate from a single curve, as in the video's second example.
`--bbox W,S,E,N` restricts the area.

## Method (maps one-to-one to the video)

| step | what | code |
|---|---|---|
| measure | Whitened STFT, then a continuity-constrained Viterbi ridge. A straight-line Doppler model with **retarded time** is fitted to it. The fit's uncertainty region gives generous ranges for speed, slant CPA distance, CPA time and f0. | `doppler.cpp`, `search.cpp:prepare_events` |
| prop filter | Known ICAO jet types are dropped and everything else is kept, since there's no reliable prop flag. | `aircraft_types.cpp` |
| A: candidates | A 200 m grid is laid around speed-matching flight segments. Each node gets the exact CPA (Hermite trajectory, golden-section refinement) using **slant** distance to the DEM ground height. The node is kept if that distance is in the measured range, and it records the **delayed CPA** time (CPA + d/c). | `gen_hits` |
| B: time offset | Node must hear event *k* at (tc_k − tc_1) ± tolerance after event 1 | join in `run_locate` |
| C: silence | Node must not hear any other aircraft between the events (radius 2 km for props, 4 km for jets) | `SilenceIndex` |
| D: full model | On a 100 m grid, the curve each node would have heard is computed from the real ADS-B trajectory. The radial velocity is taken at the emission time, not the reception time. Each curve is scored by mean spectrogram brightness under it (log-frequency image, so f0 and harmonics become row shifts). Event 1 sweeps ±3 s and f0 to absorb measurement error. The other events' timing is **locked** to event 1 by flight geometry; only their f0 is swept. | `Scorer` |
| E: rank | Best node per flight combination, ranked. Then a dense heatmap (25 m) around the winner, using only its flights. | `run_locate` |

## Keeping the ADS-B database current

```sh
audio2loc update --latest                       # newest published day (normally yesterday, UTC)
audio2loc update --from 2026-11-01 --to 2026-11-06 --bbox -125,24,-66,50
audio2loc update --list                         # what's in the database
```

`update` finds the adsb.lol release for each UTC day on GitHub and streams its split tar archive
(about 3.4 GB) through a built-in tar parser and gzip decoder straight into `adsb_YYYY-MM-DD.plb`.
Nothing is extracted to disk, and memory stays constant. A continental-US day took 8 minutes
(network-bound at about 7 MB/s) and produced a 1.6 GB file with 96k flight legs and 41M samples.

* Day files carry a `.meta` sidecar recording their region filter. A search outside that region
  prints a warning instead of silently finding nothing.
* Recency: adsb.lol publishes each UTC day at about 03:30 UTC the next day, so the newest searchable
  day is yesterday. For same-day recordings, export ADS-B yourself (for example from OpenSky) as CSV.
* Without `GITHUB_TOKEN`, the GitHub API allows 60 requests per hour. Set the variable for heavy use.

## GPU acceleration (CUDA)

| stage | CPU path | GPU path |
|---|---|---|
| spectrogram | built-in FFT | batched cuFFT |
| measure fit | (v, d, tc) grid on threads | one thread per grid cell |
| A: candidates | per flight: rasterize, sort, scan for minima | one thread per (sample, grid row). The minimum test uses neighbouring samples, so no sort is needed. Results are sorted and deduplicated on the device with Thrust. |
| D: Doppler modelling | one candidate per thread; SIMD row-vector alignment sweep | one 256-thread block per candidate: FP32 local trajectory in shared memory, emission-time sampling, coalesced frame-outer alignment sweep on harmonic-summed tables |

The GPU kernels follow the CPU algorithms step by step. On the validation scenarios, hit and
candidate counts are identical and the top rankings match to 3 decimals. On ambiguous decoy nodes,
near-tied alignment sweeps can resolve differently: about 3% of individual node scores differ by
more than 0.001.

Measured on an RTX 3090 vs 32 CPU threads (Ryzen 9 5950X). Real ADS-B for 2026-11-06 (184k flight legs), a hidden
observer near Davis, CA, and audio synthesized from two real flights. Search region 4.5° × 4°:

| | CPU | GPU |
|---|---|---|
| A: candidates (14.2M hits) | 10.7 s | 3.3 s |
| B + C: join, silence filter | 2.8 s | 0.9 s |
| D: Doppler modelling (1.53M locations) | 27.2 s | 3.8 s |
| **total** | **42.1 s** | **9.0 s** |
| result | 39 m from truth, 25% margin | identical |

The main algorithmic speedups (both paths):

* **Harmonic tables.** The model point's fractional row selects one of at most nm+1 integer offset
  patterns, so harmonics are summed once per pattern up front. This gives one table read per model
  point instead of nm, with exactly the same scores.
* **Separable sweep.** For a fixed time shift, every f0 offset shares the same per-frame base row, so
  the f0 sweep is a contiguous vector add (AVX on the CPU, coalesced loads on the GPU). The coarse
  grid uses a row-phase-split table.
* **Emission-time sampling.** The Doppler model is evaluated once per emission sample and
  interpolated to frames, instead of solving retarded time per frame.
* **CPA root-finding.** Closest approach is found by regula falsi on v·(p−o): about 6 evaluations
  instead of 32 golden-section steps.
* **Smaller host stages.** B is a parallel merge-join, and the silence index is a CSR built only
  over cells the tile can query.

## Docker and the MCP server

```sh
docker build -t audio2loc .                     # slim CPU image (~0.8 GB)
docker build -t audio2loc:cuda --target cuda .  # CUDA image
docker build -t audio2loc:qgis --target qgis .  # native QGIS .qgz export

docker run --rm -v audio2loc-db:/data audio2loc update --latest --bbox -125,24,-66,50
docker run --rm -i -v audio2loc-db:/data -v "$PWD":/work audio2loc                  # MCP server (stdio)
docker run --rm -i --device nvidia.com/gpu=all -v audio2loc-db:/data -v "$PWD":/work audio2loc:cuda
```

`audio2loc mcp` is a Model Context Protocol server on stdio (JSON-RPC 2.0) with five tools:

* `spectrogram`: returns a PNG the model can look at to find curves
* `measure_doppler`: fits the boxed curves and returns overlays
* `locate`: full search; fetches missing ADS-B days and returns ranked locations, overlays, and an
  output directory with the QGIS project
* `update_database` and `database_status`

`audio` can be a path or a URL; URLs are downloaded with yt-dlp, which the image includes. Each call
runs in a subprocess, so a long or failing search can't take the server down. Example client config:

```json
{"mcpServers": {"audio2loc": {"command": "docker",
  "args": ["run", "--rm", "-i", "-v", "audio2loc-db:/data", "-v", "/path/to/media:/work", "audio2loc"]}}}
```

## Memory

A continental-US day is ~100k flight legs and, unconstrained, hundreds of millions of intermediate
(grid node x flight pass) hits. `locate` keeps memory bounded:

* **Hard limit** `--mem-limit GB`. Default: half of the RAM *available* at start, capped at 85% of any
  cgroup limit (`docker run --memory`, systemd scopes). Also settable with `AUDIO2LOC_MEM_LIMIT`.
  A watchdog thread exits cleanly with code 3 and a message if resident memory ever exceeds it, so a
  run can never take the machine (or the agent driving it) down.
* **Streaming tiles**: stages A–D run one spatial tile at a time (`--tile-deg`, default 2°). Candidates
  are generated, joined, silence-filtered and Doppler-scored per tile, and only the best location per
  flight combination is kept. A tile whose hits would exceed the budget is split into four and
  retried, so the cap holds regardless of traffic density, and no candidate is ever dropped.
* **Bounded scratch**: each flight is rasterized in fixed-size sample chunks, scoring runs in batches,
  and large layers are streamed to disk.
* **Lean inputs**: ADS-B is thinned at load to the 3 s step the trajectories use, and raw samples are
  freed once trajectories are built.

Measured: the full continental-US day that used to grow past 14 GB and crash ran at a steady
3.0–3.2 GB under an 8 GB limit, through 7M candidates and 24M scored locations. Under a 2.2 GB limit,
the California search split 4 tiles, peaked at 1.03 GB and gave the same answer. An impossible
0.3 GB limit exits with code 3 and a clear message.

## Outputs (`--out DIR`)

* `audio2loc.qgs`: QGIS project with an OSM basemap and these layers: best location, ranked
  candidates, refinement heatmap, flight paths, and each filtering stage
* `best.geojson`, `ranked.geojson`, `refine_rank1.geojson`, `flights.geojson`, `stage{A,B,C,D}_*.geojson`
* `overlay_rank{N}_event{K}.png`: spectrogram with the modelled curve on top, the "does it cover
  the curve" check from the video
* `result.json`: fits, ranges, winning flights, and the implied absolute start time of the recording

## ADS-B inputs

* CSV (optionally `.gz`) with columns `time` (epoch or ISO), `icao24`, `lat`, `lon`, an altitude column
  (`alt_m`, `alt_ft`, `geoaltitude` (m), `baroaltitude` (m), `alt_geom` (ft), `altitude` (ft)), and
  optionally `callsign`, `typecode`, speed (`gs_kt`, `gs_mps`, `velocity` (m/s), `groundspeed` (kt)),
  `track`, and vertical rate (`vrate_fpm`, `vertrate` (m/s), `vertical_rate` (ft/min)). Override the
  unit guesses with `--csv-alt-unit`/`--csv-speed-unit`.
* readsb / adsb.lol / ADSBExchange `trace_full_*.json` files (gzip-compressed or plain), or a
  directory of them
* `.plb`: audio2loc's own binary cache (`audio2loc ingest`)

Without the `traffic`/`pyopensky` Python libraries, the OpenSky "type code" comes from your export.
`--types FILE` (`typecode,jet|turboprop|piston|heli`) extends the built-in jet list.

## Validation

`audio2loc simulate` builds a ground-truth scenario. It produces 6 hours of synthetic traffic over a
1.2° region (300 props, 150 jets) and 30 decoy pairs that match the target's speeds, distances and
920 s offset but fly straight. The audio is synthesised at a hidden observer by evaluating each
source's phase at the retarded emission time, with 1/r spreading and noise. In the target pair,
one plane turns through its CPA and the other climbs and turns before CPA.

```
audio2loc simulate --prefix sim && audio2loc $(grep suggested sim_truth.txt | cut -d' ' -f3-)
```

Measured on this machine (32 threads):

| scenario | result | error | runtime |
|---|---|---|---|
| seed 7, two events | correct flight pair, 16% ahead of the best decoy pair | 11 m | 3.3 s |
| same, as MP4 with AAC audio at 96 kbps via ffmpeg | same answer | 11 m | |
| seed 99, noisier audio, DEM tiles | correct pair after the 25 m refinement | 10 m | 50 s |
| single event, `--rec-start` ±60 s, brightness scoring | correct | 11 m | |
| single event, `score=ridge` (auto-traced) | true flight ranked **2nd** | | |
| **real ADS-B** (2026-11-06, 184k legs), observer near Davis CA, 4.5° × 4° region | correct flight pair, 25% margin | 39 m | 9 s GPU / 42 s CPU |

`audio2loc simulate --adsb DAY.plb` builds the real-traffic scenarios. It hides an observer under a
real non-jet flight, finds the next real flight to pass within `--max-slant` (default 1200 m) with no
other aircraft in between, and synthesizes the audio from the real trajectories. With faint, distant
passes the automatic curve fit can be badly wrong: a 1.7 km pass was fitted at 8 km. Check the
`measure` overlays and use the `d=`, `v=` and `tc=` overrides when the fit doesn't follow the curve.

The ridge mode needs a careful trace. The automatic one loses the faint wings, which are what tell
flights apart. Hand-traced points can be supplied with `trace=FILE`.

## Limits / notes

* A single straight flyover is symmetric: a mirror point on the other side of the track hears the
  same curve. Two events, or a manoeuvring aircraft, or a tight time constraint breaks the symmetry.
* Doppler depends on velocity relative to the air. Wind shifts the measured speed, which is part of
  why the speed margins are generous (`--v-margin`).
* Altitude is taken as MSL: geometric altitude minus `--geoid` when available, otherwise
  barometric. Barometric altitude can be off by tens of metres, which matters mostly for very close
  passes.
* The `-DAUDIO2LOC_WITH_QGIS=ON` path is built and tested in the `qgis` Docker target (QGIS 3.34). The
  OSM basemap layer is omitted from the native `.qgz` there, because the XYZ provider plugin isn't in
  `libqgis-core`. The XML `.qgs` writer needs no QGIS.
* Full continental-US searches over a whole day with loose ranges are still slow (tens of millions of
  candidates). Narrowing the region (`--bbox`), the time window (`--rec-start`), or the measured
  ranges matters far more than hardware.

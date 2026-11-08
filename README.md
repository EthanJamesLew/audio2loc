# audio2loc

audio2loc finds where an audio or video recording was made from the Doppler curves of propeller
aircraft heard in it. You box each flyover's tone on a spectrogram, and audio2loc fits its speed,
closest-approach distance and time. It then searches a day of historical ADS-B tracks from
[adsb.lol](https://adsb.lol) for the spot on the ground that would have heard exactly those curves.
It reimplements the method from colsto's *"I discovered a sound that leaks your exact location"*
([video](https://youtube.com/watch?v=UytYrhYLOKI)) as a C++20 command-line tool. Results are ranked
coordinates, GeoJSON layers and a QGIS project. It runs on the CPU and, optionally, an NVIDIA GPU.

You need at least a minute of a clearly audible piston or turboprop flyover, ideally two, and the
date. Uploads to messengers and social media keep the audio but strip the date and location
metadata, so audio2loc works on those copies too. Use it to verify or investigate material whose
subject consents or that is of public interest, not to find private individuals.
[docs/guide.md](docs/guide.md) covers the method, every option and the benchmarks.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # needs CMake 3.24+, a C++20 compiler, zlib
cmake --build build -j
ctest --test-dir build                           # includes a synthetic end-to-end locate
build/audio2loc --help
```

Add `-DAUDIO2LOC_CUDA=ON` for the GPU path, which falls back to the CPU if there's no GPU, and
`-DAUDIO2LOC_NATIVE=ON` to tune for your CPU. `ffmpeg` decodes video and compressed audio.
To try the pipeline without real data:

```sh
build/audio2loc simulate --prefix sim --seed 7        # synthetic flights and audio, hidden observer
build/audio2loc $(grep ^suggested: sim_truth.txt | cut -d' ' -f3-)
```

## MCP server

`audio2loc mcp` is a Model Context Protocol server on stdio. Its tools are `spectrogram`,
`measure_doppler`, `locate`, `update_database` and `database_status`, so an AI agent can run the
whole workflow and download ADS-B days as needed. With Claude Code:

```sh
claude mcp add audio2loc -- $PWD/build/audio2loc mcp --work $HOME/audio2loc/work --db $HOME/audio2loc/db
```

Or with Docker, which also bundles ffmpeg and yt-dlp, so the tools accept URLs:

```sh
docker build -t audio2loc .                                       # or --target cuda
claude mcp add audio2loc -- docker run --rm -i -v audio2loc-db:/data -v "$PWD":/work audio2loc
```

Other clients take the same command as JSON:

```json
{ "mcpServers": { "audio2loc": { "command": "/path/to/audio2loc", "args": ["mcp", "--db", "/path/to/db"] } } }
```

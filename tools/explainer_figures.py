#!/usr/bin/env python3
"""Draws the README's algorithm figures from audio2loc's own outputs.

    audio2loc simulate --prefix sim --seed 7
    audio2loc spectrogram --audio sim.wav --t0 0 --t1 1400 --fmin 50 --fmax 250 --hop 0.5 -o spec_full.png
    audio2loc spectrogram --audio sim.wav --t0 180 --t1 310 --fmin 75 --fmax 170 --hop 0.1 -o spec_e1.png
    audio2loc spectrogram --audio sim.wav --t0 1100 --t1 1230 --fmin 60 --fmax 120 --hop 0.1 -o spec_e2.png
    audio2loc $(grep ^suggested: sim_truth.txt | cut -d' ' -f3-) --out out
    tools/explainer_figures.py <that directory> docs/images

Needs matplotlib, numpy and pillow. Writes <name>-light.png and <name>-dark.png.
"""
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from PIL import Image, ImageFilter

SRC = Path(sys.argv[1])
DST = Path(sys.argv[2])
DST.mkdir(parents=True, exist_ok=True)

# Spectrogram extents, as `audio2loc spectrogram` reports them.
SPECS = {
    "full": (0.0, 0.5, 250.0, 50.0),  # t0, s/px, top Hz, bottom Hz
    "e1": (180.0, 0.1, 170.0, 75.0),
    "e2": (1100.0, 0.1, 120.0, 60.0),
}
C = 331.3 * math.sqrt(1 + 20 / 273.15)

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#8a8984",
                  grid="#e6e5e0", s1="#2a78d6", s2="#eb6834", s3="#1baf7a",
                  track="#c9c8c2",
                  ramp=["#f3f7fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                        "#256abf", "#184f95", "#0d366b"]),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#8f8e86",
                 grid="#2e2e2c", s1="#3987e5", s2="#d95926", s3="#199e70",
                 track="#3d3d3a",
                 ramp=["#1f2329", "#0d366b", "#184f95", "#256abf", "#3987e5",
                       "#6da7ec", "#9ec5f4", "#cde2fb"]),
}


def straight_lng(t, v, d, tc, c=C):
    """ln(f/f0) of a straight constant-speed pass; mirrors src/doppler.cpp."""
    t0 = tc - d / c
    T = t - t0
    c2 = c * c
    a = c2 - v * v
    tau = (c2 * T - np.sqrt(c2 * v * v * T * T + a * d * d)) / a
    x = v * tau
    r = np.sqrt(x * x + d * d)
    return -np.log1p(v * x / r / c)


def points(name, key=None):
    js = json.loads((SRC / "out" / name).read_text())
    xs, ys, vs = [], [], []
    for f in js["features"]:
        lon, lat = f["geometry"]["coordinates"][:2]
        xs.append(lon)
        ys.append(lat)
        if key:
            vs.append(f["properties"][key])
    return np.array(xs), np.array(ys), np.array(vs)


def truth():
    t = {}
    for line in (SRC / "sim_truth.txt").read_text().splitlines():
        k, _, v = line.partition(" ")
        t[k] = v
    return float(t["observer_lat"]), float(t["observer_lon"])


def tracks():
    """Every simulated ADS-B track, downsampled, as (icao, lons, lats)."""
    by = {}
    with open(SRC / "sim.csv") as f:
        for i, row in enumerate(csv.DictReader(f)):
            if i % 4:
                continue
            by.setdefault(row["icao24"], ([], []))
            by[row["icao24"]][0].append(float(row["lon"]))
            by[row["icao24"]][1].append(float(row["lat"]))
    return by


RESULT = json.loads((SRC / "out" / "result.json").read_text())
_BEST_FEATURE = json.loads((SRC / "out" / "best.geojson").read_text())["features"][0]
_BEST_PROPS = _BEST_FEATURE["properties"]
RESULT_D = {1: _BEST_PROPS["d1"], 2: _BEST_PROPS["d2"]}  # modelled slant CPA (m)


def haversine(la1, lo1, la2, lo2):
    p1, p2 = math.radians(la1), math.radians(la2)
    a = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lo2 - lo1) / 2) ** 2)
    return 2 * 6371000 * math.asin(math.sqrt(a))
TRUE_LAT, TRUE_LON = truth()
BEST = dict(zip(("lon", "lat"), _BEST_FEATURE["geometry"]["coordinates"]))  # full precision
TRACKS = tracks()


def style(th):
    plt.rcParams.update({
        "figure.facecolor": th["surface"], "axes.facecolor": th["surface"],
        "savefig.facecolor": th["surface"], "text.color": th["ink"],
        "axes.labelcolor": th["ink2"], "axes.edgecolor": th["grid"],
        "xtick.color": th["ink2"], "ytick.color": th["ink2"],
        "axes.titlecolor": th["ink"], "font.size": 11, "axes.titlesize": 12.5,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.spines.top": False, "axes.spines.right": False,
        "font.family": "DejaVu Sans", "legend.frameon": False,
    })


def spec_image(key, ax, thicken=1):
    t0, sp, top, bottom = SPECS[key]
    im = Image.open(SRC / f"spec_{key}.png").convert("RGB")
    if thicken > 1:  # keep one-pixel ridges visible once the image is scaled down
        im = im.filter(ImageFilter.MaxFilter(thicken))
    w, _ = im.size
    ax.imshow(im, extent=(t0, t0 + w * sp, bottom, top), aspect="auto",
              interpolation="antialiased")
    ax.set_xlim(t0, t0 + w * sp)
    ax.set_ylim(bottom, top)


def save(fig, name, theme):
    fig.savefig(DST / f"{name}-{theme}.png", dpi=150, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def fig_spectrogram(th, theme):
    fig, ax = plt.subplots(figsize=(11, 4.2))
    spec_image("full", ax, thicken=3)
    for k, ev in enumerate(RESULT["events"]):
        t0, t1, f0, f1 = ev["box"]
        ax.add_patch(Rectangle((t0, f0), t1 - t0, f1 - f0, fill=False,
                               lw=1.6, ec="#ffffff", ls=(0, (4, 2))))
        ax.text(t0, f1 + 4, f"flyover {k + 1}", color="#ffffff", fontsize=11,
                fontweight="bold", va="bottom")
    ax.set_xlabel("audio time (s)")
    ax.set_ylabel("frequency (Hz)")
    ax.set_title("1  Find the flyovers: a whitened spectrogram of the 23-minute recording")
    save(fig, "01-spectrogram", theme)


def fig_model(th, theme):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
    t = np.linspace(-90, 90, 600)
    f0 = 100.0
    for ax, vary, vals, fixed in (
        (axes[0], "d", [300, 900, 2500], 60.0),
        (axes[1], "v", [40, 60, 90], 900.0),
    ):
        for col, val in zip((th["s1"], th["s2"], th["s3"]), vals):
            v, d = (fixed, val) if vary == "d" else (val, fixed)
            f = f0 * np.exp(straight_lng(t, v, d, 0.0))
            ax.plot(t, f, color=col, lw=2.2)
            lbl = f"d = {val:,} m" if vary == "d" else f"v = {val} m/s"
            ax.plot([], [], color=col, lw=2.2, label=lbl)
        ax.legend(loc="lower left", fontsize=10, handlelength=1.6, labelcolor=th["ink"])
        ax.axvline(0, color=th["muted"], lw=1, ls=(0, (3, 3)))
        ax.axhline(f0, color=th["muted"], lw=1, ls=(0, (3, 3)))
        ax.set_xlim(-90, 90)
        ax.set_xlabel("time from closest approach, as heard (s)")
        ax.grid(axis="y", color=th["grid"], lw=0.8)
        ax.set_title(("Closer passes turn faster" if vary == "d"
                      else "Faster planes shift further") + f"  ({'v = 60 m/s' if vary == 'd' else 'd = 900 m'})",
                     fontsize=11.5, fontweight="normal", color=th["ink2"])
    axes[0].set_ylabel("heard frequency for f₀ = 100 Hz")
    axes[0].text(2, 81.5, "t_c", color=th["ink2"], fontsize=10)
    axes[0].text(-88, f0 + 0.6, "f₀", color=th["ink2"], fontsize=10)
    fig.suptitle("2  The physics: every pass draws a curve fixed by speed v, closest distance d, time t_c and tone f₀",
                 x=0.0, ha="left", fontsize=12.5, fontweight="bold", color=th["ink"])
    fig.tight_layout()
    save(fig, "02-doppler-model", theme)


def fig_measure(th, theme):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for k, (ax, key) in enumerate(zip(axes, ("e1", "e2"))):
        spec_image(key, ax)
        ev = RESULT["events"][k]
        fit, rg = ev["fit"], ev["ranges"]
        t0, sp, _, _ = SPECS[key]
        t = np.linspace(t0, t0 + 130, 800)
        g = np.exp(straight_lng(t, fit["v"], fit["d"], fit["tc"]))
        ax.plot(t, fit["f0"] * g, color="#ffffff", lw=1.1, ls=(0, (4, 4)), alpha=0.9)
        ax.set_xlabel("audio time (s)")
        if k == 0:
            ax.set_ylabel("frequency (Hz)")
        txt = (f"v  = {fit['v']:.0f} m/s   ({rg['v'][0]:.0f}–{rg['v'][1]:.0f})\n"
               f"d  = {fit['d']:.0f} m   ({rg['d'][0]:.0f}–{rg['d'][1]:.0f})\n"
               f"t_c = {fit['tc']:.1f} s   (±{rg['tc_tol']:.1f})\n"
               f"f₀ = {fit['f0']:.1f} Hz")
        ax.text(0.03, 0.05, txt, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=9.5, family="DejaVu Sans Mono", color="#ffffff",
                bbox=dict(boxstyle="round,pad=0.45", fc="#000000", ec="none", alpha=0.55))
        ax.set_title(f"flyover {k + 1}", fontsize=11.5, fontweight="normal", color=th["ink2"])
    fig.suptitle("3  Measure: fit the model to each curve (dashed); the fit and its uncertainty ranges carry into the search",
                 x=0.0, ha="left", fontsize=12.5, fontweight="bold", color=th["ink"])
    fig.tight_layout()
    save(fig, "03-measure", theme)


def map_axes(ax, th, extent):
    lo0, lo1, la0, la1 = extent
    ax.set_xlim(lo0, lo1)
    ax.set_ylim(la0, la1)
    ax.set_aspect(1 / math.cos(math.radians((la0 + la1) / 2)))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(th["grid"])


def draw_tracks(ax, th, lw=0.5, alpha=1.0):
    for lons, lats in TRACKS.values():
        ax.plot(lons, lats, color=th["track"], lw=lw, alpha=alpha, solid_capstyle="round")


def truth_marker(ax, th, s=150):
    ax.scatter([TRUE_LON], [TRUE_LAT], marker="*", s=s, color=th["s2"],
               edgecolor=th["surface"], linewidth=1.2, zorder=10)


def fig_funnel(th, theme):
    xs = [TRUE_LON - 1.25, TRUE_LON + 1.25]
    ys = [TRUE_LAT - 0.95, TRUE_LAT + 0.95]
    extent = (*xs, *ys)
    bx, by, _ = points("stageB_offset.geojson")
    cx, cy, _ = points("stageC_silence.geojson")
    dx, dy, dv = points("stageD_scored.geojson", "total")
    panels = [
        ("A  CPA hits", "1,345,423 grid nodes hear a plane pass\nat the measured distance and speed"),
        ("B  time offset", "4,869 nodes also hear a second pass\n921 s later"),
        ("C  silence", "2,802 nodes heard no other\naircraft in between"),
        ("D  Doppler score", "12,731 node × flight-pair curves\nmodelled and scored on the audio"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13, 4.3))
    for ax, (title, sub) in zip(axes, panels):
        map_axes(ax, th, extent)
        ax.set_title(title, fontsize=12)
        ax.text(0.0, -0.02, sub, transform=ax.transAxes, va="top", fontsize=9.5, color=th["ink2"])
    draw_tracks(axes[0], th, lw=0.6)
    axes[0].text(0.03, 0.97, f"all ADS-B tracks\n(6 h, {len(TRACKS)} aircraft)", transform=axes[0].transAxes,
                 va="top", fontsize=9, color=th["ink2"])
    for ax in axes[1:]:
        draw_tracks(ax, th, lw=0.4, alpha=0.6)
    axes[1].scatter(bx, by, s=4, color=th["s1"], lw=0, alpha=0.8)
    axes[2].scatter(cx, cy, s=4, color=th["s1"], lw=0, alpha=0.8)
    order = np.argsort(dv)
    cmap = LinearSegmentedColormap.from_list("seq", th["ramp"][2:])
    axes[3].scatter(dx[order], dy[order], s=4, c=dv[order], cmap=cmap, lw=0)
    for ax in axes:
        truth_marker(ax, th, 110)
    axes[3].annotate("true location", (TRUE_LON, TRUE_LAT), xytext=(TRUE_LON + 0.25, TRUE_LAT + 0.55),
                     fontsize=9.5, color=th["ink"],
                     arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=0.8))
    fig.suptitle("4  Search: each filter is cheap and removes most of what is left, so the expensive model runs on very little",
                 x=0.0, ha="left", fontsize=12.5, fontweight="bold", color=th["ink"])
    fig.tight_layout()
    save(fig, "04-search", theme)


def fig_result(th, theme):
    fx, fy, fv = points("refine_rank1.geojson", "total")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), gridspec_kw=dict(width_ratios=[1, 1]))
    # Left: the two winning flights passing the observer.
    ax = axes[0]
    pad = 0.022
    extent = (TRUE_LON - pad / 0.76, TRUE_LON + pad / 0.76, TRUE_LAT - pad, TRUE_LAT + pad)
    map_axes(ax, th, extent)
    draw_tracks(ax, th, lw=0.8, alpha=0.7)
    flights = json.loads((SRC / "out" / "flights.geojson").read_text())["features"]
    mlat = 111_000
    mlon = 111_000 * math.cos(math.radians(TRUE_LAT))
    for f in flights:
        p = f["properties"]
        if p["rank"] != 1:
            continue
        col = th["s1"] if p["event"] == 1 else th["s3"]
        xy = np.array(f["geometry"]["coordinates"])
        ax.plot(xy[:, 0], xy[:, 1], color=col, lw=2.4)
        dxm = (xy[:, 0] - TRUE_LON) * mlon
        dym = (xy[:, 1] - TRUE_LAT) * mlat
        i = int(np.argmin(dxm ** 2 + dym ** 2))
        ax.plot([TRUE_LON, xy[i, 0]], [TRUE_LAT, xy[i, 1]], color=col, lw=1, ls=(0, (3, 2)))
        ev = RESULT["events"][p["event"] - 1]
        name = p["flight"].split("/")
        alt = p["event"] == 1
        ax.annotate(f"flyover {p['event']}: {name[2]}\nclosest {RESULT_D[p['event']]:,} m slant",
                    (xy[i, 0], xy[i, 1]), xytext=(14, -26) if alt else (-14, 22),
                    textcoords="offset points", ha="left" if alt else "right",
                    fontsize=9.5, color=th["ink"],
                    arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=0.7))
    truth_marker(ax, th)
    # 1 km scale bar
    x0 = extent[0] + 0.002
    y0 = extent[2] + 0.002
    ax.plot([x0, x0 + 1000 / mlon], [y0, y0], color=th["ink"], lw=2)
    ax.text(x0 + 500 / mlon, y0 + 0.0008, "1 km", ha="center", fontsize=9, color=th["ink"])
    ax.set_title("The flights it matched")
    ax.text(0.0, -0.02, "dashed: horizontal offset at closest approach; the rest of the\nslant distance is the plane's height", transform=ax.transAxes,
            va="top", fontsize=9.5, color=th["ink2"])

    # Right: the dense 25 m score map around the winner.
    ax = axes[1]
    pad = 0.0135
    extent = (TRUE_LON - pad / 0.76, TRUE_LON + pad / 0.76, TRUE_LAT - pad, TRUE_LAT + pad)
    map_axes(ax, th, extent)
    cmap = LinearSegmentedColormap.from_list("seq", th["ramp"])
    sc = ax.scatter(fx, fy, c=fv, cmap=cmap, s=11, marker="s", lw=0)
    ax.scatter([BEST["lon"]], [BEST["lat"]], marker="o", s=70, facecolor="none",
               edgecolor=th["ink"], linewidth=1.6, zorder=9)
    truth_marker(ax, th)
    d = haversine(BEST["lat"], BEST["lon"], TRUE_LAT, TRUE_LON)
    ax.annotate(f"best score ○ is {d:.0f} m from the true location ★",
                (BEST["lon"], BEST["lat"]), xytext=(0.03, 0.95), textcoords="axes fraction",
                fontsize=9.5, color=th["ink"], va="top")
    x0 = extent[0] + 0.0015
    y0 = extent[2] + 0.0015
    ax.plot([x0, x0 + 500 / mlon], [y0, y0], color=th["ink"], lw=2)
    ax.text(x0 + 250 / mlon, y0 + 0.0005, "500 m", ha="center", fontsize=9, color=th["ink"])
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.outline.set_visible(False)
    cb.set_label("Doppler score (both flyovers)", color=th["ink2"])
    ax.set_title("Refined score at 25 m")
    ax.text(0.0, -0.02, "the score falls off within a few hundred metres", transform=ax.transAxes,
            va="top", fontsize=9.5, color=th["ink2"])
    fig.suptitle("5  Rank and refine: the best-scoring flight pair pins down the spot",
                 x=0.0, ha="left", fontsize=12.5, fontweight="bold", color=th["ink"])
    fig.tight_layout()
    save(fig, "05-result", theme)


for theme, th in THEMES.items():
    style(th)
    fig_spectrogram(th, theme)
    fig_model(th, theme)
    fig_measure(th, theme)
    fig_funnel(th, theme)
    fig_result(th, theme)
print("wrote", sorted(p.name for p in DST.glob("*.png")))

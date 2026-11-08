#!/usr/bin/env bash
# Download one day of global historical ADS-B (adsb.lol, ODbL/CC0, readsb trace format) and build a
# audio2loc cache restricted to a region / altitude.
#
#   tools/fetch_adsblol.sh 2025-09-08 data/ [W,S,E,N]
#
# A day is ~3-4 GB of tar parts; only traces/ is extracted (files stay gzip-compressed).
set -euo pipefail
DATE=${1:?usage: fetch_adsblol.sh YYYY-MM-DD OUTDIR [W,S,E,N]}
OUT=${2:?outdir}
BBOX=${3:--125,24,-66,50}   # continental US
YEAR=${DATE:0:4}
TAG="v${DATE//-/.}-planes-readsb-prod-0"
REPO="adsblol/globe_history_${YEAR}"
mkdir -p "$OUT"
echo "release $REPO@$TAG"
URLS=$(curl -fsSL "https://api.github.com/repos/$REPO/releases/tags/$TAG" |
  python3 -c 'import json,sys; [print(a["browser_download_url"]) for a in sorted(json.load(sys.stdin)["assets"], key=lambda a:a["name"]) if ".tar." in a["name"]]')
[ -n "$URLS" ] || { echo "no tar parts found for $TAG" >&2; exit 1; }
echo "$URLS"
# parts are a split tar: concatenate the streams and extract traces/ on the fly
( for u in $URLS; do curl -fL --retry 3 "$u"; done ) | tar -x -C "$OUT" --wildcards 'traces/*' || true
BIN=${AUDIO2LOC:-$(dirname "$0")/../build/audio2loc}
"$BIN" ingest --adsb "$OUT/traces" --bbox "$BBOX" --max-alt 7000 -o "$OUT/adsb_${DATE}.plb"
echo "cache: $OUT/adsb_${DATE}.plb  (pass to: audio2loc locate --adsb $OUT/adsb_${DATE}.plb ...)"

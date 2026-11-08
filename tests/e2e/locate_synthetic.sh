#!/usr/bin/env bash
# End to end: simulate a scenario with a hidden observer, run the suggested
# `locate`, and require the answer within MAX_M metres of the truth.
#
#   locate_synthetic.sh <audio2loc> <work dir> <seed> <max metres>
set -euo pipefail
BIN=$1 WORK=$2 SEED=$3 MAX_M=$4
rm -rf "$WORK" && mkdir -p "$WORK" && cd "$WORK"
"$BIN" simulate --prefix sim --seed "$SEED" -q
read -r -a cmd <<<"$(grep '^suggested:' sim_truth.txt | cut -d' ' -f3-)"
"$BIN" "${cmd[@]}" -q --out out | tee locate.txt
awk -v max="$MAX_M" '
  /^observer_lat/ { la = $2 } /^observer_lon/ { lo = $2 }
  /^BEST:/ { gsub(",", "", $2); bla = $2; blo = $3 }
  END {
    if (bla == "") { print "no BEST line"; exit 1 }
    r = 6371000; d2r = atan2(0, -1) / 180
    x = (blo - lo) * d2r * cos(la * d2r); y = (bla - la) * d2r
    d = r * sqrt(x * x + y * y)
    printf "located %.1f m from the hidden observer (limit %g m)\n", d, max
    exit d <= max ? 0 : 1
  }' sim_truth.txt locate.txt
test -s out/audio2loc.qgs && test -s out/result.json

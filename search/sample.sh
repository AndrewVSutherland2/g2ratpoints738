#!/bin/bash
# sample.sh OUTDIR [JOBS=16] [P=2] -- how much the selection loses: a sample of 87 models of data/models.txt (all seven
# with one node, 20 with two nodes and 60 with three, at regular intervals in the file) is searched again with the score
# bound lowered from 40 to 30 and the thresholds of the first two point counts lowered from 250 and 330 to 200 and 300.
# OUTDIR/search_sample.txt lists every curve closed, and says whether the search with the standard settings has it
# (data/search_curves.txt) or which cut removed it.  About fifteen minutes on a 16-core machine.
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: sample.sh OUTDIR [JOBS=16] [P=2]" >&2; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE"); export HERE MODELS=${MODELS:-$ROOT/data/models.txt}; mkdir -p "$1"; export O=$(cd "$1" && pwd) P=${3:-2}; JOBS=${2:-16}
grep -v '^#' "$MODELS" | awk -F' [|] ' '{ n[$2]++; name[$2, n[$2]] = $1 } END { for (i = 1; i <= n[1]; i++) print name[1, i]; for (i = 0; i < 20; i++) print name[2, 1 + int(i * n[2] / 20)]; for (i = 0; i < 60; i++) print name[3, 1 + int(i * n[3] / 60)] }' > "$O/models.txt"
rm -f "$O/failed.txt"
xargs -P "$JOBS" -I{} bash -c 'SMIN=30 T1=200 T2=300 "$HERE/search_model.sh" {} "$O/{}" > "$O/{}.log" 2>&1 || echo "FAILED {}" >> "$O/failed.txt"' < "$O/models.txt"
[ ! -s "$O/failed.txt" ] || { echo "sample.sh: $(wc -l < "$O/failed.txt") models failed" >&2; exit 1; }
python3 "$HERE/sample_table.py" "$O" "$ROOT/data" > "$O/search_sample.txt"
tail -n +1 "$O/sample_summary.txt"

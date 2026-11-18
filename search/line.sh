#!/bin/bash
# line.sh NAME a b c -- the curve cut out by the line a z0 + b z1 + c z2 = 0 on the model NAME of data/models.txt:
# its reduced sextic, the number of rational points of height <= 10^5, and the number after the four-point closure.
# (The same programs as search/search_model.sh, for one line of any height.)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE"); MODELS=${MODELS:-$ROOT/data/models.txt}; SRC=${SRC:-$ROOT/src}
[ $# -eq 4 ] || { echo "usage: line.sh NAME a b c" >&2; exit 2; }
ROW=$(awk -F' [|] ' -v m="$1" '$1 == m' "$MODELS"); [ -n "$ROW" ] || { echo "line.sh: no model $1 in $MODELS" >&2; exit 2; }
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
echo "$ROW" | awk -F' [|] ' '{print $3}' > "$W/sextic.txt"
for b in linesearch g2rp cubic; do [ -x "$SRC/$b" ] || { echo "line.sh: program $b not found in $SRC; run make in $ROOT" >&2; exit 2; }; done
echo "$2 $3 $4" | "$SRC/linesearch" "$W/sextic.txt" 1 -exactbatch > "$W/raw" || { echo "line.sh: linesearch failed" >&2; exit 1; }
RED_IN="$W/raw" RED_OUT="$W/red" gp -q "$HERE/reduce.gp" < /dev/null > "$W/log" 2>&1
[ "$(tail -1 "$W/red")" = "DONE 1" ] || { echo "line.sh: the line does not cut out a curve of genus 2" >&2; exit 1; }
F=$(head -1 "$W/red" | cut -d' ' -f4-)
"$SRC/g2rp" "$F" 100000 -t "${P:-2}" > "$W/p" 2> /dev/null || { echo "line.sh: g2rp failed" >&2; exit 1; }
"$SRC/cubic" "$F" "$W/p" -t "${P:-2}" > "$W/x" 2> /dev/null || { echo "line.sh: cubic failed" >&2; exit 1; }
NS=$(grep -c . "$W/p" || true); NX=$(grep -c . "$W/x" || true)
echo "model $1, line ($2, $3, $4): y^2 = g0 + g1 x + .. + g6 x^6 with g = $F"
echo "$NS points of height <= 10^5; $((NS + 2 * NX)) points after closure"

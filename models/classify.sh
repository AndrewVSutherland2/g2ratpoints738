#!/bin/bash
# classify.sh B [PAR] : canonical forms of all polarizations of lat/dp_<B>.txt -> lat/canon_<B>.txt
# (rows "id n2 index nroots nlines hash leaves"), in PAR parallel pieces.
set -euo pipefail
B=$1; P=${2:-16}; D=$(cd "$(dirname "$0")" && pwd); W=$(mktemp -d "$D/lat/cl_${B}_XXXXXX")
split -n l/"$P" -d -a 3 "$D/lat/dp_$B.txt" "$W/p."
for f in "$W"/p.???; do
  ( printf 'LAT="%s";\nDP="%s";\n' "$D/lat/lat_$B.txt" "$f" > "$f.pre.gp"
    gp -q "$f.pre.gp" "$D/config.gp" < /dev/null 2> "$f.err" | "$D/canon" "$D/lat/lat_$B.txt" > "$f.out" ) &
done; wait
cat "$W"/p.???.out > "$D/lat/canon_$B.txt"
[ "$(wc -l < "$D/lat/canon_$B.txt")" = "$(wc -l < "$D/lat/dp_$B.txt")" ] || { echo "classify.sh: $B: $(wc -l < "$D/lat/canon_$B.txt") canonical forms for $(wc -l < "$D/lat/dp_$B.txt") polarizations" >&2; grep -v Warning "$W"/p.*.err | head -5 >&2; exit 1; }
echo "$B: $(wc -l < "$D/lat/canon_$B.txt") polarizations; classes with 1, 2, 3 nodes: $(for k in 1 2 3; do awk -v k=$k '$2 == k {print $6}' "$D/lat/canon_$B.txt" | sort -u | wc -l; done | tr '\n' ' '); index != 1: $(awk '$3 != 1' "$D/lat/canon_$B.txt" | wc -l)"

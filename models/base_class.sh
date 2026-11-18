#!/bin/bash
# base_class.sh B : canonical form of the polarization D_B of the base B itself (row "B_0 ..." as in canon_<B>.txt)
set -euo pipefail
B=$1; D=$(cd "$(dirname "$0")" && pwd)
awk -v b="$B" '$1 == "D" { $1 = ""; print b "_0 1 |" $0 " | 0,0,0,0,0,0" }' "$D/lat/lat_$B.txt" > "$D/lat/self_$B.txt"
printf 'LAT="%s";\nDP="%s";\n' "$D/lat/lat_$B.txt" "$D/lat/self_$B.txt" > "$D/lat/self_$B.pre.gp"
gp -q "$D/lat/self_$B.pre.gp" "$D/config.gp" < /dev/null 2> /dev/null | "$D/canon" "$D/lat/lat_$B.txt"

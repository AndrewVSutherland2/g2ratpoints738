#!/bin/bash
# build_models.sh WORK [T] : from Elkies' sextic to the list of double-plane models of X_163 with one, two or three nodes
# that are reached by the conic-pair construction.  The only mathematical input is data/elkies163.txt (the 28
# coefficients of Elkies' sextic).  T = number of threads (default: all).  Nothing is written outside WORK: the programs
# of this directory are copied to WORK/root/models and run there.
#   1. the lines on which the sextic is a square (exhaustive search to height 150)
#   2. Elkies' model as a base: its conics, the Neron-Severi lattice, the classes of the lines and conics (both tables are
#      proved complete there), all polarizations C_a + C_b with at most three orthogonal roots, their canonical forms
#   3. the one-node classes other than Elkies' own: built as models, moved to coordinates with the node at (0:0:1)
#   4. the same as 2. on each of them
#   5. one conic pair for every class (the classes with one node must be the bases themselves)
#   6. the sextic of every class with two or three nodes: conic-pair construction, reduction, MinRedTernaryForm
#   7. the list of models (WORK/models.txt, the format of data/models.txt), with the number of nodes and the number of
#      lines of every model checked against the lattice; WORK/minimisation.txt (what the minimisation changed);
#      WORK/elkies163_lines.txt and WORK/elkies163_conics.txt (the lines and conics of Elkies' sextic)
# Times of the stages are appended to WORK/times.txt.
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: build_models.sh WORK [threads]" >&2; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd); TOP=$(dirname "$HERE"); mkdir -p "$1"; W=$(cd "$1" && pwd); T=${2:-$(nproc)}
for b in models/canon models/frames models/sqfilter models/conicmodp src/linesearch; do [ -x "$TOP/$b" ] || { echo "build_models.sh: program $b not built; run make in $TOP" >&2; exit 2; }; done
R=$W/root; [ -e "$R" ] && { echo "build_models.sh: $R exists" >&2; exit 2; }
mkdir -p "$R/data" "$R/src" "$W/tmp"; cp -r "$HERE" "$R/models"; cp "$TOP/src/linesearch" "$R/src/"; cp "$TOP/data/elkies163.txt" "$R/data/"
export TMPDIR=$W/tmp; FS=$R/models; CP=$R/models; cd "$R"; mkdir -p "$CP/out" "$CP/logs" "$FS/lat" "$FS/min"
stamp() { echo "$1 $(date +%s.%N)" >> "$W/times.txt"; times > "$W/cpu.$1"; echo "$(date +%T) $1"; }      # wall time and processor time (this shell and its finished children)
build() {      # build JOBS PAR : rows "NAME i j BASE ..." (BASE = E for Elkies' model)
  xargs -P "$2" -L 1 bash -c 'B=$3; [ "$B" = E ] && B=""; cd "$CP" && BASE="$B" FRAMES="$CP/frames" ./pipeline.sh "$0" "$1" "$2" 0 1 1 < /dev/null > "logs/pipeline_$0.log" 2>&1 || echo "$0" >> "$W/failed_builds.txt"' < "$1"; }
export CP W
: > "$W/times.txt"; stamp start
# 1.
"$CP/sqlines.sh" data/elkies163.txt 150 "$W/elkies_sq150.txt" "$T" 2> "$W/sqlines.log"; grep -v '^#' "$W/elkies_sq150.txt" > data/elkies163_lines.txt
stamp lines
# 2.
"$FS/base.sh" E "$T" | tee -a "$W/bases.log"
stamp base_E
# 3.
python3 "$FS/choose_pairs.py" "$FS/lat/sel_E.txt" E > /dev/null
H0=$(awk '{print $6}' "$FS/lat/selfhash_E.txt")
awk -v h="$H0" '$1 == 1 && $2 != h { k++; print "b" k, $4, $5, "E" }' "$FS/lat/sel_E.txt" > "$FS/lat/jobs_bases.txt"
NB=$(wc -l < "$FS/lat/jobs_bases.txt"); build "$FS/lat/jobs_bases.txt" "$NB"
[ ! -s "$W/failed_builds.txt" ] || { echo "build_models.sh: a base failed to build: $(tr '\n' ' ' < "$W/failed_builds.txt")" >&2; exit 1; }
BASES=""; for k in $(seq 1 "$NB"); do printf 'SRC="b%d";\nDST="B%d";\n' "$k" "$k" > "$W/tmp/nf_$k.gp"; gp -q "$W/tmp/nf_$k.gp" "$FS/node_frame.gp" < /dev/null 2>&1 | grep -v Warning | tee -a "$W/bases.log"; [ -s "data/x163_B$k.txt" ] || { echo "build_models.sh: no base B$k" >&2; exit 1; }; BASES="$BASES B$k"; done
stamp bases_built
# 4.
for B in $BASES; do "$FS/base.sh" "$B" $(( (T + NB - 1) / NB )) >> "$W/bases.log.$B" 2>&1 & done; wait
for B in $BASES; do cat "$W/bases.log.$B" | tee -a "$W/bases.log"; [ -s "$FS/lat/selfhash_$B.txt" ] || { echo "build_models.sh: base $B failed" >&2; exit 1; }; done
stamp bases_classified
# 5.
python3 "$FS/choose_pairs.py" "$FS/lat/sel_all.txt" E $BASES | tee -a "$W/bases.log"
diff <(awk '$1 == 1 {print $2}' "$FS/lat/sel_all.txt" | sort) <(cat "$FS"/lat/selfhash_*.txt | awk '{print $6}' | sort) > /dev/null || { echo "build_models.sh: the one-node classes are not the bases" >&2; exit 1; }
awk '$1 >= 2 { k[$1]++; printf "%s%0*d %d %d %s %d %d %s %d\n", ($1 == 2 ? "p" : "t"), ($1 == 2 ? 3 : 4), k[$1], $4, $5, $3, $8, $9, $2, $11 }' "$FS/lat/sel_all.txt" > "$FS/lat/jobs_all.txt"
stamp selected
# 6.
build "$FS/lat/jobs_all.txt" "$T"
[ ! -s "$W/failed_builds.txt" ] || { echo "build_models.sh: $(wc -l < "$W/failed_builds.txt") models failed to build: $(head -5 "$W/failed_builds.txt" | tr '\n' ' ')" >&2; exit 1; }
stamp models_built
# 7.
"$CP/minred_file.sh" E "$R/data/elkies163.txt" "$R/data/elkies163_lines.txt" "0,0,1" "$FS/min" > "$W/minred_E.log"
( echo "E $(awk '{print $5}' "$FS/lat/selfhash_E.txt") 1"; for k in $(seq 1 "$NB"); do echo "b$k $(awk '{print $5}' "$FS/lat/selfhash_B$k.txt") 1"; done; awk '{print $1, $8, ($1 ~ /^p/ ? 2 : 3)}' "$FS/lat/jobs_all.txt" ) > "$W/expect.txt"
python3 "$FS/write_models.py" "$W/models.txt" "E=$FS/min" -l <(awk 'NR > 1 {print $1}' "$W/expect.txt") "$R/data" | tee -a "$W/bases.log"
awk -F' [|] ' 'NR == FNR { split($0, a, " "); nl[a[1]] = a[2]; nn[a[1]] = a[3]; next } /^#/ { next } { n = split($5, l, ";"); m = split($4, s, ";"); if (2 * n != nl[$1] || m != nn[$1]) { print "build_models.sh: " $1 ": " n " lines, " m " singular points; the lattice gives " nl[$1] / 2 " and " nn[$1]; bad = 1 } } END { exit bad }' "$W/expect.txt" "$W/models.txt" || { echo "build_models.sh: a model does not have the lines or nodes of its class" >&2; exit 1; }
echo "all $(grep -vc '^#' "$W/models.txt") models have the number of nodes and of lines of their class" | tee -a "$W/bases.log"
# the lines and conics of Elkies' sextic, and what the minimisation did: name, largest coefficient before and after,
# determinant of the change of coordinates (rows in the order of models.txt)
cp data/elkies163_lines.txt data/elkies163_conics.txt "$W/"
{ echo "# Minimisation and reduction (MinRedTernaryForm) of the sextics of models.txt: name, largest coefficient before and after,"
  echo "# determinant of the change of coordinates y = N x (1 or -1: the sextic was minimal already).  The first row is Elkies' sextic."
  tail -1 "$W/minred_E.log" | awk '{ print $1, $2, $3, $4 }'
  awk 'NR > 1 { print $1 }' "$W/expect.txt" | while read -r n; do sed -n "s/^minred.m: max |coeff| \([0-9]*\) -> \([0-9]*\), det N = \([-0-9\/]*\),.*/$n \1 \2 \3/p" "$CP/logs/minredm_$n.log"; done
} > "$W/minimisation.txt"
[ "$(grep -vc '^#' "$W/minimisation.txt")" = "$(grep -vc '^#' "$W/models.txt")" ] || { echo "build_models.sh: minimisation.txt is incomplete" >&2; exit 1; }
stamp done
while read -r name t; do echo "$name $t $(awk '{ for (i = 1; i <= NF; i++) { split($i, a, "m"); sub(/s/, "", a[2]); s += a[1] * 60 + a[2] } } END { printf "%.1f", s }' "$W/cpu.$name")"; done < "$W/times.txt" > "$W/times_cpu.txt"
awk 'NR > 1 { printf "%-18s %8.1f s wall %9.1f s processor\n", $1, $2 - p, $3 - c } { p = $2; c = $3 } NR == 1 { s = $2; d = $3 } END { printf "%-18s %8.1f s wall %9.1f s processor\n", "total", p - s, c - d }' "$W/times_cpu.txt" | tee "$W/times_summary.txt"

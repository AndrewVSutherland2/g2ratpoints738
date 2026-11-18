#!/bin/bash
# search_model.sh NAME [OUTDIR] -- the search on one model of data/models.txt.
#
#   1. the points of the plane of height <= HP at which S is a rational square ("surface points"), off the lines on
#      which S is a square                                                                      (k3pts)
#   2. the number of these points on every line of height <= L, the score of the line          (k3inc)
#   3. lines through a singular point are dropped; the lines with score >= SMIN are kept
#   4. for each line kept: the sextic it cuts out (linesearch -exactbatch) and a reduced model  (PARI/GP, hyperellred)
#   5. rational points of height <= H1, 10^4, 10^5 on these curves, keeping those with at least T1, T2, T3 (g2rp)
#   6. for the curves with at least CMIN points of height <= 10^5: the closure of these points under the relation
#      given by a cubic through four points                                                     (cubic)
#
# Settings (environment; the defaults are those of the note): L=400 HP=4000 SMIN=40 H1=2000 T1=250 T2=330 T3=400 CMIN=440.
# With SMIN=0 the NTOP lines of highest score are kept instead, with ties (NTOP, or NTOP1/NTOP2/NTOP3 by number of nodes).
# P (default 2) is the number of threads.  MODELS (default data/models.txt) is the list of models, SRC (default src/) the
# directory of the programs.  HIST=1 also writes score_hist.txt (the number of lines, not through a singular point, with
# each score from KMIN = 15 up).
# Results in OUTDIR (default ./out_NAME): summary.txt; scores.txt ("a b c score" for the lines kept, by score);
# counts_first.txt, counts_1e4.txt, counts_1e5.txt ("points a b c" for the curves counted to H1, 10^4, 10^5);
# candidates.txt ("points of height <= 10^5, line, reduced sextic g0..g6" for the curves with at least T3 points);
# closed.txt ("points after closure, points of height <= 10^5, line, sextic").
# Needs the programs k3pts, k3inc, linesearch, g2rp and cubic (make, in the top directory), gp (PARI/GP 2.17 or later
# for the reduced equations of the stored run: see README.md) and the usual shell tools of a Linux system.  A run that
# fails writes no summary.txt and has exit status 1.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE")
MODELS=${MODELS:-$ROOT/data/models.txt}
SRC=${SRC:-$ROOT/src}
[ $# -ge 1 ] || { echo "usage: search_model.sh NAME [OUTDIR]" >&2; exit 2; }
M=$1; O=${2:-out_$M}; L=${L:-400}; HP=${HP:-4000}; KMIN=${KMIN:-15}; H1=${H1:-2000}; T1=${T1:-250}; T2=${T2:-330}; T3=${T3:-400}; CMIN=${CMIN:-440}; P=${P:-2}; SMIN=${SMIN:-40}
for b in k3pts k3inc linesearch g2rp cubic; do [ -x "${SRC:-}/$b" ] || { echo "search_model.sh: program $b not found in $SRC; run make in $ROOT" >&2; exit 2; }; done
ROW=$(awk -F' [|] ' -v m="$M" '$1 == m' "$MODELS"); [ -n "$ROW" ] || { echo "search_model.sh: no model $M in $MODELS" >&2; exit 2; }
mkdir -p "$O"; O=$(cd "$O" && pwd); W=$(mktemp -d "$O/work.XXXXXX"); trap 'rm -rf "$W"' EXIT
rm -f "$O/summary.txt"            # the summary is written last: a run that fails leaves none
K=$(echo "$ROW" | awk -F' [|] ' '{print $2}'); NODES=$(echo "$ROW" | awk -F' [|] ' '{print $4}')
echo "$ROW" | awk -F' [|] ' '{print $3}' > "$W/sextic.txt"
echo "$ROW" | awk -F' [|] ' '{n = split($5, a, ";"); for (i = 1; i <= n; i++) { gsub(",", " ", a[i]); print a[i] } }' > "$W/sqlines.txt"
case $K in 1) D=${NTOP1:-30000};; 2) D=${NTOP2:-12000};; *) D=${NTOP3:-6000};; esac; NTOP=${NTOP:-$D}
tm() { date +%s.%N; }
# processor time (user + system, this shell and the children it has waited for) between two of the time marks
cpud() { awk '{ for (i = 1; i <= NF; i++) { split($i, a, "m"); sub(/s/, "", a[2]); s[FILENAME] += a[1] * 60 + a[2] } } END { printf "%.1f", s[ARGV[2]] - s[ARGV[1]] }' "$W/cpu.$1" "$W/cpu.$2"; }
t0=$(tm); times > "$W/cpu.0"
# 1, 2: surface points and their numbers on the lines of height <= L
if [ -n "${PTS:-}" ]; then awk -v h="$HP" 'function a(x) { return x < 0 ? -x : x } a($1) <= h && a($2) <= h && a($3) <= h' "$PTS" > "$W/pts.txt"; else "$SRC/k3pts" "$W/sextic.txt" "$HP" -t "$P" -skip "$W/sqlines.txt" > "$W/pts.txt" 2> "$W/pts.err"; fi; t1=$(tm); times > "$W/cpu.1"
if [ "${SLABS:-1}" -gt 1 ]; then : > "$W/inc.raw"
  for ((i = 0; i < SLABS; i++)); do lo=$((i * L / SLABS)); hi=$(( i == SLABS - 1 ? L : (i + 1) * L / SLABS - 1 ))
    "$SRC/k3inc" "$W/pts.txt" "$L" "$HP" "$KMIN" "$KMIN" -t "$P" -slots "${SLOTS:-28}" -A "$lo" "$hi" >> "$W/inc.raw" 2>> "$O/inc.err" || { echo "k3inc failed on slab $lo..$hi" >&2; exit 1; }; done
  sort -k1,1nr -k2,2nr -k3,3n -k4,4n -k5,5n "$W/inc.raw" > "$W/inc.txt"; rm -f "$W/inc.raw"
else "$SRC/k3inc" "$W/pts.txt" "$L" "$HP" "$KMIN" "$KMIN" -t "$P" -slots "${SLOTS:-28}" > "$W/inc.txt" 2> "$O/inc.err" || { echo "k3inc failed" >&2; exit 1; }; fi; t2=$(tm); times > "$W/cpu.2"
# 3: selection.  inc.txt is sorted by score: "score score A B C"
awk -v nodes="$NODES" -v ntop="$NTOP" -v lmin="${LMIN:-0}" -v smin="$SMIN" 'BEGIN { n = split(nodes, a, ";"); for (i = 1; i <= n; i++) { split(a[i], v, ","); x[i] = v[1]; y[i] = v[2]; z[i] = v[3] } }
  { h = ($3 < 0 ? -$3 : $3); t = ($4 < 0 ? -$4 : $4); if (t > h) h = t; t = ($5 < 0 ? -$5 : $5); if (t > h) h = t; if (h <= lmin) next;
    for (i = 1; i <= n; i++) if ($3*x[i] + $4*y[i] + $5*z[i] == 0) next; k++; if (smin > 0) { if ($1 < smin) exit } else if (k > ntop && $1 < last) exit; last = $1; print $3, $4, $5, $1 }' "$W/inc.txt" > "$O/scores.txt"
[ -z "${HIST:-}" ] || awk -v nodes="$NODES" 'BEGIN { n = split(nodes, a, ";"); for (i = 1; i <= n; i++) { split(a[i], v, ","); x[i] = v[1]; y[i] = v[2]; z[i] = v[3] } }
  { for (i = 1; i <= n; i++) if ($3*x[i] + $4*y[i] + $5*z[i] == 0) next; h[$1]++ } END { for (s in h) print s, h[s] }' "$W/inc.txt" | sort -n > "$O/score_hist.txt"
SBOUND=$SMIN; [ "$SBOUND" -gt 0 ] && SEL="lines with score >= $SBOUND" || SEL="NTOP=$NTOP"
NP=$(wc -l < "$W/pts.txt"); NI=$(wc -l < "$W/inc.txt"); NK=$(wc -l < "$O/scores.txt"); SMIN=$(tail -1 "$O/scores.txt" | cut -d' ' -f4)
cut -d' ' -f1-3 "$O/scores.txt" > "$W/lines.txt"
# 4: restriction and reduction, in P parallel pieces
split -n l/"$P" -d -a 3 "$W/lines.txt" "$W/l."
for f in "$W"/l.*; do ( "$SRC/linesearch" "$W/sextic.txt" 1 -exactbatch < "$f" > "$f.raw" && RED_IN="$f.raw" RED_OUT="$f.red" gp -q "$HERE/reduce.gp" < /dev/null > "$f.log" 2>&1 ) & done; wait      # a piece that failed has no DONE row (checked next)
for f in "$W"/l.???; do tail -1 "$f.red" 2>/dev/null | grep -q '^DONE ' || { echo "search_model.sh: the reduction failed on $f: $(tail -2 "$f.log" | tr '\n' ' ')" >&2; trap - EXIT; exit 1; }; done
cat "$W"/l.???.red | grep -v '^DONE ' > "$W/curves.txt"; t3=$(tm); times > "$W/cpu.3"
[ -z "${KEEPCURVES:-}" ] || cp "$W/curves.txt" "$O/curves.txt"
# 5: point counts.  g2rp -batch H reads "A B C g0..g6" and prints "count A B C g0..g6"
count() { local IN=$1 H=$2 OUT=$3 OPT=${4:-}; : > "$OUT"; [ -s "$IN" ] || return 0
  rm -f "$W"/c.*; split -n l/$((P * 8 < $(wc -l < "$IN") ? P * 8 : $(wc -l < "$IN"))) -d -a 4 "$IN" "$W/c."
  ls "$W" | grep '^c\.[0-9]*$' | xargs -P "$P" -I{} bash -c '"$0/g2rp" -batch "$1" -t 1 -q $3 < "$2/{}" > "$2/{}.out" 2> /dev/null || exit 255' "$SRC" "$H" "$W" "$OPT" \
    || { echo "search_model.sh: g2rp failed in the point count to $H" >&2; exit 1; }
  cat "$W"/c.*.out | sort -k1,1nr -k2,2n -k3,3n -k4,4n > "$OUT"
  [ "$(wc -l < "$OUT")" = "$(wc -l < "$IN")" ] || { echo "search_model.sh: point count to $H failed" >&2; exit 1; }; }
keep() { awk -v t="$2" '$1 >= t { $1 = ""; sub(/^ /, ""); print }' "$1"; }
count "$W/curves.txt" "$H1" "$O/counts_first.txt" "${G1OPT--Pmax 130 -Mmax 800 -nvec 14}";   keep "$O/counts_first.txt" "$T1" > "$W/k1.txt"; t3a=$(tm); times > "$W/cpu.3a"
count "$W/k1.txt" 10000 "$O/counts_1e4.txt";      keep "$O/counts_1e4.txt" "$T2" > "$W/k2.txt"; t3b=$(tm); times > "$W/cpu.3b"
count "$W/k2.txt" 100000 "$O/counts_1e5.txt";     awk -v t="$T3" '$1 >= t' "$O/counts_1e5.txt" > "$O/candidates.txt"; t4=$(tm); times > "$W/cpu.4"
# 6: closure
: > "$W/closed.txt"
close_one() { local n=$1 A=$2 B=$3 C=$4; shift 4; local f="$*" tag="${A}_${B}_${C}" ns nx
  "$SRC/g2rp" "$f" 100000 -t 1 > "$W/p.$tag" 2> /dev/null || { echo "search_model.sh: g2rp failed on the line $A $B $C" >&2; return 1; }
  "$SRC/cubic" "$f" "$W/p.$tag" -t "$TC" > "$W/x.$tag" 2> /dev/null || { echo "search_model.sh: cubic failed on the line $A $B $C" >&2; return 1; }
  ns=$(grep -c . "$W/p.$tag" || true); nx=$(grep -c . "$W/x.$tag" || true)
  [ "$ns" = "$n" ] || { echo "search_model.sh: $ns points of height <= 1e5 on the line $A $B $C, counted $n before" >&2; return 1; }
  echo "$((ns + 2 * nx)) $ns $A $B $C $f"; }
NC=$(awk -v c="$CMIN" '$1 >= c' "$O/candidates.txt" | wc -l); TC=$(( NC > 0 && NC < P ? P / NC : 1 ))       # threads for each closure
export -f close_one; export SRC W TC
# a closure that fails stops the run (exit status 255 makes xargs stop), and no summary is written
awk -v c="$CMIN" '$1 >= c' "$O/candidates.txt" | xargs -r -P "$P" -L 1 bash -c 'close_one "$@" || exit 255' _ > "$W/closed.txt" || { echo "search_model.sh: a closure failed" >&2; exit 1; }
[ "$(wc -l < "$W/closed.txt")" = "$NC" ] || { echo "search_model.sh: $NC curves to close, $(wc -l < "$W/closed.txt") closed" >&2; exit 1; }
sort -k1,1nr -k2,2nr "$W/closed.txt" > "$O/closed.txt"; t5=$(tm); times > "$W/cpu.5"
cut -d' ' -f1-4 "$O/counts_first.txt" > "$W/c3"; mv "$W/c3" "$O/counts_first.txt"       # "count A B C" is enough for the two long lists
cut -d' ' -f1-4 "$O/counts_1e4.txt" > "$W/c4"; mv "$W/c4" "$O/counts_1e4.txt"
{ echo "model $M: $K nodes, $(wc -l < "$W/sqlines.txt") lines on which S is a square; L=$L HP=$HP $SEL first count to $H1, cuts $T1/$T2/$T3 closure from $CMIN; $P threads"
  echo "surface points of height <= $HP off those lines: $NP"
  echo "lines of height <= $L with at least $KMIN of them: $NI"
  echo "lines kept: $NK (score >= $SMIN); best score $(head -1 "$O/scores.txt" | awk '{print $4 " on the line (" $1 ", " $2 ", " $3 ")"}')"
  echo "curves counted to $H1: $(wc -l < "$O/counts_first.txt"); with at least $T1 points: $(wc -l < "$W/k1.txt")"
  echo "curves counted to 10^4: $(wc -l < "$O/counts_1e4.txt"); with at least $T2 points: $(wc -l < "$W/k2.txt")"
  echo "curves counted to 10^5: $(wc -l < "$O/counts_1e5.txt"); with at least $T3 points: $(wc -l < "$O/candidates.txt")"
  echo "curves closed: $(wc -l < "$O/closed.txt")$(head -1 "$O/closed.txt" | awk 'NF { print "; best " $1 " points (" $2 " of height <= 10^5) on the line (" $3 ", " $4 ", " $5 ")" }')"
  echo "seconds: points $(awk -v a="$t1" -v b="$t0" 'BEGIN { printf "%.1f", a - b }'), scores $(awk -v a="$t2" -v b="$t1" 'BEGIN { printf "%.1f", a - b }'), reduction $(awk -v a="$t3" -v b="$t2" 'BEGIN { printf "%.1f", a - b }'), counts $(awk -v a="$t4" -v b="$t3" 'BEGIN { printf "%.1f", a - b }'), closure $(awk -v a="$t5" -v b="$t4" 'BEGIN { printf "%.1f", a - b }'), total $(awk -v a="$t5" -v b="$t0" 'BEGIN { printf "%.1f", a - b }')"
  echo "seconds in the counts: first $(awk -v a="$t3a" -v b="$t3" 'BEGIN { printf "%.1f", a - b }'), to 10^4 $(awk -v a="$t3b" -v b="$t3a" 'BEGIN { printf "%.1f", a - b }'), to 10^5 $(awk -v a="$t4" -v b="$t3b" 'BEGIN { printf "%.1f", a - b }')"
  echo "processor seconds: points $(cpud 0 1), scores $(cpud 1 2), reduction $(cpud 2 3), counts $(cpud 3 4), closure $(cpud 4 5), total $(cpud 0 5)"
  echo "processor seconds in the counts: first $(cpud 3 3a), to 10^4 $(cpud 3a 3b), to 10^5 $(cpud 3b 4)"
} > "$O/summary.txt"
cat "$O/summary.txt"

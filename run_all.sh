#!/bin/bash
# run_all.sh WORK [MODELS] : the whole computation, from Elkies' sextic (data/elkies163.txt) to the curves.
#   models/build_models.sh   the models with one, two and three nodes                       -> WORK/models.txt
#   search/search_all.sh     the search of all lines of height <= 400 on each of them       -> WORK/search/NAME
#   search/collect.py        the results as two tables, compared with data/                 -> WORK/data
#   search/time_by_part.py   where the time went                                            -> WORK/time_by_part.txt
# With a second argument (a list of models: data/models.txt) the models are not built and the search runs on that list;
# this needs neither Magma nor numpy.
# WORK must not contain an earlier run.  Nothing outside WORK is written.  The settings of the search are the defaults
# of search/search_model.sh (they can be changed through the environment); JOBS models are searched at a time, with P
# threads each (default: half the number of processors, and 2).  About an hour on a machine with 16 cores (32 threads)
# and AVX-512; WORK grows to about 1 GB.  MARCH (default native) is the -march option of the compiler.
set -euo pipefail
[ $# -eq 1 ] || [ $# -eq 2 ] || { echo "usage: run_all.sh WORK [MODELS]   (WORK: a directory for the run, which must not hold an earlier run; MODELS: search this list of models instead of building them)" >&2; exit 2; }
TOP=$(cd "$(dirname "$0")" && pwd); GIVEN=""
if [ $# -eq 2 ]; then [ -s "$2" ] || { echo "run_all.sh: no list of models $2" >&2; exit 2; }; GIVEN=$(readlink -f "$2"); fi
for b in gcc make gp python3 bc awk xargs; do command -v $b > /dev/null || { echo "run_all.sh: $b not found" >&2; exit 2; }; done
[ -n "$GIVEN" ] || command -v magma > /dev/null || { echo "run_all.sh: magma not found (it is needed to build the models; \"run_all.sh WORK data/models.txt\" runs the search on the stored models without it)" >&2; exit 2; }
# hyperellred, which gives the reduced equation of each curve, reduces a binary quadratic form exactly since PARI 2.17;
# before that it returned another equation for about 0.6% of the curves, and the points of bounded height are counted on
# it.  With PARI 2.15 the construction also gives other coordinates for 78 of the models.
GPV=$(gp --version 2>&1 | sed -n 's/.*Version \([0-9]*\)\.\([0-9]*\)\..*/\1 \2/p' | head -1); OLDGP=""
if [ -n "$GPV" ]; then read -r GMAJ GMIN <<< "$GPV"; [ "$GMAJ" -gt 2 ] || { [ "$GMAJ" -eq 2 ] && [ "$GMIN" -ge 17 ]; } || OLDGP="PARI/GP $GMAJ.$GMIN is older than 2.17: some models get other coordinates and some curves another reduced equation than in the stored run, so the results will differ from data/ (see README.md)"; fi
[ -z "$OLDGP" ] || echo "run_all.sh: warning: $OLDGP" >&2
mkdir -p "$1"; W=$(cd "$1" && pwd)
[ ! -e "$W/root" ] && [ ! -e "$W/search" ] && [ ! -e "$W/all.log" ] || { echo "run_all.sh: $W holds an earlier run" >&2; exit 2; }
J=${JOBS:-$(( $(nproc) / 2 ))}; [ "$J" -ge 1 ] || J=1
make -s -C "$TOP"
echo "START $(date +%T) $(date +%s.%N)" > "$W/all.log"
if [ -n "$GIVEN" ]; then cp "$GIVEN" "$W/models.txt"
else "$TOP/models/build_models.sh" "$W" "$(nproc)" >> "$W/all.log" 2>&1 || { echo "run_all.sh: the construction of the models failed; see $W/all.log" >&2; exit 1; }
fi
echo "MODELS $(date +%T) $(date +%s.%N)" >> "$W/all.log"
"$TOP/search/search_all.sh" "$W/models.txt" "$W/search" "$J" "${P:-2}" >> "$W/all.log" 2>&1 || { echo "run_all.sh: the search failed; see $W/all.log and $W/search/failed.txt" >&2; exit 1; }
echo "DONE $(date +%T) $(date +%s.%N)" >> "$W/all.log"
times > "$W/cpu.all"          # processor time of the whole run (this shell and its children)
awk '$1 == "START" { s = $3 } $1 == "MODELS" { m = $3 } $1 == "DONE" { printf "models %.0f s, search %.0f s, total %.0f s = %.1f minutes\n", m - s, $3 - m, $3 - s, ($3 - s) / 60 }' "$W/all.log" | tee -a "$W/all.log"
python3 "$TOP/search/time_by_part.py" "$W" > "$W/time_by_part.txt"; cat "$W/time_by_part.txt"
mkdir -p "$W/data"
{ echo "# The run of run_all.sh that produced the tables of this directory: machine, programs, settings and times."
  echo "date: $(awk '$1 == "START" { print strftime("%Y-%m-%d %H:%M", $3) }' "$W/all.log")"
  echo "machine: $(lscpu | sed -n 's/^Model name: *//p' | head -1), $(nproc) threads; $(uname -sr)"
  echo "programs: $(gcc --version | head -1), 512-bit operations: $(printf '#include "simd.h"\nSIMD_NAME\n' | gcc -march="${MARCH:-native}" -E -P -I "$TOP/src" -x c - | tail -1 | tr -d '"'); $(gp --version 2>&1 | head -1 | sed 's/^ *//')$([ -n "$GIVEN" ] || echo "; Magma V$(magma -b <<< 'a, b, c := GetVersion(); printf "%o.%o-%o", a, b, c; quit;' 2> /dev/null | tail -1)")"
  [ -z "$GIVEN" ] || echo "models: not built in this run; the search ran on $(grep -vc '^#' "$W/models.txt") models of a given list"
  echo "settings: $(head -1 "$(ls "$W"/search/*/summary.txt | head -1)" | sed 's/^[^;]*; //'); $J models at a time"
  cat "$W/time_by_part.txt"; } > "$W/data/search_summary.txt"
[ -z "$OLDGP" ] || echo "note: $OLDGP"
python3 "$TOP/search/collect.py" "$W" "$W/data" -compare "$TOP/data"

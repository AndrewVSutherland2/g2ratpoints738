#!/bin/bash
# run_all.sh WORK : the whole computation, from Elkies' sextic (data/elkies163.txt) to the curves.
#   models/build_models.sh   the models with one, two and three nodes                       -> WORK/models.txt
#   search/search_all.sh     the search of all lines of height <= 400 on each of them       -> WORK/search/NAME
#   search/collect.py        the results as two tables, compared with data/                 -> WORK/data
#   search/time_by_part.py   where the time went                                            -> WORK/time_by_part.txt
# WORK must not contain an earlier run.  Nothing outside WORK is written.  The settings of the search are the defaults
# of search/search_model.sh (they can be changed through the environment); JOBS models are searched at a time, with P
# threads each (default: half the number of processors, and 2).  About an hour on a machine with 16 cores (32 threads);
# WORK grows to about 1 GB.
set -euo pipefail
TOP=$(cd "$(dirname "$0")" && pwd); mkdir -p "$1"; W=$(cd "$1" && pwd)
make -s -C "$TOP"
echo "START $(date +%T) $(date +%s.%N)" > "$W/all.log"
"$TOP/models/build_models.sh" "$W" "$(nproc)" >> "$W/all.log" 2>&1 || { echo "run_all.sh: the construction of the models failed; see $W/all.log" >&2; exit 1; }
echo "MODELS $(date +%T) $(date +%s.%N)" >> "$W/all.log"
"$TOP/search/search_all.sh" "$W/models.txt" "$W/search" "${JOBS:-$(( $(nproc) / 2 ))}" "${P:-2}" >> "$W/all.log" 2>&1 || { echo "run_all.sh: the search failed; see $W/all.log and $W/search/failed.txt" >&2; exit 1; }
echo "DONE $(date +%T) $(date +%s.%N)" >> "$W/all.log"
times > "$W/cpu.all"          # processor time of the whole run (this shell and its children)
awk '$1 == "START" { s = $3 } $1 == "MODELS" { m = $3 } $1 == "DONE" { printf "models %.0f s, search %.0f s, total %.0f s = %.1f minutes\n", m - s, $3 - m, $3 - s, ($3 - s) / 60 }' "$W/all.log" | tee -a "$W/all.log"
python3 "$TOP/search/time_by_part.py" "$W" > "$W/time_by_part.txt"; cat "$W/time_by_part.txt"
mkdir -p "$W/data"
{ echo "# The run of run_all.sh that produced the tables of this directory: machine, programs, settings and times."
  echo "date: $(awk '$1 == "START" { print strftime("%Y-%m-%d %H:%M", $3) }' "$W/all.log")"
  echo "machine: $(lscpu | sed -n 's/^Model name: *//p' | head -1), $(nproc) threads; $(uname -sr)"
  echo "programs: $(gcc --version | head -1); $(gp --version 2>&1 | head -1 | sed 's/^ *//'); Magma V$(magma -b <<< 'a, b, c := GetVersion(); printf "%o.%o-%o", a, b, c; quit;' 2> /dev/null | tail -1)"
  echo "settings: $(head -1 "$(ls "$W"/search/*/summary.txt | head -1)" | sed 's/^[^;]*; //'); ${JOBS:-$(( $(nproc) / 2 ))} models at a time"
  cat "$W/time_by_part.txt"; } > "$W/data/search_summary.txt"
python3 "$TOP/search/collect.py" "$W" "$W/data" -compare "$TOP/data"

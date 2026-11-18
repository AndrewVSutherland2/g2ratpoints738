#!/bin/bash
# search_all.sh MODELS OUTDIR [JOBS=16] [P=2] : the search (search/search_model.sh) on every model of MODELS, JOBS models
# at a time with P threads each.  The settings are taken from the environment (see search_model.sh).  The output of a model
# is in OUTDIR/NAME; OUTDIR.log holds the times.  Exit status 1 if a model failed.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); export HERE MODELS=$(readlink -f "$1"); mkdir -p "$2"; export O=$(cd "$2" && pwd) P=${4:-2}; JOBS=${3:-16}
export SRC=${SRC:-$(dirname "$HERE")/src}
N=$(grep -vc '^#' "$MODELS"); rm -f "$O/failed.txt"
t0=$(date +%s.%N); echo "START $(date +%T) $t0" > "$O.log"
grep -v '^#' "$MODELS" | awk -F' [|] ' '{print $1}' | xargs -P "$JOBS" -I{} bash -c '"$HERE/search_model.sh" {} "$O/{}" > "$O/{}.log" 2>&1 || echo "FAILED {}" >> "$O/failed.txt"'
t1=$(date +%s.%N); NF=0; [ ! -s "$O/failed.txt" ] || NF=$(wc -l < "$O/failed.txt")
echo "$JOBS jobs x $P threads: $(ls "$O"/*/summary.txt 2>/dev/null | wc -l) of $N models in $(echo "$t1 - $t0" | bc | cut -d. -f1) s; failed: $NF" | tee -a "$O.log"
echo "DONE $(date +%T) $t1" >> "$O.log"
[ "$NF" = 0 ] && [ "$(ls "$O"/*/summary.txt 2>/dev/null | wc -l)" = "$N" ]

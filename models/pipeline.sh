#!/bin/bash
# pipeline.sh NAME i j sign [H=150] [WRITE=1] : explicit model of D' = C_i + sign C_j, from the sections to the data files:
#   model.gp (sections, branch sextic, images of the known curves) -> reduce.gp (coordinate reduction; with FRAMES=<path of
#   the program frames> in the environment its frame search is done by that program, with the same result) -> minred.m, minred.gp
#   (Magma MinRedTernaryForm: minimisation at every prime and reduction; skipped with MINRED=0) -> sing.m (Magma:
#   singular points) -> sqlines.sh (exhaustive square lines of height <= H) -> finalize.gp (checks, data/x163_NAME.txt,
#   data/x163_NAME_lines.txt).  Logs in logs/<step>_NAME.log; a one-line summary is appended to out/summary.txt.
# Every step must succeed (exit status 0 and its "written" line): any gp error (run.sh scans the log), any failed witness
# check (error() in model.gp / finalize.gp) or Magma failure aborts the pipeline with status 1 and no summary line.
set -uo pipefail
cd "$(dirname "$0")" || exit 1
NAME=$1; I=$2; J=$3; SG=$4; H=${5:-150}; WRITE=${6:-1}; BASE=${BASE:-}   # BASE="<M>": build on the x163_<M> base surface
export TMPDIR=${TMPDIR:-/tmp}
T0=$(date +%s)
# The lifts of the conics and lines of the base are a prerequisite, made once for each base.  Many pipelines of the same
# base run at the same time, so the test and the computation are done under a lock: without it a pipeline could read the
# file while another was writing it.
mkdir -p out logs
LIFTS="out/lifts.gp"; [ -n "$BASE" ] && LIFTS="out/lifts_$BASE.gp"
( flock 9; [ -f "$LIFTS" ] || ./run.sh lifts.gp "logs/lifts${BASE:+_$BASE}.log" "BASE=\"$BASE\"" > /dev/null ) 9> "$LIFTS.lock" || { echo "$NAME: lifts.gp FAILED (missing $LIFTS)"; exit 1; }
[ -f "$LIFTS" ] || { echo "$NAME: $LIFTS still missing after lifts.gp"; exit 1; }
./run.sh model.gp "logs/model_$NAME.log" "PAIR=[$I,$J,$SG]" "NAME=\"$NAME\"" "BASE=\"$BASE\"" > /dev/null || { echo "$NAME ($I,$J,$SG): model.gp FAILED"; grep '^  \*\*\*' "logs/model_$NAME.log" | grep -v Warning | head -3; exit 1; }
grep -q "^written out/model_${NAME}_raw.gp" "logs/model_$NAME.log" || { echo "$NAME ($I,$J,$SG): model.gp FAILED"; grep -v Warning "logs/model_$NAME.log" | tail -3; exit 1; }
./run.sh reduce.gp "logs/reduce_$NAME.log" "NAME=\"$NAME\"" "BASE=\"$BASE\"" ${FRAMES:+"FRAMES=\"$FRAMES\""} > /dev/null || { echo "$NAME: reduce.gp FAILED"; grep '^  \*\*\*' "logs/reduce_$NAME.log" | grep -v Warning | head -3; exit 1; }
grep -q "^written out/model_${NAME}_red.gp" "logs/reduce_$NAME.log" || { echo "$NAME: reduce.gp FAILED"; grep -v Warning "logs/reduce_$NAME.log" | tail -3; exit 1; }
if [ "${MINRED:-1}" = 1 ]; then   # minimisation and reduction by Magma's MinRedTernaryForm (minred.m finds the change of coordinates, minred.gp applies it to the model, its lines and singular points); MINRED=0 gives the coordinates of the models built before this step existed
  rm -f "out/model_${NAME}_minred.gp"
  timeout 3600 magma -b "file:=out/model_${NAME}_red_coeffs.txt" "out:=out/model_${NAME}_minred.gp" minred.m > "logs/minredm_$NAME.log" 2>&1 < /dev/null || { echo "$NAME: minred.m FAILED (magma exit $?)"; tail -3 "logs/minredm_$NAME.log"; exit 1; }
  { grep -q "^written" "logs/minredm_$NAME.log" && [ -s "out/model_${NAME}_minred.gp" ]; } || { echo "$NAME: minred.m FAILED"; tail -3 "logs/minredm_$NAME.log"; exit 1; }
  ./run.sh minred.gp "logs/minred_$NAME.log" "NAME=\"$NAME\"" > /dev/null || { echo "$NAME: minred.gp FAILED"; grep '^  \*\*\*' "logs/minred_$NAME.log" | grep -v Warning | head -3; exit 1; }
  grep -q "^written out/model_${NAME}_red.gp (minimised)" "logs/minred_$NAME.log" || { echo "$NAME: minred.gp FAILED"; grep -v Warning "logs/minred_$NAME.log" | tail -3; exit 1; }
fi
rm -f "out/model_${NAME}_sing.gp"   # a stale singular-locus file must not satisfy a rerun
timeout 3600 magma -b "file:=out/model_${NAME}_red_coeffs.txt" "out:=out/model_${NAME}_sing.gp" sing.m > "logs/sing_$NAME.log" 2>&1 || { echo "$NAME: sing.m FAILED (magma exit $?)"; tail -3 "logs/sing_$NAME.log"; exit 1; }
grep -q "^written" "logs/sing_$NAME.log" || { echo "$NAME: sing.m FAILED"; tail -3 "logs/sing_$NAME.log"; exit 1; }
[ -s "out/model_${NAME}_sing.gp" ] || { echo "$NAME: sing.m wrote no output"; exit 1; }
./sqlines.sh "out/model_${NAME}_red_coeffs.txt" "$H" "out/model_${NAME}_sq$H.txt" 16 2> "logs/sqlines_$NAME.log" || { echo "$NAME: sqlines.sh FAILED"; exit 1; }
./run.sh finalize.gp "logs/finalize_$NAME.log" "NAME=\"$NAME\"" "H=$H" "WRITE=$WRITE" "BASE=\"$BASE\"" > /dev/null || { echo "$NAME: finalize.gp FAILED"; grep '^  \*\*\*' "logs/finalize_$NAME.log" | grep -v Warning | head -3; exit 1; }
{ grep -q "^TOTAL" "logs/finalize_$NAME.log" && grep -q "^written" "logs/finalize_$NAME.log"; } || { echo "$NAME: finalize.gp FAILED"; grep -v Warning "logs/finalize_$NAME.log" | tail -3; exit 1; }
T1=$(date +%s)
SUMMARY="$NAME ($I,$J,$SG): $(grep '^Magma' "logs/finalize_$NAME.log" | sed 's/Magma: //') | $(grep '^TOTAL' "logs/finalize_$NAME.log" | sed 's/TOTAL: //') | max coeff $(grep -o 'max |coeff| [0-9]*' "logs/finalize_$NAME.log" | head -1 | cut -d' ' -f3) | $(grep '^a_p' "logs/finalize_$NAME.log") | $((T1-T0)) s"
echo "$SUMMARY" | tee -a out/summary.txt

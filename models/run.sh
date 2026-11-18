#!/bin/bash
# run.sh SCRIPT.gp [LOGFILE] [gp assignments...] : run a gp script of this directory (all scripts assume this working
# directory), with timing, output to logs/<script>.log (or LOGFILE) and echoed.  Extra arguments are gp
# statements executed before the script, e.g.  ./run.sh model.gp logs/model_A.log 'NAME="A"' 'PAIR=[5,19,1]'.
# Exit status 1 if gp reported any error (a line starting with "  ***" that is not a warning): gp continues with the next
# statement of the file after an error, which is why the log is scanned (Codex round 1 of PR #3: the E2 validation had
# failed inside reduce.gp without the pipeline noticing).
cd "$(dirname "$0")" || exit 1
mkdir -p out logs
S=$1; shift
L=${1:-logs/$(basename "$S" .gp).log}; [ $# -gt 0 ] && shift
PRE=$(mktemp "${TMPDIR:-/tmp}/cp_pre_XXXXXX.gp")
echo "default(parisizemax, ${PARISIZEMAX:-8000000000});" > "$PRE"
for a in "$@"; do echo "$a;" >> "$PRE"; done
{ time gp -q "$PRE" "$S" < /dev/null; } > "$L" 2>&1
RC=$?   # Codex round 4 of PR #3: a gp that crashes or is killed without printing a *** diagnostic must fail the wrapper too
rm -f "$PRE"
cat "$L"
if [ "$RC" -ne 0 ]; then echo "run.sh: gp exited with status $RC in $S (see $L)" >&2; exit 1; fi
if grep '^  \*\*\*' "$L" | grep -qv 'Warning'; then echo "run.sh: gp reported an error in $S (see $L)" >&2; exit 1; fi
exit 0

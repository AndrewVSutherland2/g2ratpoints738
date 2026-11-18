#!/bin/bash
# minred_file.sh NAME SEXTIC LINESFILE "a,b,c;a,b,c" OUTDIR : minimised and reduced copy of an existing model, for searches that
# need comparable line heights across models (the pipeline does this for new models; see minred.m).  Writes
# OUTDIR/x163_NAME.txt (28 coefficients), OUTDIR/x163_NAME_lines.txt (the square lines and, in the header, the singular
# points and the change of coordinates) and prints "NAME maxcoeff_before maxcoeff_after detN maxline_before maxline_after
# singular_points".  Nothing in data/ is touched.  Exit status 1 on any failure.
set -uo pipefail
cd "$(dirname "$0")" || exit 1
NAME=$1; SX=$(readlink -f "$2"); LF=$(readlink -f "$3"); NODES=$4; OUT=$(readlink -f "$5"); mkdir -p "$OUT"
MRF=$OUT/minred_$NAME.gp; rm -f "$MRF" "$OUT/x163_$NAME.txt" "$OUT/x163_${NAME}_lines.txt"
timeout 3600 magma -b "file:=$SX" "out:=$MRF" minred.m > "$OUT/minred_$NAME.log" 2>&1 < /dev/null || { echo "$NAME: minred.m FAILED"; exit 1; }
[ -s "$MRF" ] || { echo "$NAME: minred.m wrote nothing"; exit 1; }
NV="[$(echo "$NODES" | sed 's/;/],[/g; s/^/[/; s/$/]/')]"
PRE=$(mktemp "$OUT/pre_XXXXXX.gp"); printf 'default(parisizemax, 2000000000);\nNAME = "%s"; SXF = "%s"; LNF = "%s"; NODES = %s; MRF = "%s"; OUTS = "%s"; OUTL = "%s";\n' "$NAME" "$SX" "$LF" "$NV" "$MRF" "$OUT/x163_$NAME.txt" "$OUT/x163_${NAME}_lines.txt" > "$PRE"
RES=$(gp -q "$PRE" minred_file.gp < /dev/null 2>&1); RC=$?; rm -f "$PRE"
if [ $RC -ne 0 ] || echo "$RES" | grep '^  \*\*\*' | grep -qv Warning || [ ! -s "$OUT/x163_$NAME.txt" ] || [ ! -s "$OUT/x163_${NAME}_lines.txt" ]; then echo "$NAME: minred_file.gp FAILED: $(echo "$RES" | grep -v Warning | tail -2 | tr '\n' ' ')"; exit 1; fi
echo "$RES" | grep -v Warning | tail -1

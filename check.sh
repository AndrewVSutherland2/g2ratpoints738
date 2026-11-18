#!/bin/bash
# check.sh : a short test of the installation (about three minutes on eight cores).  It runs the search on three models
# of data/models.txt, one of each kind, and compares the numbers and the curves with the stored results:
#   b4     one node; it carries the curve with 642 points of Elkies and Stoll
#   p144   two nodes; the curve with 638 points
#   t0453  three nodes; the curve C_738
# Needs what the search needs (the C programs, PARI/GP, Python 3), not Magma.  The run is left in a new directory
# under ${TMPDIR:-/tmp}, whose name is printed.
set -euo pipefail
TOP=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d "${TMPDIR:-/tmp}/check738.XXXXXX")
{ grep '^#' "$TOP/data/models.txt"; grep -E '^(b4|p144|t0453) ' "$TOP/data/models.txt"; } > "$W/models.txt"
echo "check.sh: searching b4, p144 and t0453 in $W/run"
JOBS=3 "$TOP/run_all.sh" "$W/run" "$W/models.txt"

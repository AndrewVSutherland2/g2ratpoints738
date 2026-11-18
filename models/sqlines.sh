#!/bin/bash
# sqlines.sh SEXTIC.txt H OUT.txt [threads=16] : exhaustive list of the primitive lines A x0 + B x1 + C x2 = 0 with
# max(|A|,|B|,|C|) <= H on which the sextic (28 int64 coefficients, linesearch order) restricts to a perfect square
# (tritangents and node lines; lines contained in the branch curve would also appear).  Enumeration by
# src/linesearch -enum -nosym, exact restriction by src/linesearch -exactbatch, square test by sqfilter (GMP), in parallel.
set -euo pipefail
SX=$(readlink -f "$1"); H=$2; OUT=$(readlink -f "$3"); P=${4:-16}
DIR=$(cd "$(dirname "$0")" && pwd); SRC=$(readlink -f "$DIR/../src")
[ -x "$DIR/sqfilter" ] || { echo "sqlines.sh: $DIR/sqfilter not built" >&2; exit 1; }
[ -x "$SRC/linesearch" ] || { echo "sqlines.sh: $SRC/linesearch not built" >&2; exit 1; }
T=$(mktemp -d "${TMPDIR:-/tmp}/sqlines_XXXXXX")
"$SRC/linesearch" "$SX" "$H" -enum -nosym > "$T/lines.txt"
NL=$(wc -l < "$T/lines.txt")
split -n l/$((P*4)) -d -a 3 "$T/lines.txt" "$T/chunk."
export SX SRC DIR
ls "$T" | grep '^chunk\.' | xargs -P "$P" -I{} bash -c 'set -o pipefail; "$SRC/linesearch" "$SX" 1 -exactbatch < "'"$T"'/{}" | "$DIR/sqfilter" > "'"$T"'/{}.sq" 2> "'"$T"'/{}.err"'
cat "$T"/chunk.*.sq | sort -t' ' -k1,1n -k2,2n -k3,3n > "$OUT"
echo "sqlines: $NL lines with max coefficient <= $H enumerated, $(wc -l < "$OUT") square lines written to $OUT" >&2
rm -rf "$T"

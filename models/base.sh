#!/bin/bash
# base.sh B [T] : everything about a one-node base B (E = Elkies' model: data/elkies163*.txt; else data/x163_<B>*.txt,
# node at (0:0:1)): the table of conics (if not there yet), the lattice, the conic pairs with at most three orthogonal
# roots, and their canonical forms.  T threads / pieces.  Run from the root; results in models/lat.
set -euo pipefail
B=$1; T=${2:-16}; D=$(cd "$(dirname "$0")" && pwd); R=$(pwd); L=$D/lat; C=$D/cw; mkdir -p "$L" "$C"
if [ "$B" = E ]; then F=data/elkies163; else F=data/x163_$B; fi
t0=$(date +%s.%N)
if [ ! -s "${F}_conics.txt" ]; then
  : > "$C/${B}_lifted.txt"
  for p in 41 43 47; do                 # one prime is enough when the lattice confirms the table (base.gp); else add a prime
    "$D/conicmodp" "$F.txt" $p -t "$T" > "$C/${B}_mod$p.txt"
    python3 "$D/lift_conics.py" "$F.txt" $p "$C/${B}_mod$p.txt" -K 80 -o "$C/${B}_lifted_$p.txt" > "$C/${B}_lift_$p.log" 2>&1
    sort -t'|' -k1,1 -u "$C/${B}_lifted.txt" "$C/${B}_lifted_$p.txt" > "$C/${B}_lifted.new"; mv "$C/${B}_lifted.new" "$C/${B}_lifted.txt"
    echo "read(\"$D/verify_conics.gp\"); main(\"$F.txt\", \"$C/${B}_lifted.txt\", \"$C/${B}_param.txt\")" | gp -q > "$C/${B}_verify.log" 2>&1
    python3 "$D/conic_table.py" "$C/${B}_param.txt" "${F}_conics.txt" > "$C/${B}_mkconics.log"
    printf 'BNAME="%s";\nOUTD="%s";\n' "$B" "$L" > "$L/pre_$B.gp"
    gp -q "$L/pre_$B.gp" "$D/base.gp" < /dev/null > "$L/base_$B.log" 2>&1 || true
    grep -q "^BASEDONE $B" "$L/base_$B.log" && break
    grep -q "conic table is not complete\|is not a conic of the table" "$L/base_$B.log" || { echo "base.sh: $B: base.gp failed: $(grep -v Warning "$L/base_$B.log" | tail -2)" >&2; exit 1; }
    echo "$B: the conics found modulo $p do not fill the lattice count; adding a prime" >&2
  done
else
  printf 'BNAME="%s";\nOUTD="%s";\n' "$B" "$L" > "$L/pre_$B.gp"
  gp -q "$L/pre_$B.gp" "$D/base.gp" < /dev/null > "$L/base_$B.log" 2>&1 || true
fi
grep -q "^BASEDONE $B" "$L/base_$B.log" || { echo "base.sh: $B: no complete lattice: $(grep -v Warning "$L/base_$B.log" | tail -2)" >&2; exit 1; }
t1=$(date +%s.%N)
python3 "$D/pairs.py" "$L/lat_$B.txt" "$L/dp_$B.txt" 3 > "$L/pairs_$B.log"
t2=$(date +%s.%N)
"$D/classify.sh" "$B" "$T" > "$L/classify_$B.log"
t3=$(date +%s.%N)
"$D/base_class.sh" "$B" > "$L/selfhash_$B.txt"
echo "$B: $(grep -c . "${F}_conics.txt") conics; $(tail -1 "$L/pairs_$B.log" | sed 's/^[^:]*: //'); $(sed 's/^[^:]*: //' "$L/classify_$B.log"); seconds: conics and lattice $(echo "$t1 - $t0" | bc), pairs $(echo "$t2 - $t1" | bc), canonical forms $(echo "$t3 - $t2" | bc)"

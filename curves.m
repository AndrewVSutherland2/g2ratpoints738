// curves.m -- the curves found by the search (data/search_curves.txt), checked against the models (data/models.txt).
//   magma -b curves.m                 (all rows)
//   magma -b lo:=1 hi:=50 curves.m    (a range of rows; the summary is skipped)
//   magma -b file:=F classes:=G curves.m    (rows of F without the last two columns, as written by search/collect.py:
//                                      the two columns are computed and written to G, one row "Galois group, class" per curve)
// For each row "points, points of height <= 10^5, model, line a b c, sextic g0 .. g6, score, rank, points of height
// <= 2000, <= 10^4, Galois group, class" it proves:
//  1. the curve y^2 = g(x) is isomorphic over Q to the section of the model by the line: with p1, p2 a basis of the
//     points of the line, y^2 = S(x p1 + p2).  (So g is that binary sextic up to a change of variable in GL_2(Q) and
//     a square factor: the row describes the curve cut out by the line, not a twist of it.)
//  2. the Galois group of g is the one stated.  For the rows marked S6 it is the symmetric group; by Zarhin's theorem
//     the Jacobian then has endomorphism ring Z over an algebraic closure, so it is absolutely simple and the curve has
//     no automorphism besides the hyperelliptic involution.  For any other row the script computes the automorphism
//     group over an algebraic closure and checks that it contains an involution other than the hyperelliptic one, so
//     that the Jacobian is isogenous to a product of two elliptic curves.
//  3. two rows have the same class number exactly when their curves are isomorphic over an algebraic closure of Q
//     (equal absolute invariants).
// It does not recount the points, the scores or the ranks: these are the output of the search programs
// (search/search_model.sh regenerates them).  The summary gives the numbers of classes quoted in the note.
SetColumns(0); SetQuitOnError(true);
Q := Rationals(); P<z0,z1,z2> := PolynomialRing(Q, 3); R<x> := PolynomialRing(Q);
mons := []; for d in [0..6] do for j in [0..d] do Append(~mons, z0^(6-d)*z1^(d-j)*z2^j); end for; end for;
ints := func< s | [StringToInteger(t) : t in Split(s, " ,") | #t gt 0] >;
model := AssociativeArray(); nodes := AssociativeArray();
for r in Split(Read("data/models.txt"), "\n") do
  if #r eq 0 or r[1] eq "#" then continue; end if;
  F := Split(r, "|"); name := &cat Split(F[1], " "); c := ints(F[3]); assert #c eq 28;
  model[name] := &+[c[i]*mons[i] : i in [1..28]]; nodes[name] := StringToInteger(&cat Split(F[2], " "));
end for;
if not assigned file then file := "data/search_curves.txt"; end if;
if not assigned classes then classes := ""; end if;
writing := classes ne "";
rows := [Split(r, " ") : r in Split(Read(file), "\n") | #r gt 0 and r[1] ne "#"];
assert forall{r : r in rows | #r eq (writing select 17 else 19)};
full := not assigned lo and not assigned hi;
a := assigned lo select StringToInteger(lo) else 1; b := assigned hi select StringToInteger(hi) else #rows;
inv := []; nS6 := 0; other := []; gal := [];
for i in [a..b] do
  r := rows[i]; S := model[r[3]]; l := [StringToInteger(r[k]) : k in [4..6]]; g := R![StringToInteger(r[k]) : k in [7..13]];
  K := Basis(Kernel(Matrix(Q, 3, 1, l))); p1 := Eltseq(K[1]); p2 := Eltseq(K[2]);
  f := Evaluate(S, [p1[k]*x + p2[k] : k in [1..3]]);
  error if Degree(f) lt 5 or Discriminant(f) eq 0, "row", i, ": the line does not cut out a curve of genus 2";
  Cf := HyperellipticCurve(f); Cg := HyperellipticCurve(g);
  error if not IsIsomorphic(Cf, Cg), "row", i, ": the sextic is not the section of the model by the line";
  G := GaloisGroup(g);
  if writing then
    if #G eq 720 then Append(~r, "S6"); else n, d := TransitiveGroupIdentification(G); assert d eq 6; Append(~r, "6T" cat IntegerToString(n)); end if;
    Append(~gal, r[18]);
  end if;
  if r[18] eq "S6" then
    error if #G ne 720, "row", i, ": the Galois group is not S6";
    nS6 +:= 1;
  else
    n, d := TransitiveGroupIdentification(G);
    error if d ne 6 or r[18] ne "6T" cat IntegerToString(n), "row", i, ": the Galois group is not the one stated";
    A := GeometricAutomorphismGroup(Cg);
    error if #[t : t in A | Order(t) eq 2] lt 2, "row", i, ": no involution besides the hyperelliptic one";
    Append(~other, <r[1], r[3], r[18], #G, #A>);
  end if;
  Append(~inv, G2Invariants(Cg));
end for;
if not full then printf "rows %o..%o: no failure\n", a, b; quit; end if;
// 3. the classes
if writing then          // classes numbered in the order of the rows
  num := AssociativeArray(); cls := [];
  for i in [1..#rows] do if not IsDefined(num, inv[i]) then num[inv[i]] := #Keys(num) + 1; end if; Append(~cls, num[inv[i]]); end for;
  PrintFile(classes, &cat[Sprintf("%o %o\n", gal[i], cls[i]) : i in [1..#rows]] : Overwrite := true);
else
  cls := [StringToInteger(r[19]) : r in rows];
end if;
byinv := AssociativeArray();
for i in [1..#rows] do
  if IsDefined(byinv, inv[i]) then assert byinv[inv[i]] eq cls[i]; else byinv[inv[i]] := cls[i]; end if;
end for;
assert #Keys(byinv) eq #Set(cls);
printf "%o rows: each sextic is the section of its model by its line; %o have Galois group S6\n", #rows, nS6;
for t in other do printf "the curve with %o points on %o: Galois group %o of order %o, automorphism group of order %o with an extra involution\n", t[1], t[2], t[3], t[4], t[5]; end for;
printf "%o classes of curves up to isomorphism over an algebraic closure\n", #Set(cls);
for k in [1..3] do
  best := AssociativeArray();
  for i in [1..#rows] do if nodes[rows[i][3]] eq k then
    n := StringToInteger(rows[i][1]); if not IsDefined(best, cls[i]) or best[cls[i]] lt n then best[cls[i]] := n; end if;
  end if; end for;
  v := [best[c] : c in Keys(best)];
  printf "%o node%o: %o curves closed, %o classes; with at least 500, 550, 600 points: %o, %o, %o; most points %o\n", k, k eq 1 select "" else "s",
    #[i : i in [1..#rows] | nodes[rows[i][3]] eq k], #v, #[n : n in v | n ge 500], #[n : n in v | n ge 550], #[n : n in v | n ge 600], Max(v);
end for;
onmodels := AssociativeArray(); nrows := AssociativeArray();
for i in [1..#rows] do
  if not IsDefined(onmodels, cls[i]) then onmodels[cls[i]] := {}; nrows[cls[i]] := 0; end if;
  Include(~onmodels[cls[i]], rows[i][3]); nrows[cls[i]] +:= 1;
end for;
printf "classes on more than one model: %o; classes that occur more than once on their model: %o (%o of them on Elkies' sextic E)\n",
  #[c : c in Keys(onmodels) | #onmodels[c] gt 1], #[c : c in Keys(nrows) | nrows[c] gt 1], #[c : c in Keys(nrows) | nrows[c] gt 1 and onmodels[c] eq {"E"}];
print "CURVES CHECK PASSED";
quit;

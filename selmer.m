// selmer.m -- the 2-Selmer group of the Jacobian of C_738, with class group computations under GRH: it has rank 21, so
// rank J(Q) <= 21 under GRH (J(Q)[2] = 0).      magma -b selmer.m      (about four hours on one core; output in selmer.log)
SetColumns(0); SetQuitOnError(true);
SetClassGroupBounds("GRH");
R<x> := PolynomialRing(Rationals());
f := R![6838225,-46160230,40133975,172162750,60026500,-16144220,743425]; C := HyperellipticCurve(f); J := Jacobian(C);
t := Cputime();
S, m := TwoSelmerGroup(J);
printf "2-Selmer group: %o generators, invariants %o (GRH), %o s\n", Ngens(S), Invariants(S), Cputime(t);
printf "RankBound %o\n", RankBound(J);
quit;

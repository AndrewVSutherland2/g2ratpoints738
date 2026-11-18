// mass.m -- the masses quoted in the section "Models as lattice classes" (magma -b mass.m; a second)
// L0 = <2> + E8 + E8 + T : the genus of D^perp(-1) for the classes D with D^2 = 2 (rank 19, det 326)
// W0 = E8 + E8 + T      : the genus of M(-1) for the excluded classes, D^perp = <r> + M (rank 18, det 163)
T := LatticeWithGram(Matrix(Integers(), 2, 2, [2,1,1,82]));
E8 := Lattice("E", 8);
W0 := DirectSum(DirectSum(E8, E8), T);
L0 := DirectSum(LatticeWithGram(Matrix(Integers(), 1, 1, [2])), W0);
assert Rank(W0) eq 18 and Determinant(W0) eq 163 and Rank(L0) eq 19 and Determinant(L0) eq 326;
mW := Mass(W0); mL := Mass(L0); R := RealField(12);
printf "mass of the genus of L0 = %o = %o\n", mL, R!mL;
printf "mass of the genus of W0 = %o = %o\n", mW, R!mW;
printf "mass of the polarizations = mass(L0) - mass(W0)/2 = %o\n", R!(mL - mW/2);
printf "lower bound for the number of models = 2 * that = %o\n", R!(2*(mL - mW/2));
assert Abs(R!mL - 16417.77) lt 0.01 and Abs(R!mW - 2.22) lt 0.01 and Abs(R!(mL - mW/2) - 16416.66) lt 0.01 and 2*(mL - mW/2) gt 32800;
print "MASS CHECK PASSED";
quit;

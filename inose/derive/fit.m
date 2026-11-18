load "pencil.m";
function LoadPencil(fn)
  d := eval Read(fn); MA := d[1]; MB := d[2]; res := [];
  for v in d[3] do
    A := &+[R | v[i]*x0^MA[i][1]*x1^MA[i][2]*x2^MA[i][3] : i in [1..#MA]];
    B := &+[R | v[#MA+i]*x0^MB[i][1]*x1^MB[i][2]*x2^MB[i][3] : i in [1..#MB]];
    Append(~res, <A, B>);
  end for;
  return res;
end function;
// value of A + B w at a sample of a lifted line
function AtLift(N, i, sgn, vv)
  l := LIN[i]; P0 := [l[2][m]*vv[1] + l[3][m]*vv[2] : m in [1..3]]; W0 := sgn*Evaluate(l[4], vv);
  return Evaluate(N[1], P0) + Evaluate(N[2], P0)*W0;
end function;
// the member of the pencil vanishing on a lifted line
function Member(pen, i, sgn)
  vv := [101, 37];                    // a point of the line that is on no other curve used here
  v1 := AtLift(pen[1], i, sgn, vv); v2 := AtLift(pen[2], i, sgn, vv);
  assert v1 ne 0 or v2 ne 0;
  N := <v2*pen[1][1] - v1*pen[2][1], v2*pen[1][2] - v1*pen[2][2]>;
  assert AtLift(N, i, sgn, [89, -53]) eq 0 and AtLift(N, i, sgn, [71, 113]) eq 0;
  c := LCM([Denominator(x) : x in Coefficients(N[1]) cat Coefficients(N[2])]); g := GCD([Integers() ! (c*x) : x in Coefficients(N[1]) cat Coefficients(N[2])]);
  return <N[1]*c/g, N[2]*c/g>;
end function;
tp := LoadPencil("t_pencil.txt"); up := LoadPencil("u_pencil.txt");
Nt := Member(tp, 31, -1);      // vanishes on a1: the fibre F_a
Gb := Member(tp, 32, 1);       // vanishes on b1: the fibre F_b
print "G_b has no w-part:", Gb[2] eq 0;
lf := func<i | LIN[i][1][1]*x0 + LIN[i][1][2]*x1 + LIN[i][1][3]*x2>;
qb7 := 768*x0^2 + 1960*x0*x1 - 100*x0*x2 + 1152*x1^2 - 117*x1*x2 + 3*x2^2;
GbP := lf(52)*lf(32)^2*lf(29)^3*lf(39)^4*lf(17)^5*lf(33)^6*lf(37)^4*qb7^2*lf(15)^3;
print "G_b is the product of the components of F_b:", IsDivisibleBy(Gb[1], GbP) and Degree(Gb[1]) eq Degree(GbP), Gb[1] div GbP;
qa7 := 768*x0^2 - 1960*x0*x1 + 100*x0*x2 + 1152*x1^2 - 117*x1*x2 + 3*x2^2;
GaP := lf(51)*lf(31)^2*lf(28)^3*lf(38)^4*lf(16)^5*lf(30)^6*lf(36)^4*qa7^2*lf(14)^3;
nrm := Nt[1]^2 - Nt[2]^2*S;
print "norm of N_t = c * G_a * G_b:", IsDivisibleBy(nrm, GaP*GbP), nrm div (GaP*GbP);
N0 := Member(up, 36, -1);      // vanishes on a6: the III* fibre of the second pencil
Ni := Member(up, 32, 1);       // vanishes on b1: the I6* fibre
PrintFile("maps_raw.txt", <Nt, GbP, N0, Ni> : Overwrite := true);
print "digits:", [#IntegerToString(Max([Abs(Integers() ! c) : c in Coefficients(f)])) : f in [Nt[1], Nt[2], N0[1], N0[2], Ni[1], Ni[2]]];
quit;

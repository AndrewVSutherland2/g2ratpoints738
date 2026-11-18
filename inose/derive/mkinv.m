// writes inverse_data.m from the reduced solutions of the interpolation (first vector)
V := eval Read("inv_lll.txt"); v := V[1];
R2<t,u> := PolynomialRing(Rationals(), 2);
KK := 34; EA := [[i, j] : j in [0..KK-i], i in [0..KK]]; EB := [[i, j] : j in [0..KK-3-i], i in [0..KK-3]]; n := #EA + #EB;
assert #v eq 3*n;
SetColumns(0); fn := "inverse_data.m";
PrintFile(fn, "// inverse map: (x0 : x1 : x2) = (X0 : X1 : X2), X_j = XA[j] + XB[j]*w' on w'^2 = t(u^3(t+1) - 3 a u^2 t - t^2); see README.txt" : Overwrite := true);
for j in [0..2] do
  A := &+[v[j*n + i]*t^EA[i][1]*u^EA[i][2] : i in [1..#EA]]; B := &+[v[j*n + #EA + i]*t^EB[i][1]*u^EB[i][2] : i in [1..#EB]];
  PrintFile(fn, Sprintf("X%oA := %o;", j, A)); PrintFile(fn, Sprintf("X%oB := %o;", j, B));
  print j, Degree(A), Degree(B), #Terms(A), #Terms(B), #IntegerToString(Max([Abs(Integers() ! c) : c in Coefficients(A) cat Coefficients(B)]));
end for;
quit;

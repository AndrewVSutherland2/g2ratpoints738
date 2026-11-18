// writes map_data.m (the forms of the map, readable by Magma) from the two pencils computed by t_q.m and u_q.m
load "pencil.m";
d := eval Read("maps_raw.txt"); Nt := d[1]; N0 := d[3]; Ni := d[4];
// t = -N_t/G_b and u = N_0/(2 N_inf) satisfy the equation; absorb the constants
At := -Nt[1]; Bt := -Nt[2]; A0 := N0[1]; B0 := N0[2]; Ai := 2*Ni[1]; Bi := 2*Ni[2];
SetColumns(0);
fn := "map_data.m";
PrintFile(fn, "// forms on Elkies' sextic w^2 = S(x0,x1,x2): t = (At + Bt*w)/Gb, u = (A0 + B0*w)/(Ai + Bi*w); see README.txt" : Overwrite := true);
PrintFile(fn, "At := " cat Sprint(At) cat ";"); PrintFile(fn, "Bt := " cat Sprint(Bt) cat ";");
PrintFile(fn, "A0 := " cat Sprint(A0) cat ";"); PrintFile(fn, "B0 := " cat Sprint(B0) cat ";");
PrintFile(fn, "Ai := " cat Sprint(Ai) cat ";"); PrintFile(fn, "Bi := " cat Sprint(Bi) cat ";");
print [Degree(f) : f in [At, Bt, A0, B0, Ai, Bi]], [#Terms(f) : f in [At, Bt, A0, B0, Ai, Bi]];
print [#IntegerToString(Max([Abs(Integers() ! c) : c in Coefficients(f)])) : f in [At, Bt, A0, B0, Ai, Bi]];
quit;

// verify.m : the map from Elkies' sextic model of X_163 to the Inose surface.     magma -b verify.m   (a few seconds)
//
//   X : w^2 = S(x0,x1,x2)                      Elkies' sextic (Clay note), node at (0:0:1)
//   Y : w'^2 = t (u^3 (t+1) - 3 a u^2 t - t^2)  double plane of the Inose surface, a = 53360
//   Inose:  y^2 = x^3 - 3 a^2 x - (t + 1/t + 2 a^3 + 2),   x = u (t+1)/t - a,   y = w' (t+1)/t^2
//
//   t = (At + Bt w)/Gb,   u = (A0 + B0 w)/(Ai + Bi w),   w' = w J/8,   J = d(t,u)/d(x1,x2) on the chart x0 = 1.
//
// (1) the fibre t = 0 is the divisor F_a: eight lifted lines and a lifted conic, meeting as the components of a fibre of
//     type II*, with the multiplicities of that fibre; t = infinity is its image F_b under (x1,x2) -> (-x1,-x2);
// (2) the equation of Y holds identically on X;
// (3) the Weierstrass equation follows from the equation of Y;
// (4) the map is birational: the fibres of t are connected (F_a is connected with a reduced component), and u has
//     degree at most F.(24 D - Z) = 2 on them, where Z is the base divisor of the pencil of u (fifteen lifted lines).
//     Since w' is not in Q(t,u), the function field of X is Q(t,u,w').
Q := Rationals();
R<x0,x1,x2> := PolynomialRing(Q, 3);
S := 3686400*x0^6 + (-8697600*x1^2 + 2193664*x2*x1 + 11520*x2^2)*x0^4 + (5130225*x1^4 - 1860100*x2*x1^3 + 138414*x2^2*x1^2 - 4180*x2^3*x1 + 9*x2^4)*x0^2 + (-45144*x2*x1^5 + 19960*x2^2*x1^4 - 1672*x2^3*x1^3 + 40*x2^4*x1^2);
a := 53360;
load "map_data.m";
ok := true;
procedure Check(name, cond, ~ok) printf "%o: %o\n", name, cond select "ok" else "FAILED"; ok := ok and cond; end procedure;
PV<v> := PolynomialRing(Q);
lin := func<c | c[1]*x0 + c[2]*x1 + c[3]*x2>;
// a line l = [l0,l1,l2] (l2 != 0 or l1 != 0): points P0 + v P1, a transversal vector n, and the cubic c(v) with S = c^2
function OnLine(l)
  if l[3] ne 0 then P0 := [l[3], 0, -l[1]]; P1 := [0, l[3], -l[2]]; n := [0, 0, 1];
  else P0 := [l[2], -l[1], 0]; P1 := [0, 0, 1]; n := [0, 1, 0]; end if;
  fl, c := IsSquare(Evaluate(S, [P0[m] + v*P1[m] : m in [1..3]])); assert fl;
  return <l, P0, P1, n, c>;
end function;
// order of vanishing >= m of A + B w along the lift w = e c(v) of a line (exact, v an indeterminate)
function VanishesOnLift(A, B, L, e, m)
  K := FieldOfFractions(PV); Pe<eps> := PolynomialRing(K); Te := quo<Pe | eps^m>;
  X := [Te | K ! (L[2][i] + v*L[3][i]) + L[4][i]*eps : i in [1..3]]; c := K ! L[5];
  del := Evaluate(S, X)/c^2 - 1; sq := Te ! 1; pw := Te ! 1; bc := K ! 1;
  for j in [1..m-1] do bc := bc*(3 - 2*j)/(2*j); pw := pw*del; sq := sq + bc*pw; end for;
  return Evaluate(A, X) + Evaluate(B, X)*e*c*sq eq 0;
end function;
// value of w on the lift (L, e) at the point s0 P0 + s1 P1
wAt := func<L, e, s | e*&+[Coefficient(L[5], m)*s[2]^m*s[1]^(3-m) : m in [0..3]]>;
// intersection number of two lifted lines
function LL(L1, e1, L2, e2)
  if L1[1] eq L2[1] then return e1 eq e2 select -2 else (L1[1][3] eq 0 select 2 else 3); end if;
  pt := Eltseq(Basis(Kernel(Transpose(Matrix(Q, 2, 3, [L1[1], L2[1]]))))[1]);
  if pt[1] eq 0 and pt[2] eq 0 then return 0; end if;       // the node: the two lifts are separated by the blow-up
  ws := [wAt(L, e, Solution(Matrix(Q, 2, 3, [L[2], L[3]]), Vector(Q, pt))) : L in [L1, L2]] where e is 1;
  w1 := e1*ws[1]; w2 := e2*ws[2]; assert w1^2 eq w2^2 and w1 ne 0;
  return w1 eq w2 select 1 else 0;
end function;

// ---- (1) the fibre F_a ----
la := [ <[9,-1,0], 1>, <[80,-11,1], 2>, <[72,-27,1], 3>, <[112,-49,3], 4>, <[33,-22,1], 5>, <[80,-57,3], 6>, <[108,-99,5], 4>, <[24,-19,1], 3> ];
qa := 768*x0^2 - 1960*x0*x1 + 100*x0*x2 + 1152*x1^2 - 117*x1*x2 + 3*x2^2;      // the conic, multiplicity 2
sig := hom<R -> R | x0, -x1, -x2>;
Ga := &*[lin(f[1])^f[2] : f in la]*qa^2; Gb := sig(Ga);
Check("S is invariant under (x1,x2) -> (-x1,-x2)", sig(S) eq S, ~ok);
Check("Ga is a form of degree 32", Degree(Ga) eq 32 and IsHomogeneous(Ga), ~ok);
Check("norm of the numerator of t: At^2 - Bt^2 S = Ga Gb", At^2 - Bt^2*S eq Ga*Gb, ~ok);
FA := [];        // the lifted lines of F_a: <line data, sign, multiplicity>
for f in la do
  L := OnLine(f[1]); ep := [e : e in [1, -1] | VanishesOnLift(At, Bt, L, e, 1)];
  Check(Sprintf("the numerator of t vanishes on exactly one lift of the line %o", f[1]), #ep eq 1, ~ok);
  Append(~FA, <L, ep[1], f[2]>);
end for;
Pc := [-34*v^2 - 67*v - 34, -16*v^2 - 60*v - 44, 416*v^2 + 264*v - 140];
fl, cc := IsSquare(Evaluate(S, Pc)); epc := [e : e in [1, -1] | Evaluate(At, Pc) + e*Evaluate(Bt, Pc)*cc eq 0];
Check("the conic is parametrised, S is a square on it, and the numerator of t vanishes on exactly one lift of it", Evaluate(qa, Pc) eq 0 and fl and #epc eq 1, ~ok);
// intersection number of the lifted conic with a lifted line
function CL(L, e)
  h := Evaluate(lin(L[1]), Pc); n := 0;
  for fc in Factorisation(h) do
    K := Degree(fc[1]) eq 1 select Q else NumberField(fc[1]); rt := Roots(fc[1], K)[1][1];
    pt := [Evaluate(Pc[m], rt) : m in [1..3]];
    sol := Solution(Matrix(K, 2, 3, L[2] cat L[3]), Vector(K, pt));
    wl := wAt(L, e, sol); wc := epc[1]*Evaluate(cc, rt); assert wl^2 eq wc^2 and wl ne 0;
    if wl eq wc then n +:= Degree(fc[1])*fc[2]; end if;
  end for;
  return n;
end function;
// order: a0, a1, ..., a6, a7 (conic), a8
I := ZeroMatrix(Integers(), 9, 9); pos := [1,2,3,4,5,6,7,9];
for i in [1..8] do for j in [1..8] do I[pos[i], pos[j]] := LL(FA[i][1], FA[i][2], FA[j][1], FA[j][2]); end for;
  I[pos[i], 8] := CL(FA[i][1], FA[i][2]); I[8, pos[i]] := I[pos[i], 8]; end for; I[8,8] := -2;
E8t := ZeroMatrix(Integers(), 9, 9);
for e in [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[6,9]] do E8t[e[1],e[2]] := 1; E8t[e[2],e[1]] := 1; end for;
for i in [1..9] do E8t[i,i] := -2; end for;
Check("the nine components of F_a meet as in the diagram of a fibre of type II*", I eq E8t, ~ok);
mult := Vector(Integers(), [1,2,3,4,5,6,4,2,3]);
Check("F_a . (each of its components) = 0, and F_a has a component of multiplicity 1", IsZero(mult*I) and mult[1] eq 1, ~ok);

// ---- (2) the equation of Y, on 260 lines x1 = c (both sides are polynomials of degree <= 244 in x1) ----
Py<y> := PolynomialRing(Q);
D1 := func<f | Derivative(f, 2)>; D2 := func<f | Derivative(f, 3)>;
forms := [* S, At, Bt, Gb, A0, B0, Ai, Bi *];
ders := [* <D1(f), D2(f)> : f in forms *];
nbad := 0; nlines := 0;
for c in [-130..129] do
  ev := func<f | Evaluate(f, [Py | 1, c, y])>;
  s := ev(S); sx := ev(ders[1][1]); sy := ev(ders[1][2]);
  Kw<w> := quo<PolynomialRing(Py) | $.1^2 - s>;
  pr := func<i | ev(forms[i]) + ev(forms[i+1])*w>;                       // the pair A + B w
  dx := func<i, k | 2*s*ev(ders[i][k]) + (2*s*ev(ders[i+1][k]) + ev(forms[i+1])*(k eq 1 select sx else sy))*w>;   // 2 S d/dx_k (A + B w)
  T := pr(2); U := pr(5); V := pr(7); G := Kw ! ev(Gb);
  Gx := Kw ! (2*s*ev(ders[4][1])); Gy := Kw ! (2*s*ev(ders[4][2]));
  Jn := (dx(2,1)*G - T*Gx)*(dx(5,2)*V - U*dx(7,2)) - (dx(2,2)*G - T*Gy)*(dx(5,1)*V - U*dx(7,1));    // J = Jn/(4 S^2 Gb^2 V^2)
  lhs := Jn^2;
  rhs := 1024*s^3*G*V*T*(U^3*(T + G)*G - 3*a*U^2*T*V*G - T^2*V^3);
  nlines +:= 1; if lhs ne rhs then nbad +:= 1; end if;
  if c eq 0 then Check("the Jacobian of (t,u) is not zero", Jn ne 0, ~ok); end if;
end for;
Check(Sprintf("w'^2 = t (u^3 (t+1) - 3 a u^2 t - t^2) on %o lines x1 = c", nlines), nbad eq 0, ~ok);

// ---- (3) from the double plane to the Weierstrass equation ----
F3<tt,uu> := FunctionField(Q, 2);
xx := uu*(tt+1)/tt - a;
Check("x^3 - 3 a^2 x - (t + 1/t + 2 a^3 + 2) = ((t+1)/t^2)^2 t (u^3 (t+1) - 3 a u^2 t - t^2)",
      xx^3 - 3*a^2*xx - (tt + 1/tt + 2*a^3 + 2) eq ((tt+1)/tt^2)^2*tt*(uu^3*(tt+1) - 3*a*uu^2*tt - tt^2), ~ok);
Check("t (u^3 (t+1) - 3 a u^2 t - t^2) is not a square in Q(t,u)", not IsSquare(Numerator(tt*(uu^3*(tt+1) - 3*a*uu^2*tt - tt^2))), ~ok);

// ---- (4) degree of u on the fibres of t ----
// base divisor Z of the pencil spanned by A0 + B0 w and Ai + Bi w: lifted lines with multiplicities
lz := [ <[24,-19,1], 1>, <[24,19,-1], 1>, <[33,-22,1], 2>, <[33,22,-1], 3>, <[34,-11,3], 1>, <[58,57,-1], 1>, <[72,-27,1], 1>, <[72,27,-1], 2>,
        <[80,-57,3], 2>, <[80,11,-1], 1>, <[80,57,-3], 3>, <[108,-99,5], 1>, <[108,99,-5], 2>, <[112,-49,3], 1>, <[112,49,-3], 3> ];
ZU := []; good := true;
for f in lz do
  L := OnLine(f[1]); ep := [e : e in [1, -1] | VanishesOnLift(A0, B0, L, e, f[2]) and VanishesOnLift(Ai, Bi, L, e, f[2])];
  good := good and #ep ge 1; if #ep ge 1 then Append(~ZU, <L, ep[1], f[2]>); end if;
end for;
Check("both forms of u vanish along the fifteen lifted lines of Z to the stated orders", good, ~ok);
Check("Z has degree 25 = 2*24 - 23", &+[z[3] : z in ZU] eq 25, ~ok);
FZ := &+[FA[i][3]*z[3]*LL(FA[i][1], FA[i][2], z[1], z[2]) : i in [1..8], z in ZU] + &+[2*z[3]*CL(z[1], z[2]) : z in ZU];
Check(Sprintf("F_a . (24 D - Z) = 24*32 - %o = 2", FZ), 24*32 - FZ eq 2, ~ok);
printf "%o\n", ok select "INOSE MAP CHECK PASSED" else "INOSE MAP CHECK FAILED";
quit;

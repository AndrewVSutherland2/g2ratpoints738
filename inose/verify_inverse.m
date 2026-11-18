// verify_inverse.m : the forms of inverse_data.m give the inverse of the map of verify.m.
//   magma -b LO:=1 HI:=100 verify_inverse.m      (fibres number LO..HI of the list; the whole list is 1..100)
//
// Let phi : X -> Y be the map of verify.m (birational) and psi = (X0 : X1 : X2), forms of degree 34 in (t, u, w').
// For a point q of Y the program computes P = psi(q) exactly and checks phi(P) = q.  Then psi(q) = phi^-1(q).
// It does this at 102 points (51 values of u and both values of w') on each of 100 smooth fibres t = 2, ..., 101.
// The sections X_i xi_j - X_j xi_i (xi = phi^-1, sections of D_E) have class 34 D_I + D_E; on a fibre F this has degree
// 34*2 + 32 = 100 < 102, so they vanish on each of these fibres; and (34 D_I + D_E - 100 F).D_I <= 68 + 112 - 200 < 0
// (D_I.D_E <= 112 because the map to Y is given by forms of degree 56 on X), so they vanish identically: psi = phi^-1.
Q := Rationals();
R<x0,x1,x2> := PolynomialRing(Q, 3);
S := 3686400*x0^6 + (-8697600*x1^2 + 2193664*x2*x1 + 11520*x2^2)*x0^4 + (5130225*x1^4 - 1860100*x2*x1^3 + 138414*x2^2*x1^2 - 4180*x2^3*x1 + 9*x2^4)*x0^2 + (-45144*x2*x1^5 + 19960*x2^2*x1^4 - 1672*x2^3*x1^3 + 40*x2^4*x1^2);
a := 53360;
load "map_data.m";
la := [ <9*x0 - x1, 1>, <80*x0 - 11*x1 + x2, 2>, <72*x0 - 27*x1 + x2, 3>, <112*x0 - 49*x1 + 3*x2, 4>, <33*x0 - 22*x1 + x2, 5>,
        <80*x0 - 57*x1 + 3*x2, 6>, <108*x0 - 99*x1 + 5*x2, 4>, <768*x0^2 - 1960*x0*x1 + 100*x0*x2 + 1152*x1^2 - 117*x1*x2 + 3*x2^2, 2>, <24*x0 - 19*x1 + x2, 3> ];
Gb := Evaluate(&*[f[1]^f[2] : f in la], [x0, -x1, -x2]);
R2<t,u> := PolynomialRing(Q, 2);
load "inverse_data.m";
forms := [S, At, Bt, Gb, A0, B0, Ai, Bi];
d1 := [Derivative(f, 2) : f in forms]; d2 := [Derivative(f, 3) : f in forms];
lo := StringToInteger(LO); hi := StringToInteger(HI);
// phi(psi(q)) = q for q = (t1, u1, wp), exactly, in the field K of q; -1 if psi or phi is not defined there
function TestPoint(K, t1, u1, wp)
  P := [K | Evaluate(X0A, [t1,u1]) + Evaluate(X0B, [t1,u1])*wp, Evaluate(X1A, [t1,u1]) + Evaluate(X1B, [t1,u1])*wp, Evaluate(X2A, [t1,u1]) + Evaluate(X2B, [t1,u1])*wp];
  if P[1] eq 0 then return -1; end if;
  P := [K | 1, P[2]/P[1], P[3]/P[1]];
  val := [Evaluate(f, P) : f in forms]; vx := [Evaluate(f, P) : f in d1]; vy := [Evaluate(f, P) : f in d2];
  if val[3] eq 0 or val[4] eq 0 then return -1; end if;
  w := (t1*val[4] - val[2])/val[3];                         // from t = (At + Bt w)/Gb
  if w^2 ne val[1] or w eq 0 then return 0; end if;         // P is a point of X
  V := val[7] + val[8]*w; U := val[5] + val[6]*w;
  if V eq 0 or U ne u1*V then return 0; end if;             // u(P) = u1
  T := val[2] + val[3]*w; G := val[4];
  dxx := func<i | vx[i] + vx[i+1]*w + val[i+1]*vx[1]/(2*w)>; dyy := func<i | vy[i] + vy[i+1]*w + val[i+1]*vy[1]/(2*w)>;
  tx := (dxx(2)*G - T*vx[4])/G^2; ty := (dyy(2)*G - T*vy[4])/G^2;
  ux := (dxx(5)*V - U*dxx(7))/V^2; uy := (dyy(5)*V - U*dyy(7))/V^2;
  return w*(tx*uy - ty*ux)/8 eq wp select 1 else 0;         // w'(P) = wp
end function;
npts := 0; nbad := 0; minpts := 1000;
for t1 in [lo+1..hi+1] do
  nf := 0; u1 := 0;                                           // nf: points of this fibre at which psi = phi^-1 is proved
  while nf lt 102 do
    u1 +:= 1;
    m := t1*(u1^3*(t1+1) - 3*a*u1^2*t1 - t1^2);
    if m eq 0 then continue; end if;
    fl, rt := IsSquare(Q ! m);
    if fl then
      for wp in [rt, -rt] do r := TestPoint(Q, t1, u1, wp); if r eq 1 then nf +:= 1; npts +:= 1; elif r eq 0 then nbad +:= 1; end if; end for;
    else
      D, f := SquarefreeFactorization(Integers() ! m); K := QuadraticField(D); wp := f*K.1; assert wp^2 eq m;
      r := TestPoint(K, t1, u1, wp);                          // the conjugate point follows by conjugation
      if r eq 1 then nf +:= 2; npts +:= 2; elif r eq 0 then nbad +:= 1; end if;
    end if;
    if u1 gt 200 then break; end if;
  end while;
  minpts := Min(minpts, nf);
end for;
printf "fibres %o..%o: %o points, at least %o on each fibre, %o failures: %o\n", lo, hi, npts, minpts, nbad, (nbad eq 0 and minpts ge 101) select "ok" else "FAILED";
quit;

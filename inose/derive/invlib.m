// inverse map by interpolation: forms X_j = A'_j(t,u) + B'_j(t,u) w' of degree KK (and KK-3) on the Inose double plane with
// (X_0 : X_1 : X_2) = (x0 : x1 : x2).  Everything modulo a prime p.
Q := Rationals();
R<x0,x1,x2> := PolynomialRing(Q, 3);
S := 3686400*x0^6 + (-8697600*x1^2 + 2193664*x2*x1 + 11520*x2^2)*x0^4 + (5130225*x1^4 - 1860100*x2*x1^3 + 138414*x2^2*x1^2 - 4180*x2^3*x1 + 9*x2^4)*x0^2 + (-45144*x2*x1^5 + 19960*x2^2*x1^4 - 1672*x2^3*x1^3 + 40*x2^4*x1^2);
a := 53360;
load "map_data.m";
la := [ <9*x0 - x1, 1>, <80*x0 - 11*x1 + x2, 2>, <72*x0 - 27*x1 + x2, 3>, <112*x0 - 49*x1 + 3*x2, 4>, <33*x0 - 22*x1 + x2, 5>,
        <80*x0 - 57*x1 + 3*x2, 6>, <108*x0 - 99*x1 + 5*x2, 4>, <768*x0^2 - 1960*x0*x1 + 100*x0*x2 + 1152*x1^2 - 117*x1*x2 + 3*x2^2, 2>, <24*x0 - 19*x1 + x2, 3> ];
Gb := Evaluate(&*[f[1]^f[2] : f in la], [x0, -x1, -x2]);
forms := [S, At, Bt, Gb, A0, B0, Ai, Bi];
D1 := func<f | Derivative(f, 2)>; D2 := func<f | Derivative(f, 3)>;
// a random point of X over GF(p) and its image (t, u, w')
function ImagePoint(Fp, formsp, d1p, d2p)
  repeat
    x := Random(Fp); y := Random(Fp); P := [Fp | 1, x, y];
    s := Evaluate(formsp[1], P); fl, w := IsSquare(s);
    if fl and s ne 0 then
      val := [Evaluate(f, P) : f in formsp]; vx := [Evaluate(f, P) : f in d1p]; vy := [Evaluate(f, P) : f in d2p];
      pr := func<i | val[i] + val[i+1]*w>;
      dxx := func<i | vx[i] + vx[i+1]*w + val[i+1]*vx[1]/(2*w)>; dyy := func<i | vy[i] + vy[i+1]*w + val[i+1]*vy[1]/(2*w)>;
      T := pr(2); U := pr(5); V := pr(7); G := val[4];
      if G ne 0 and V ne 0 and T ne 0 then
        t := T/G; u := U/V;
        tx := (dxx(2)*G - T*vx[4])/G^2; ty := (dyy(2)*G - T*vy[4])/G^2;
        ux := (dxx(5)*V - U*dxx(7))/V^2; uy := (dyy(5)*V - U*dyy(7))/V^2;
        wp := w*(tx*uy - ty*ux)/8;
        assert wp^2 eq t*(u^3*(t+1) - 3*a*u^2*t - t^2);
        return x, y, w, t, u, wp;
      end if;
    end if;
  until false;
end function;
function Exps2(k) return [[i, j] : j in [0..k-i], i in [0..k]]; end function;
function InvKernel(p, KK, npts : Seed := 1)
  SetSeed(Seed);
  Fp := GF(p); Rp := ChangeRing(R, Fp);
  formsp := [Rp ! f : f in forms]; d1p := [Rp ! D1(f) : f in forms]; d2p := [Rp ! D2(f) : f in forms];
  EA := Exps2(KK); EB := Exps2(KK - 3); n := #EA + #EB;
  rows := [];
  for it in [1..npts] do
    x, y, w, t, u, wp := ImagePoint(Fp, formsp, d1p, d2p);
    tp := [Fp ! 1]; up := [Fp ! 1]; for i in [1..KK] do Append(~tp, tp[i]*t); Append(~up, up[i]*u); end for;
    m := [tp[e[1]+1]*up[e[2]+1] : e in EA] cat [tp[e[1]+1]*up[e[2]+1]*wp : e in EB];
    z := [Fp ! 0 : i in [1..n]];
    Append(~rows, [x*c : c in m] cat [-c : c in m] cat z);
    Append(~rows, [y*c : c in m] cat z cat [-c : c in m]);
  end for;
  M := Matrix(Fp, rows);
  return Kernel(Transpose(M)), EA, EB;
end function;

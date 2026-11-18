// pencil.m : sections N = A + B*w of k*D on Elkies' sextic (w^2 = S) vanishing along given lifted curves to given orders.
Q := Rationals();
R<x0,x1,x2> := PolynomialRing(Q, 3);
S := 3686400*x0^6 + (-8697600*x1^2 + 2193664*x2*x1 + 11520*x2^2)*x0^4 + (5130225*x1^4 - 1860100*x2*x1^3 + 138414*x2^2*x1^2 - 4180*x2^3*x1 + 9*x2^4)*x0^2 + (-45144*x2*x1^5 + 19960*x2^2*x1^4 - 1672*x2^3*x1^3 + 40*x2^4*x1^2);
PV<v0,v1> := PolynomialRing(Q, 2);
PS<s> := PolynomialRing(Q);
load "elkdata.m";
// a sample point of a lifted curve: <point P0 in Q^3, value W0 of w, transversal vector n>
function LineSamples(i, sgn, npts)
  l := LIN[i]; p := l[2]; q := l[3]; cu := l[4]; res := [];
  n := [0,0,0]; j := 1; while l[1][j] eq 0 do j +:= 1; end while; n[j] := 1;
  a := 0;
  while #res lt npts do
    a +:= 1; pr := [[1, a], [1, -a]];
    for vv in pr do
      P0 := [p[m]*vv[1] + q[m]*vv[2] : m in [1..3]]; W0 := sgn*Evaluate(cu, vv);
      if W0 ne 0 then Append(~res, <P0, W0, n>); end if;
    end for;
  end while;
  return res;
end function;
function ConicSamples(i, sgn, npts)
  X := CONX[i]; g := CONG[i]; res := []; a := 0;
  while #res lt npts do
    a +:= 1;
    for s0 in [a, -a] do
      P0 := [Evaluate(X[m], s0) : m in [1..3]]; W0 := sgn*Evaluate(g, s0);
      T := [Evaluate(Derivative(X[m]), s0) : m in [1..3]];
      // transversal: a coordinate vector not in the span of P0 and the tangent T
      n := [0,0,0];
      for j in [1..3] do nn := [0,0,0]; nn[j] := 1; if Determinant(Matrix(Q, 3, 3, P0 cat T cat nn)) ne 0 then n := nn; break; end if; end for;
      if W0 ne 0 and n ne [0,0,0] then Append(~res, <P0, W0, n>); end if;
    end for;
  end while;
  return res;
end function;
// conditions: for each sample, the coefficients of e^0..e^(c-1) of N(P0 + e n, W(e))
function Conditions(k, samples, c, MA, MB : Field := Q)
  K := Field;
  Pe<e> := PolynomialRing(K); Te := quo<Pe | e^c>;
  rows := [];
  SK := ChangeRing(S, K);
  for smp in samples do
    X := [Te ! (K ! smp[1][m] + (K ! smp[3][m])*e) : m in [1..3]]; W0 := K ! smp[2];
    SX := Evaluate(SK, X);
    del := SX/(W0^2) - 1;                       // in e*Te
    // sqrt(1 + del) = sum binom(1/2, j) del^j
    sq := Te ! 1; pw := Te ! 1; bc := K ! 1;
    for j in [1..c-1] do bc := bc*(K ! (3 - 2*j))/(K ! (2*j)); pw := pw*del; sq := sq + bc*pw; end for;
    W := W0*sq;
    pws := [[Te ! 1] : m in [1..3]];
    for m in [1..3] do for j in [1..k] do Append(~pws[m], pws[m][j]*X[m]); end for; end for;
    vals := [pws[1][ex[1]+1]*pws[2][ex[2]+1]*pws[3][ex[3]+1] : ex in MA] cat [pws[1][ex[1]+1]*pws[2][ex[2]+1]*pws[3][ex[3]+1]*W : ex in MB];
    cf := [Eltseq(Pe ! v) cat [K ! 0 : jj in [1..c]] : v in vals];
    for j in [1..c] do Append(~rows, [cf[t][j] : t in [1..#vals]]); end for;
  end for;
  return rows;
end function;
function Exps(k) return [[a, b, k-a-b] : b in [0..k-a], a in [0..k]]; end function;
// Zlist: list of <"L" or "C", index, sign, multiplicity>; emult = multiplicity along the node curve E (node (0:0:1))
function PencilRows(k, Zlist, emult : Field := Q)
  MA := [ex : ex in Exps(k) | ex[1] + ex[2] ge emult];          // mult_node(A) >= emult
  MB := [ex : ex in Exps(k-3) | ex[1] + ex[2] ge emult - 1];    // mult_node(B) >= emult - 1
  rows := [];
  for z in Zlist do
    c := z[4];
    if z[1] eq "L" then smp := LineSamples(z[2], z[3], k + 6*c + 3); else smp := ConicSamples(z[2], z[3], 2*(k + 6*c) + 3); end if;
    rows cat:= Conditions(k, smp, c, MA, MB : Field := Field);
  end for;
  return rows, MA, MB;
end function;

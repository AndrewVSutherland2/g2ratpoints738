/* interlib.gp : intersection numbers on X (the minimal resolution of Elkies' double plane) between the lifts of the conics
   and lines of Elkies' model (needs cplib.gp and out/lifts.gp).
   Lift C_k of conic k: {q_k = 0, w = +G_k(x)}, C_k' = {q_k = 0, w = -G_k(x)}; lift l+ of line l: {l = 0, w = +G_l(x)},
   l- = {l = 0, w = -G_l(x)}; for a curve through the node N = (0:0:1) the lift means the strict transform on X (for a node
   line, l+ and l- are the two components E, E' of its pullback E + E' + r; r = exceptional curve of the node, r^2 = -2).
   1. Curve 1 not through N: the plane intersection points with curve 2 are away from N and
        C_1.C_2 = deg gcd( F_2(x_1(s,t)), g_1 - G_2(x_1) ),  C_1.C_2' = deg gcd( F_2(x_1), g_1 + G_2(x_1) )
      as binary forms (root (1:0) included); exact unless an intersection point lies on the branch sextic (never happened;
      the third value returned is the multiplicity of such points).
   2. Curve 1 through N, curve 2 not: by the covering involution C_1.C_2 = C_2.C_1 and C_1.C_2' = C_2'.C_1 = C_2.C_1', i.e.
      the same two numbers computed from curve 2.
   3. Both through N (exceptional-curve rule): the points away from N are counted as in 1 (the parameter of N is divided out
      of F_2(x_1)), and N contributes (m - 1) (m = plane intersection multiplicity at N) to C_1.C_2 if the two lifts pass
      through the same point of r, to C_1.C_2' if C_1 and C_2' do: the point of r of a lift is (alpha : beta : gamma) with
      (alpha, beta) the tangent direction at N in the chart x2 = 1 and gamma the derivative of w/x2^3, on the conic
      r = {gamma^2 = T(alpha, beta)}, T = tangent cone of S at N (T = 9 a^2 + 40 b^2 for Elkies' sextic).
      Consistency: C_1.C_2 + C_1.C_2' = (2h - r).C_1 = 3 for two conics through N, (h - r).C = 1 for a conic and a line
      through N, 0 for two node lines (checked by checkN()).
   C.r = 1 for a curve through N, 0 otherwise. */
read(liftsfile());
NC = #CON; NL = #LIN;
isnodeconic(k) = CON[k][1][6] == 0;
haspar(k) = CON[k][2] != 0;
\\ tangent cone of S at N (coefficients of x0^2 x2^4, x0 x1 x2^4, x1^2 x2^4)
TC = [cE[15], cE[20], cE[26]];
Tcone(al, be) = TC[1]*al^2 + TC[2]*al*be + TC[3]*be^2;
if(Sx(NODE) != 0 || gradat(SE, NODE) != [0,0,0], error("interlib: (0:0:1) is not a singular point of S"));
\\ curve records [xv, gv, F, G, e, thruN]: parametrisation xv (3 binary forms of degree e in s,t), lift value gv = w on the
\\ lift (binary form of degree 3e), plane equation F (form in x), cubic lift G (form in x), flag "through N"
CV = vector(NC, k, if(!haspar(k), 0, [cpar(CON[k][2]), CON[k][4], qpol(CON[k][1]), cubic(CON[k][5]), 2, isnodeconic(k)]));
LV = vector(NL, l, [lpar(LIN[l][1]), LIN[l][3], LIN[l][1][1]*x0 + LIN[l][1][2]*x1 + LIN[l][1][3]*x2, cubic(LIN[l][4]), 1, LIN[l][1][3] == 0]);
\\ fast path data (dehomogenised, t = 1) for the conics
X1 = vector(NC); MV = vector(NC); G1 = vector(NC);
{
for(k=1,NC, if(!haspar(k), next); my(xv = subst(cpar(CON[k][2]), t, 1));
  X1[k] = xv; MV[k] = vector(10, n, monval(MON3[n], xv)); G1[k] = subst(CON[k][4], t, 1));
}
\\ degree of the gcd of dehomogenised binary forms A (degree dA), B (degree dB), the root (1:0) included; B = 0 allowed
gdeg(A, dA, B, dB) = { if(B == 0, return(dA)); my(d = poldegree(gcd(A, B), s)); if(d < 0, d = 0); d + min(dA - poldegree(A, s), dB - poldegree(B, s)) }
\\ [C_i.C_j, C_i.C_j', branch multiplicity] for conics, conic i not through the node (fast path)
interCC(i, j) = { my(qj = qform(CON[j][1], X1[i]), Gj = sum(n=1,10, CON[j][5][n]*MV[i][n]));
  if(isnodeconic(i), error("interCC: conic ", i, " passes through the node"));
  if(qj == 0, error("conics ", i, " and ", j, " coincide"));
  [gdeg(qj, 4, G1[i] - Gj, 6), gdeg(qj, 4, G1[i] + Gj, 6), gdeg(qj, 4, G1[i], 6)] }
\\ [C_i.l+, C_i.l-, branch multiplicity], conic i not through the node, l = index in LIN
interCL(i, l) = { my(lv = LIN[l][1]*X1[i]~, Gl = sum(n=1,10, LIN[l][4][n]*MV[i][n]));
  if(isnodeconic(i), error("interCL: conic ", i, " passes through the node"));
  if(lv == 0, error("interCL: zero form"));
  [gdeg(lv, 2, G1[i] - Gl, 6), gdeg(lv, 2, G1[i] + Gl, 6), gdeg(lv, 2, G1[i], 6)] }
\\ ---- general intersection of two curve records
\\ dehomogenised data of a curve, in the chart (s, t=1) or, if N is at the parameter (1:0), in the swapped chart:
\\ [x (3 polynomials in s), g (polynomial in s), s0 (parameter of N, or 0 if the curve does not pass through N)]
cdata(cv) = { my(xv = cv[1], gv = cv[2], x, g, s0 = 0);
  if(cv[6], my(h = gcd(subst(xv[1], t, 1), subst(xv[2], t, 1)));
    if(poldegree(h, s) < 1, xv = substvec(xv, [s, t], [t, s]); gv = substvec(gv, [s, t], [t, s]); h = gcd(subst(xv[1], t, 1), subst(xv[2], t, 1)));
    if(poldegree(h, s) != 1, error("cdata: the curve does not pass through N exactly once"));
    s0 = -polcoeff(h, 0, s)/polcoeff(h, 1, s));
  [subst(xv, t, 1), subst(gv, t, 1), s0] }
\\ point (alpha : beta : gamma) of r through which the lift of a curve through N passes
rpoint(cv) = { my(d = cdata(cv), x = d[1], g = d[2], s0 = d[3], aa = x[1]/x[3], bb = x[2]/x[3], yy = g/x[3]^3, v);
  if(!cv[6], error("rpoint: curve not through N"));
  v = [subst(deriv(aa, s), s, s0), subst(deriv(bb, s), s, s0), subst(deriv(yy, s), s, s0)];
  if(v[1] == 0 && v[2] == 0, error("rpoint: singular parametrisation at N"));
  if(v[3]^2 != Tcone(v[1], v[2]), error("rpoint: gamma^2 != T(alpha, beta)"));
  if(v[3] == 0, error("rpoint: the curve is tangent to a branch of the node (gamma = 0): rule not applicable"));
  prim(v) }
rconj(p) = [p[1], p[2], -p[3]];
\\ [C_1.C_2, C_1.C_2', branch multiplicity away from N]
inter(c1, c2) = { my(d1, xx, gg, s0, F2x, G2x, D, m = 0, e1 = c1[5], q);
  if(c1[6] && !c2[6], return(inter(c2, c1)));
  d1 = cdata(c1); xx = d1[1]; gg = d1[2]; s0 = d1[3];
  F2x = substvec(c2[3], [x0,x1,x2], xx); G2x = substvec(c2[4], [x0,x1,x2], xx);
  if(F2x == 0, error("inter: the curves coincide"));
  D = e1 * c2[5];                       \\ degree of F_2(x_1(s,t)): (degree of F_2 = 2 for a conic, 1 for a line) x e_1
  q = F2x;
  if(c1[6], while(subst(q, s, s0) == 0, q = q/(s - s0); m++));
  my(dp = gdeg(q, D - m, gg - G2x, 3*e1), dm = gdeg(q, D - m, gg + G2x, 3*e1), br = gdeg(q, D - m, gg, 3*e1));
  if(c1[6] && m > 1, my(p1 = rpoint(c1), p2 = rpoint(c2)); if(p1 == p2, dp += m - 1); if(p1 == rconj(p2), dm += m - 1));
  [dp, dm, br] }
selfint(cv) = if(cv[6], [-2, 5], [-2, 6]);           \\ [C.C, C.C'] ((2h - r).C = 3 or 2h.C = 4)
rint(cv) = if(cv[6], 1, 0);                           \\ C.r
\\ consistency check of the exceptional-curve rule on all pairs of curves through N: the sums C_1.C_2 + C_1.C_2'
checkN() = { my(bad = 0, cnt = 0, nodeC = select(k -> haspar(k) && isnodeconic(k), vector(NC, k, k)), nodeL = select(l -> LV[l][6], vector(NL, l, l)));
  for(a=1,#nodeC, for(b=a+1,#nodeC, my(r = inter(CV[nodeC[a]], CV[nodeC[b]])); cnt++; if(r[1] + r[2] != 3, bad++; if(bad <= 5, print("conics ", nodeC[a], ",", nodeC[b], ": ", r)))));
  for(a=1,#nodeC, for(b=1,#nodeL, my(r = inter(CV[nodeC[a]], LV[nodeL[b]])); cnt++; if(r[1] + r[2] != 1, bad++; if(bad <= 5, print("conic ", nodeC[a], ", line ", nodeL[b], ": ", r)))));
  for(a=1,#nodeL, for(b=a+1,#nodeL, my(r = inter(LV[nodeL[a]], LV[nodeL[b]])); cnt++; if(r[1] + r[2] != 0, bad++; if(bad <= 5, print("lines ", nodeL[a], ",", nodeL[b], ": ", r)))));
  print("exceptional-curve rule: ", cnt, " pairs of curves through the node checked, ", bad, " violations of C_1.C_2 + C_1.C_2' = (2h - r).C_1 - [self]");
  bad == 0 }
\\ mirror x0 -> -x0 (a symmetry of Elkies' model): index of the mirror conic / line
mirrorq(q) = [q[1], -q[2], -q[3], q[4], q[5], q[6]];
MIRC = vector(NC, k, my(mq = mirrorq(CON[k][1]), r = 0); for(j=1,NC, if(CON[j][1] == mq, r = j; break)); r);
primline(l) = { my(v = l/content(l)); if(v[1] < 0 || (v[1] == 0 && v[2] < 0) || (v[1] == 0 && v[2] == 0 && v[3] < 0), v = -v); v }
MIRL = vector(NL, k, my(ml = primline([-LIN[k][1][1], LIN[k][1][2], LIN[k][1][3]]), r = 0); for(j=1,NL, if(LIN[j][1] == ml, r = j; break)); r);
\\ degree profile of D' = C_i + sign*C_j (sign = +1: C_j, -1: C_j') on all known curves: returns
\\ [degC (vector NC: D'.C_k), degCp (D'.C_k'), degL (D'.l+), degLp (D'.l-), degr (D'.r)]; -1 = not computable
\\ (with the lift C~ = C_j or C_j': [C~.C, C~.C'] = inter(c_j, C) for sign +1 and its reverse for sign -1, by the involution)
lifted(cj, sg, C) = { my(r = inter(cj, C)); if(sg == 1, [r[1], r[2]], [r[2], r[1]]) }
profile(i, j, sg) = { my(dC = vector(NC), dCp = vector(NC), dL = vector(NL), dLp = vector(NL), ci = CV[i], cj = CV[j]);
  for(k=1,NC, if(!haspar(k), dC[k] = dCp[k] = -1; next);
    if(k == i || k == j, my(o = if(k == i, j, i), so = if(k == i, sg, 1), sf = selfint(CV[k]), r = lifted(CV[o], if(k == i, sg, 1), CV[k]));
      \\ D'.C_k = C_k^2 + (lift of the other).C_k ; for k = j with sign -1 the chosen lift is C_j'
      if(k == i, dC[k] = sf[1] + r[1]; dCp[k] = sf[2] + r[2],
        my(rr = inter(ci, cj)); if(sg == 1, dC[k] = sf[1] + rr[1]; dCp[k] = sf[2] + rr[2], dC[k] = sf[2] + rr[1]; dCp[k] = sf[1] + rr[2]));
      next);
    my(ri = inter(ci, CV[k]), rj = lifted(cj, sg, CV[k]));
    dC[k] = ri[1] + rj[1]; dCp[k] = ri[2] + rj[2]);
  for(l=1,NL, my(ri = inter(ci, LV[l]), rj = lifted(cj, sg, LV[l]));
    dL[l] = ri[1] + rj[1]; dLp[l] = ri[2] + rj[2]);
  [dC, dCp, dL, dLp, rint(ci) + rint(cj)] }

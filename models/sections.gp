/* sections.gp : basis of L(D') for D' = C_i + sign*C_j (sign = +1: D' = C_i + C_j with C_i.C_j = 3; -1: D' = C_i + C_j'
   with C_i.C_j' = 3).  H^0(X, O(D')) = { f = (P(x) + Q(x) w) / (q_i q_j) : P quartic form, Q linear form, P + Q w vanishing
   on the other lifts C_i' = {q_i = 0, w = -G_i} and C_j~' = {q_j = 0, w = -sign G_j} }, i.e. P - Q G_i = 0 mod q_i and
   P - sign Q G_j = 0 mod q_j: 9 + 9 linear conditions on 15 + 3 unknowns, of rank 15 (the 3 common points of C_i' and
   C_j~' give shared conditions), kernel of dimension 3 = h^0(D') (D' nef and big, D'^2 = 2).  Conics through the node:
   the divisor of q_k on X is C_k + C_k' + r, so the numerator must also vanish on r to the order (number of node conics);
   for one node conic this is P(N) = 0, already implied; for two, P must vanish to order 2 at N and Q(N) = 0 (4 more rows).
   sections(i, j, sg) returns [[P0,Q0],[P1,Q1],[P2,Q2]] (integral, LLL-reduced) in the convention of modellib.gp
   (function (A + y B)/x2^N with N = 4; the common denominator q_i q_j is dropped: it only rescales the map).
   With PAIR = [i, j, sign] set (./run.sh sections.gp logs/sections_i_j.log 'PAIR=[i,j,sign]') the basis is computed,
   checked and written to out/sections_i_j.gp. */
\r cplib.gp
\r interlib.gp
MON4 = mons(4); MON1 = mons(1);
polcoefs(P, M) = vector(#M, n, polcoeff(polcoeff(polcoeff(P, M[n][1], x0), M[n][2], x1), M[n][3], x2));   \\ coefficients of a form on the monomial list M
row8(F) = vector(9, n, polcoeff(polcoeff(F, 9-n, s), n-1, t));   \\ the 9 coefficients of a binary form of degree 8
condrows(xv, g) = { my(cols = vector(18)); for(n=1,15, cols[n] = row8(monval(MON4[n], xv))); for(n=1,3, cols[15+n] = row8(-g*monval(MON1[n], xv))); matconcat(cols~)~ }
sections(i, j, sg) = { my(xi = cpar(CON[i][2]), xj = cpar(CON[j][2]), gi = CON[i][4], gj = sg*CON[j][4], M, K, f, qq, v);
  M = matconcat([condrows(xi, gi); condrows(xj, gj)]);
  if(matsize(M) != [18, 18], error("sections: matrix size"));
  \\ conics through the node: div(q_k) = C_k + C_k' + r on X, so P + Qw must vanish on r to the order (number of node
  \\ conics); w vanishes on r to order 1: one node conic needs P(N) = 0 (implied by the conditions on C_k' at the parameter
  \\ of N), two need P vanishing to order 2 at N and Q(N) = 0: 4 more rows (only some of them new)
  if(isnodeconic(i) && isnodeconic(j),
    my(rows = List()); for(n=1,15, if(MON4[n][3] >= 3, my(row = vector(18)); row[n] = 1; listput(rows, row)));
    for(n=1,3, if(MON1[n][3] == 1, my(row = vector(18)); row[15+n] = 1; listput(rows, row)));
    M = matconcat([M; matconcat(Vec(rows)~)]));
  print("condition matrix: rank ", matrank(M), " (expected 15)", if(isnodeconic(i) || isnodeconic(j), " [conic(s) through the node]", ""));
  K = matker(M, 1); if(#K == 0, error("sections: empty kernel"));
  K = K*qflll(K);
  if(#K != 3, error("sections: dim L(D') = ", #K, ", not 3"));
  f = vector(3, k, [sum(n=1,15, K[n,k]*monval(MON4[n], [x0,x1,x2])), sum(n=1,3, K[15+n,k]*monval(MON1[n], [x0,x1,x2]))]);
  \\ the constant function 1 = (q_i q_j + 0 w)/(q_i q_j) must be in the space
  qq = qpol(CON[i][1])*qpol(CON[j][1]);
  v = concat(vector(15, n, polcoeff(polcoeff(polcoeff(qq, MON4[n][1], x0), MON4[n][2], x1), MON4[n][3], x2)), [0,0,0])~;
  if(M*v != 0, error("sections: q_i q_j does not satisfy the conditions"));
  if(matrank(matconcat([K, v])) != 3, error("sections: q_i q_j is not in the span of the basis"));
  \\ direct re-check of the vanishing conditions on the basis
  for(k=1,3, if(substvec(f[k][1], [x0,x1,x2], xi) - gi*substvec(f[k][2], [x0,x1,x2], xi) != 0, error("basis element ", k, " does not vanish on C_i'"));
    if(substvec(f[k][1], [x0,x1,x2], xj) - gj*substvec(f[k][2], [x0,x1,x2], xj) != 0, error("basis element ", k, " does not vanish on C_j~'")));
  f }
{
if(type(PAIR) == "t_VEC",
  my(i = PAIR[1], j = PAIR[2], sg = PAIR[3], t0 = getabstime(), f);
  \\ sign 0 = determine from the exact intersection numbers (as model.gp does; needed because this block runs before model.gp
  \\ gets to normalise the sign -- with sg = 0 the conic-j condition degenerated to P = 0 mod q_j and dim L came out 1)
  if(sg == 0, my(rr = inter(CV[i], CV[j])); sg = if(rr[1] == 3, 1, if(rr[2] == 3, -1, error("sections: ", [i,j], " is not a 3-pair: ", rr))); PAIR = [i, j, sg]);
  print("pair ", PAIR, ": conic ", i, " = ", CON[i][1], ", conic ", j, " = ", CON[j][1], ", C_i.C_j = ", inter(CV[i], CV[j])[1], ", lift of conic j: ", if(sg == 1, "+G_j", "-G_j"));
  f = sections(i, j, sg);
  print("basis of L(D') (P, Q):"); for(k=1,3, print("  ", f[k]));
  print("max |coefficient|: ", vecmax(abs(concat(vector(3, k, concat(polcoefs(f[k][1], MON4), polcoefs(f[k][2], MON1)))))), "  (", (getabstime()-t0)/1000., " s)");
  writefile(concat(["out/sections_", Str(i), "_", Str(j), ".gp"]), Str("\\\\ basis of L(C_i + sign C_j) for PAIR = ", PAIR, " (models/sections.gp): functions (P + w Q)/(q_i q_j)\nPAIR = ", PAIR, ";\nf = ", f, ";"));
  print("written out/sections_", i, "_", j, ".gp"));
}

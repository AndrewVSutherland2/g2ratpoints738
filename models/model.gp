/* model.gp : raw double-plane model of X_163 for the polarization D' = C_i + sign C_j (PAIR = [i, j, sign], NAME = "...").
   1. basis f = (f0, f1, f2) of L(D') (sections.gp), numerators N_m = P_m + w Q_m (the common denominator q_i q_j is dropped);
   2. w' = y det(N, d_a N, d_b N) (modellib.gp: anti-invariant for the new covering
      involution, w'^2 = S_new(N0, N1, N2) for a unique sextic form S_new), fitted exactly from rational points on other
      conics and checked on fresh points;
   3. images of the known curves of X (conic lifts, line lifts, node curve r) of new degree 1 (lines of the new model: S_new is
      a square on them) and of new degree 0 (contracted to singular points of S_new);
   Output: out/model_<NAME>_raw.gp (PAIR, f, c (28 coefficients, linesearch order), S, newlines (with their origin), newpoints
   (images of the contracted known curves), nodeimg (image of r)).
   usage: ./run.sh model.gp logs/model_<NAME>.log 'PAIR=[i,j,sign]' 'NAME="<NAME>"' */
\r cplib.gp
\r interlib.gp
\r sections.gp
\r modellib.gp
if(type(PAIR) != "t_VEC", error("set PAIR = [i, j, sign]"));
if(type(NAME) != "t_STR", NAME = concat(["p", Str(PAIR[1]), "_", Str(PAIR[2])]));
i = PAIR[1]; j = PAIR[2]; sg = PAIR[3];
T0 = getabstime();
\\ the sign is determined by the exact intersection numbers: C_i.C_j = 3 -> +1 (D' = C_i + C_j), C_i.C_j' = 3 -> -1 (D' = C_i + C_j');
\\ sign 0 = choose automatically; a given sign is checked (the lattice file's eps uses another sign convention for sqrt)
if(isnodeconic(i) && !isnodeconic(j), my(tmp = i); i = j; j = tmp; PAIR = [i, j, sg]; print("NOTE: conics swapped so that conic i does not pass through the node"));
r0 = inter(CV[i], CV[j]);
sgx = if(r0[1] == 3, 1, if(r0[2] == 3, -1, error("model: C_i.C_j = ", r0[1], ", C_i.C_j' = ", r0[2], ": not a 3-pair")));
if(sg == 0, sg = sgx; PAIR = [i, j, sg]);
if(sg != sgx, print("NOTE: with the sign conventions of lifts.gp the 3-pair is (", i, ", ", j, ", ", sgx, "), not ", PAIR, "; using ", sgx); sg = sgx; PAIR = [i, j, sg]);
print("C_i.C_j = ", r0[1], ", C_i.C_j' = ", r0[2], " -> D' = C_", i, " + ", if(sg == 1, "C_", "C'_"), j);
f = sections(i, j, sg);
qij = qpol(CON[i][1])*qpol(CON[j][1]);
\\ ---- 2. the branch sextic
w = wsym(f);
print("w' = y det(N, N_a, N_b) computed (", (getabstime()-T0)/1000., " s)");
\\ rational points of X on conic k (parameter v, both signs of w), as [x, w]; skipped: x2 = 0, w = 0, points on q_i q_j = 0
conicpts(k, vs) = { my(L = List(), xv = X1[k], gv = G1[k]); for(n=1,#vs, my(x = substvec(xv, [s], [vs[n]]), y = subst(gv, s, vs[n]));
    if(x[3] == 0 || y == 0 || substvec(qij, [x0,x1,x2], x) == 0, next); listput(L, [x, y]); listput(L, [x, -y])); Vec(L) }
evalpts(ks, vs) = { my(L = List()); for(n=1,#ks, my(k = ks[n]); if(k == i || k == j || !haspar(k), next); my(P = conicpts(k, vs)); for(m=1,#P, listput(L, affpt(P[m])))); Vec(L) }
KS = select(k -> k != i && k != j && haspar(k), vector(60, n, n))[1..40];
VS = [0, 1, -1, 2, -2, 1/2, -1/2, 3, -3];
pts = evalpts(KS[1..20], VS);
print(#pts, " evaluation points on ", 20, " conics");
c0 = fitsextic(f, w, pts);
c = normsq(c0); S = sextic(c);
print("S_new fitted: max |coeff| has ", #Str(vecmax(abs(c))), " digits, ", 28 - #select(x -> x != 0, c), " zero coefficients (", (getabstime()-T0)/1000., " s)");
\\ consistency at fresh points: S_new(N(P)) / w'(P)^2 must be one constant (the square scaling absorbed by normsq)
g = vector(3, m, aff(f[m]));
pts2 = evalpts(KS[21..40], [5, -5, 1/3, -1/3, 7]);
\\ w' vanishes at the points of the contracted curves (the map is constant there) and on the new ramification curve; at such
\\ a fresh point S_new(N) must vanish too (the old check divided by w'^2 there, and gp went on after
\\ the error in 45 of the 99 runs, so that the consistency check had silently not been evaluated for them)
rat = List(); nz = 0;
{
for(n=1,#pts2, my(pt = pts2[n], fv = vector(3, m, geval(g[m], pt)), wv = geval(w, pt), sv = substvec(S, [x0,x1,x2], fv));
  if(wv == 0, nz++; if(sv != 0, error("model: w' = 0 but S_new(N) != 0 at the fresh point ", pt)); next);
  listput(rat, sv/wv^2));
}
rat = Set(rat);
print("S_new(N)/w'^2 at ", #pts2 - nz, " fresh points: ", rat, " (", nz, " fresh points with w' = 0 and S_new(N) = 0: on a contracted curve or on the new ramification curve)", if(#rat == 1 && #pts2 - nz >= 50, "  OK", "  INCONSISTENT"));
if(#rat != 1 || #pts2 - nz < 50, error("model: inconsistent fit"));
\\ ---- 3. images of the known curves
Nmap(P) = vector(3, m, substvec(f[m][1], [x0,x1,x2], P[1]) + P[2]*substvec(f[m][2], [x0,x1,x2], P[1]));
\\ points of a known curve cv = ["C", k, e] or ["L", l, e] at parameters vs (homogeneous [x, w])
curvepts(cv, vs) = { my(L = List());
  if(cv[1] == "C", my(k = cv[2]); for(n=1,#vs, my(x = substvec(X1[k], [s], [vs[n]]), y = cv[3]*subst(G1[k], s, vs[n])); listput(L, [x, y])),
    my(l = cv[2], PQ = LIN[l][2]); for(n=1,#vs, my(x = vs[n]*PQ[1] + PQ[2], y = cv[3]*substvec(LIN[l][3], [t, s], [vs[n], 1])); listput(L, [x, y])));
  Vec(L) }
\\ image of a known curve of degree d (0 or 1): a point or a line, from its image points (zero images skipped)
imgcurve(cv, d) = { my(P = curvepts(cv, [0, 1, -1, 2, -2, 1/2, -1/2, 3, -3, 1/3, -1/3, 4, -4, 5, -5]), Q = List());
  for(n=1,#P, my(q = Nmap(P[n])); if(q != [0,0,0], listput(Q, prim(q))));
  if(d == 0, if(#Set(Q) != 1, error("image of the contracted curve ", cv, " is not a point: ", Set(Q))); return(Q[1]));
  if(#Q < 3, error("too few image points for ", cv));
  my(l = prim(cross(Q[1], Q[2]))); for(n=3,#Q, if(l*Q[n]~ != 0, error("image of ", cv, " is not a line")));
  l }
pr = profile(i, j, sg); dC = pr[1]; dCp = pr[2]; dL = pr[3]; dLp = pr[4]; degr = pr[5];
known = List();
{
for(k=1,NC, if(dC[k] < 0, next); listput(known, [["C", k, 1], dC[k]]); listput(known, [["C", k, -1], dCp[k]]));
for(l=1,NL, listput(known, [["L", l, 1], dL[l]]); listput(known, [["L", l, -1], dLp[l]]));
}
known = Vec(known);
newlines = List(); newpoints = List();
{
for(n=1,#known, my(cv = known[n][1], d = known[n][2]);
  if(d == 1, listput(newlines, [cv, imgcurve(cv, 1)]));
  if(d == 0, listput(newpoints, [cv, imgcurve(cv, 0)])));
}
\\ the exceptional curve r (D'.r = degr = number of the two conics through the node): its rational points (alpha:beta:gamma),
\\ gamma^2 = T(alpha,beta), are mapped by the leading term of N_m along a curve through N with that tangent direction:
\\ (P_m(N))_m if nonzero (then r is contracted), else (dP_m(N).(alpha,beta) + gamma Q_m(N))_m
rpts = List(); for(al=-12,12, for(be=-12,12, if(gcd(al,be) != 1, next); my(g2 = Tcone(al, be), g); if(g2 > 0 && issquare(g2, &g), listput(rpts, [al,be,g]); listput(rpts, [al,be,-g]))));
rpts = Vec(rpts);
v0 = vector(3, m, substvec(f[m][1], [x0,x1,x2], NODE));
rimg(p) = if(v0 != [0,0,0], v0, vector(3, m, substvec(deriv(f[m][1], x0), [x0,x1,x2], NODE)*p[1] + substvec(deriv(f[m][1], x1), [x0,x1,x2], NODE)*p[2] + p[3]*substvec(f[m][2], [x0,x1,x2], NODE)));
\\ with two conics through the node (D'.r = 2) the sections vanish to order 2 at N: the first-order term is 0 and r maps to a
\\ conic (not computed; nodeimg = 0); the images are only needed for D'.r <= 1
\\ zero first-order images (a point of r on a contracted curve) are skipped, as in imgcurve; the
\\ line fitted from the nonzero ones is verified on all of them and again downstream (S_new a square on it)
rimgs = if(degr <= 1, my(Q = List()); for(n=1,#rpts, my(q = rimg(rpts[n])); if(q != [0,0,0], listput(Q, prim(q)))); Set(Q), []);
{
if(degr == 0, if(#rimgs != 1, error("r has D'.r = 0 but its image is not a point")); nodeimg = rimgs[1]; listput(newpoints, [["r"], nodeimg]),
  if(degr == 1, if(#rimgs < 3, error("too few image points of r")); nodeimg = prim(cross(rimgs[1], rimgs[2]));
    for(n=3,#rimgs, if(nodeimg*rimgs[n]~ != 0, error("r has D'.r = 1 but its image is not a line"))); listput(newlines, [["r"], nodeimg]),
    nodeimg = 0));      \\ D'.r = 2: r maps to a conic of the new model
}
newlines = Vec(newlines); newpoints = Vec(newpoints);
print(#newlines, " known curves of degree 1 -> lines; ", #newpoints, " known curves of degree 0 -> points: ", newpoints, "; D'.r = ", degr, if(degr == 0, concat(", r contracted to ", Str(nodeimg)), if(degr == 1, concat(", r -> the line ", Str(nodeimg)), ", r -> a conic")));
\\ cross-check when r is contracted: the images of the known curves through r (node-line components, strict transforms of
\\ conics through the node) of degree 1 are lines through the image of r
{
if(degr == 0, my(thru_r = select(nl -> (nl[1][1] == "L" && LIN[nl[1][2]][1][3] == 0) || (nl[1][1] == "C" && isnodeconic(nl[1][2])), newlines));
  for(n=1,#thru_r, if(thru_r[n][2]*nodeimg~ != 0, error("the image of a curve through r does not pass through the image of r")));
  print("image of r: ", nodeimg, " (", #Set(vector(#thru_r, n, thru_r[n][2])), " distinct known lines through it)"));
}
\\ checks: S_new is a square on every new line, singular at every new point
nsq = sum(n=1,#newlines, issqline(S, newlines[n][2]));
print("S_new is a square on ", nsq, " of the ", #newlines, " lines", if(nsq == #newlines, "  OK", "  FAILURE"));
if(nsq != #newlines, error("S_new is not a square on ", #newlines - nsq, " of the images of the known curves of degree 1: the model is wrong (a failed witness check aborts, nothing is written)"));
{
my(ok = 1); for(n=1,#newpoints, if(gradat(S, newpoints[n][2]) != [0,0,0], ok = 0; print("NOT singular at ", newpoints[n])));
print("S_new singular at the images of the contracted curves: ", if(ok, "OK", "FAILURE"));
if(!ok, error("S_new is not singular at the image of a contracted curve: the model is wrong (nothing is written)"));
}
\\ distinct lines, and which pass through a singular point (node lines) or not (tritangents)
sing = Set(vector(#newpoints, n, newpoints[n][2]));
LNEW = Set(vector(#newlines, n, newlines[n][2]));
nodel = select(l -> vecmin(vector(#sing, m, abs(l*sing[m]~))) == 0, LNEW); tritl = setminus(LNEW, Set(nodel));
print(#LNEW, " distinct new lines: ", #tritl, " tritangents (through no known singular point) and ", #nodel, " lines through a singular point");
writefile(concat(["out/model_", NAME, "_raw.gp"]), Str("\\\\ raw model of D' = C_i + sign C_j (models/model.gp) for PAIR = ", PAIR, "\nNAME = \"", NAME, "\";\nPAIR = ", PAIR, ";\nf = ", f, ";\nc = ", c, ";\nS = ", S, ";\nnewlines = ", newlines, ";\nnewpoints = ", newpoints, ";\nnodeimg = ", nodeimg, ";\nsing = ", sing, ";\ntrit_raw = ", tritl, ";\nnodelines_raw = ", nodel, ";"));
writefile(concat(["out/model_", NAME, "_raw_coeffs.txt"]), Str(c));
print("written out/model_", NAME, "_raw.gp and out/model_", NAME, "_raw_coeffs.txt  (total ", (getabstime()-T0)/1000., " s)");
\q

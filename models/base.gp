/* base.gp : the lattice side of a one-node model B of X_163, from its sextic, its square lines and its conics alone.
   1. Neron-Severi lattice from the lifts of B's square lines and the exceptional curve r (rank 20, |det| = 163 asserted);
   2. the class of each conic of B's table (every intersection number with the line lifts is reproduced);
   3. all classes v with v^2 = -2 and v.D_B = 0, 1, 2 (enumeration); checks: no class E with E^2 = 0, E.D_B = 1, the classes of
      degree 1 are the 2 NT + 4 NN line classes, the irreducible classes of degree 2 are exactly the 2 NC conic lifts
      (so the table of lines and the table of conics are complete);
   4. writes OUTD/lat_<B>.txt for pairs.py: GZ, D_B, r, the 2 NC conic lifts ("X i sign thru height z1..z20") and the
      classes of 3. ("R degree z1..z20").
   usage (from the root): gp -q with BNAME = "E" (Elkies: data/elkies163*.txt) or "<M>" (data/x163_<M>*.txt), OUTD = "<dir>";
   the node of B must be (0:0:1). */
default(parisizemax, 8000000000);
if(type(BNAME) != "t_STR", error("set BNAME"));
if(type(OUTD) != "t_STR", error("set OUTD"));
BF = if(BNAME == "E", "data/elkies163", concat("data/x163_", BNAME));
tokens(str) = { my(t = strsplit(str, " "), v = []); for(i=1,#t, if(#t[i], v = concat(v, [eval(t[i])]))); v }
readtab(f) = { my(l = readstr(f), R = List()); for(i=1,#l, if(#l[i] && Vecsmall(l[i])[1] != 35, listput(R, tokens(l[i])))); Vec(R) }
\\ 1. ---- the sextic and its lines ----
cB = readtab(concat(BF, ".txt"))[1];
if(#cB != 28, error("sextic record: ", #cB, " fields"));
MON = []; for(i=0,6, my(ii=6-i); for(jj=0,6-ii, my(j=6-ii-jj, k=6-ii-j); MON = concat(MON, [[ii,j,k]])));
SB = sum(n=1,28, cB[n]*x0^MON[n][1]*x1^MON[n][2]*x2^MON[n][3]);
Sx(v) = substvec(SB, [x0,x1,x2], v);
node = [0,0,1];
if(substvec(SB, [x0,x1,x2], node) != 0, error("node not on the sextic"));
{ my(g = [substvec(deriv(SB,x0),[x0,x1,x2],node), substvec(deriv(SB,x1),[x0,x1,x2],node), substvec(deriv(SB,x2),[x0,x1,x2],node)]);
  if(g != [0,0,0], error("(0:0:1) is not singular on S_", BNAME)); }
LNS = readtab(concat(BF, "_lines.txt"));
{ my(tr = List(), nd = List()); for(i=1,#LNS, if(#LNS[i] != 3, error("line record"));
    if(LNS[i][3] == 0, listput(nd, LNS[i]), listput(tr, LNS[i])));
  trit = Vec(tr); nodel = Vec(nd); }
NT = #trit; allL = concat(trit, nodel); NL = #allL;
print("base ", BNAME, ": ", NT, " tritangents, ", NL - NT, " node lines");
param(l) = { my(p, q); if(l[3]!=0, p = [1,0,-l[1]/l[3]]; q = [0,1,-l[2]/l[3]], if(l[2]!=0, p=[1,-l[1]/l[2],0]; q=[0,0,1], p=[0,1,0]; q=[0,0,1]));
  my(d=1); d = lcm(d, denominator(p)); p = p*d; d=lcm(1,denominator(q)); q=q*d; [p,q] }
cubic(pq) = { my(f = Sx(pq[1]*v0 + pq[2]*v1), c); if(!issquare(f, &c), error("not a square")); c }
PQ = vector(NL, i, param(allL[i]));
CU = vector(NL, i, cubic(PQ[i]));
for(i=1,NL, if(CU[i]^2 != Sx(PQ[i][1]*v0+PQ[i][2]*v1), error("cubic check")));
cval(i, x) = { my(pq=PQ[i], M = matconcat([pq[1]~, pq[2]~]), sol, ok=0);
  for(a=1,3, for(b=a+1,3, my(Ms = matrix(2,2,u,t,M[[a,b][u],t])); if(matrank(Ms)==2, sol = matsolve(Ms, [x[a],x[b]]~); ok=1; break(2))));
  if(!ok, error("param")); if(M*sol != x~, error("point not on line")); substvec(CU[i],[v0,v1],sol~) }
cross(l,m) = [l[2]*m[3]-l[3]*m[2], l[3]*m[1]-l[1]*m[3], l[1]*m[2]-l[2]*m[1]];
for(i=NT+1, NL, if(cval(i, node) != 0, error("node line cubic does not vanish at node")));
\\ Gram of the 2 + 2 NL classes (1 = D_B, 2 = r, 2+2i-1 = l_i+, 2+2i = l_i-)
nc = 2 + 2*NL;
G = matrix(nc, nc); G[1,1]=2; G[2,2]=-2;
{ for(i=1,NL, my(ip = 2+2*i-1, im = 2+2*i, isnode = (i > NT));
    G[1,ip]=1; G[ip,1]=1; G[1,im]=1; G[im,1]=1;
    G[2,ip]=isnode; G[ip,2]=isnode; G[2,im]=isnode; G[im,2]=isnode;
    G[ip,ip]=-2; G[im,im]=-2; G[ip,im]=if(isnode,2,3); G[im,ip]=G[ip,im]);
  for(i=1,NL, for(k=i+1,NL, my(ip=2+2*i-1, im=2+2*i, kp=2+2*k-1, km=2+2*k);
    my(P = cross(allL[i], allL[k])); P = P/content(P);
    if(P == node || P == -node,
       G[ip,kp]=0;G[ip,km]=0;G[im,kp]=0;G[im,km]=0; G[kp,ip]=0;G[km,ip]=0;G[kp,im]=0;G[km,im]=0; next);
    if(Sx(P)==0, error("intersection point on the branch curve: ", allL[i], allL[k]));
    my(ci = cval(i,P), ck = cval(k,P));
    if(ci^2 != ck^2, error("values"));
    if(ci == ck, G[ip,kp]=1; G[im,km]=1; G[ip,km]=0; G[im,kp]=0, G[ip,kp]=0; G[im,km]=0; G[ip,km]=1; G[im,kp]=1);
    G[kp,ip]=G[ip,kp]; G[km,im]=G[im,km]; G[km,ip]=G[ip,km]; G[kp,im]=G[im,kp])); }
print("rank of the Gram matrix of the ", nc, " classes: ", matrank(G));
if(matrank(G) != 20, error("span has rank ", matrank(G), ", not 20"));
{ for(i=1,NL, my(ip=2+2*i-1, im=2+2*i); if(i<=NT, if(G[ip,]+G[im,] != G[1,], error("relation l+ + l- = D fails ", i)), if(G[ip,]+G[im,]+G[2,] != G[1,], error("relation E+E'+r = D fails ", i)))); }
\\ Z-basis
{ my(ix = List()); for(j=1,nc, my(t = concat(Vec(ix),[j])); if(matrank(vecextract(G,t,t)) == #t, listput(ix, j)); if(#ix==20, break));
  idx = Vec(ix); }
GB = vecextract(G,idx,idx);
{ my(C = GB^-1 * vecextract(G, idx, vector(nc,j,j)), d = 1);
  for(j=1,nc, d = lcm(d, denominator(C[,j])));
  Hb = mathnf(C*d)/d; }
GZ = Hb~*GB*Hb;
print("det GZ = ", matdet(GZ));
if(abs(matdet(GZ)) != 163, error("the lines of ", BNAME, " span det ", matdet(GZ), ", not +-163"));
Hi = Hb^-1;
Zcls = matrix(20, nc);
{ my(Cc = GB^-1 * vecextract(G, idx, vector(nc, j, j))); Zcls = Hi * Cc;
  if(denominator(Zcls) != 1, error("class coordinates not integral")); }
DB = Zcls[,1]; EB = Zcls[,2];
if(DB~*GZ*DB != 2 || EB~*GZ*EB != -2 || DB~*GZ*EB != 0, error("D_B/E_B sanity"));
\\ 2. ---- conic classes ----
CONraw = readtab(concat(BF, "_conics.txt"));
{ CON = List(); for(i=1,#CONraw, my(v = CONraw[i]); if(#v < 15, error("conic record ", i, ": ", #v, " fields"));
    if(v[7] == 0 && v[10] == 0, for(k=0,2, my(a = v[7+3*k], c = v[9+3*k]); v[7+3*k] = c; v[9+3*k] = a));
    listput(CON, v)); CON = Vec(CON); }
NC = #CON; print(NC, " conics");
xC(i) = { my(m = CON[i]); vector(3, k, m[6+3*(k-1)+1]*s^2 + m[6+3*(k-1)+2]*s + m[6+3*(k-1)+3]) }
qC(i) = { my(c = CON[i]); c[1]*x0^2 + c[2]*x0*x1 + c[3]*x0*x2 + c[4]*x1^2 + c[5]*x1*x2 + c[6]*x2^2 }
gC = vector(NC); thrunode = vector(NC);
{ for(i=1,NC, my(x = xC(i), f = Sx(x), g); if(!issquare(f, &g), error("S not a square on conic ", i)); gC[i] = g;
    thrunode[i] = (CON[i][6] == 0);
    if(substvec(qC(i), [x0,x1,x2], x) != 0, error("parametrisation not on the conic ", i))); }
print("lifts of all conics found; conics through the node: ", vecsum(thrunode));
\\ line cubic FORMS Gl[k] (cubic in x restricting to CU[k] on line k) for the conic-line intersections
CMON = [x0^3, x0^2*x1, x0^2*x2, x0*x1^2, x0*x1*x2, x0*x2^2, x1^3, x1^2*x2, x1*x2^2, x2^3];
Gl = vector(NL);
{ for(i=1,NL, my(pq = PQ[i], M = matrix(4, 10), rhs = vector(4)~);
    for(m=1,10, my(r = substvec(CMON[m], [x0,x1,x2], pq[1]*v0 + pq[2]*v1)); for(d=0,3, M[d+1,m] = polcoef(polcoef(r, d, v0), 3-d, v1)));
    for(d=0,3, rhs[d+1] = polcoef(polcoef(CU[i], d, v0), 3-d, v1));
    my(sol = matinverseimage(M, rhs)); if(type(sol) != "t_COL", error("no cubic form for line ", i));
    Gl[i] = sum(m=1,10, sol[m]*CMON[m])); }
\\ node machinery (node at (0:0:1)): r-point of a curve through the node
nodept(x, g, s0) = { my(a = deriv(x[1], s), b = deriv(x[2], s), c = x[3], gp = deriv(g, s), P);
  P = [subst(a, s, s0)*subst(c, s, s0)^2, subst(b, s, s0)*subst(c, s, s0)^2, subst(gp, s, s0)];
  if(P[1] == 0 && P[2] == 0, error("degenerate node point")); P / content(P) * sign(if(P[1] != 0, P[1], P[2])) }
nodeinfo = vector(NC);
{ for(i=1,NC, if(!thrunode[i], next); my(x = xC(i), r0 = gcd(x[1], x[2]));
    if(poldegree(r0) != 1, if(poldegree(r0) == 0 && polcoef(x[1],2,s) == 0 && polcoef(x[2],2,s) == 0, error("node at infinity for conic ", i, " (not handled)"), error("conic ", i, " meets the node ", poldegree(r0), " times")));
    my(s0 = -polcoef(r0,0,s)/polcoef(r0,1,s)); if(subst(x[3], s, s0) == 0, error("bad node parameter"));
    nodeinfo[i] = [s0, nodept(x, gC[i], s0)]); }
nodeline_pt = vector(NL);
{ for(k=NT+1,NL, my(pq = PQ[k], M = matconcat([pq[1]~, pq[2]~]), sol, ok = 0);
    for(a=1,3, for(b=a+1,3, my(Ms = matrix(2,2,u,t,M[[a,b][u],t])); if(matrank(Ms)==2, sol = matsolve(Ms, [[0,0,1][a], [0,0,1][b]]~); ok=1; break(2))));
    if(!ok || M*sol != [0,0,1]~, error("node not on node line ", k));
    my(d = if(sol[1] != 0, [0,1], [1,0]), x = vector(3, m, pq[1][m]*(sol[1] + s*d[1]) + pq[2][m]*(sol[2] + s*d[2])), g = substvec(CU[k], [v0,v1], [sol[1] + s*d[1], sol[2] + s*d[2]]));
    if(subst(x[1], s, 0) != 0 || subst(x[2], s, 0) != 0, error("node parameter"));
    nodeline_pt[k] = nodept(x, g, 0)); }
rv(f, d) = s^(d - poldegree(f)) * polrecip(f);
ramflag = 0;
inter_line(i, k) = { my(x = xC(i), h = allL[k][1]*x[1] + allL[k][2]*x[2] + allL[k][3]*x[3], gl = substvec(Gl[k], [x0,x1,x2], x), g = gC[i]);
  my(np = 0, nm = 0, dh, h0 = h, dh0 = poldegree(h), nodeadd = [0, 0]); if(h == 0, error("conic on a line"));
  if(thrunode[i] && k > NT,
     my(s0 = nodeinfo[i][1], m = 0); while(subst(h, s, s0) == 0, h = h \ (s - s0); m++);
     if(m >= 2, my(P = nodeinfo[i][2], Q = nodeline_pt[k], Qm = [Q[1], Q[2], -Q[3]] / content([Q[1], Q[2], -Q[3]]) * sign(if(Q[1] != 0, Q[1], Q[2])));
        if(P == Q, nodeadd[1] = m - 1, if(P == Qm, nodeadd[2] = m - 1))));
  dh = poldegree(h);
  if(poldegree(gcd(h, g)) > 0, ramflag++);
  np = poldegree(gcd(h, g - gl)); nm = poldegree(gcd(h, g + gl));
  if(dh0 < 2,
     my(hr = rv(h0, 2), gr = rv(g, 6), glr = rv(gl, 6), vh = valuation(hr, s));
     if(valuation(gr, s) > 0 && valuation(glr, s) > 0, ramflag++);
     np += min(vh, valuation(gr - glr, s)); nm += min(vh, valuation(gr + glr, s)));
  np += nodeadd[1]; nm += nodeadd[2];
  my(tot = if(thrunode[i] && k > NT, 1, 2)); if(np + nm != tot, error("C.l+ + C.l- != ", tot, " for conic ", i, " line ", k, ": ", np, " ", nm)); [np, nm] }
Zcon = matrix(20, NC);
{ for(i=1,NC,
    my(v = vector(nc)~); v[1] = 2; v[2] = thrunode[i];
    for(k=1,NL, my(r = inter_line(i, k)); v[2+2*k-1] = r[1]; v[2+2*k] = r[2]);
    my(c = GB^-1 * vecextract(v, idx), z = Hi * c);
    if(denominator(z) != 1, error("non-integral class for conic ", i));
    if(z~*GZ*z != -2, error("C^2 != -2 for conic ", i, ": ", z~*GZ*z));
    my(chk = (z~*GZ*Zcls)~); if(chk != v, error("intersection numbers not reproduced for conic ", i));
    Zcon[,i] = z); }
print("classes of all ", NC, " conics computed and verified; ramification flags: ", ramflag);
{ my(fh = fileopen(concat([OUTD, "/conic_classes_", BNAME, ".txt"]), "w"));
  filewrite(fh, concat(["# i thrunode z1..z20 (class of the lift w = +sqrt of conic i of data/x163_", BNAME, "_conics.txt in the Z-basis of this run of baseB_classes.gp)"]));
  for(i=1,NC, filewrite(fh, Str(i, " ", thrunode[i], " ", concat(vector(20, k, Str(Zcon[k,i], " ")))))); fileclose(fh); }
\\ 3. ---- all (-2)-classes of degree 0, 1, 2 ----
T0 = getabstime();
gD = GZ*DB;
Phi = -GZ + gD*gD~;                                   \\ Phi(v) = -v^2 + (v.D)^2, positive definite
{ my(q = qfminim(Phi, 6), V = q[3], n = #V, R0 = List(), R1 = List(), R2 = List(), nell = 0);
  for(k = 1, n, my(v = V[,k], t = gD~*v, s = v~*GZ*v);
    if(t < 0, v = -v; t = -t);
    if(s == 0 && t == 1, nell++);
    if(s != -2 || t > 2, next);
    if(t == 0, listput(R0, v), if(t == 1, listput(R1, v), listput(R2, v))));
  if(nell, error(nell, " classes E with E^2 = 0, E.D = 1: not a polarization without base points"));
  ROOT0 = Vec(R0); ROOT1 = Vec(R1); ROOT2 = Vec(R2);
  print("classes with v^2 = -2: degree 0: ", #ROOT0, " (up to sign), degree 1: ", #ROOT1, ", degree 2: ", #ROOT2, "  (", n, " vectors enumerated, ", getabstime() - T0, " ms)"); }
if(#ROOT0 != 1 || (ROOT0[1] != EB && ROOT0[1] != -EB), error("the classes of degree 0 are not +-r"));
if(#ROOT1 != 2*NT + 4*(NL - NT), error("classes of degree 1: ", #ROOT1, ", lines give ", 2*NT + 4*(NL - NT)));
\\ irreducible classes of degree 2: non-negative on r and on every class of degree 1
{ my(M1 = matconcat(ROOT1), nirr = 0, isconic = Map());
  for(i = 1, NC, mapput(isconic, Str(Zcon[,i]~), 1); mapput(isconic, Str((2*DB - thrunode[i]*EB - Zcon[,i])~), 1));
  for(k = 1, #ROOT2, my(v = ROOT2[k], w = (v~*GZ*M1), ok = (v~*GZ*EB >= 0));
    if(ok, for(j = 1, #w, if(w[j] < 0, ok = 0; break)));
    if(ok, nirr++; if(!mapisdefined(isconic, Str(v~)), error("an irreducible class of degree 2 is not a conic of the table: ", v~))));
  print("irreducible classes of degree 2: ", nirr, "; conic lifts: ", 2*NC);
  if(nirr != 2*NC, error("the conic table is not complete")); }
\\ 4. ---- output ----
hq(i) = vecmax(abs(CON[i][1..6]));
{ my(fh = fileopen(concat([OUTD, "/lat_", BNAME, ".txt"]), "w"), row(v) = concat(vector(#v, k, Str(v[k], " "))));
  filewrite(fh, concat("GZ ", row(concat(Vec(GZ~)))));
  filewrite(fh, concat("D ", row(DB~)));
  filewrite(fh, concat("E ", row(EB~)));
  for(i = 1, NC,
    filewrite(fh, concat(["X ", Str(i), " 1 ", Str(thrunode[i]), " ", Str(hq(i)), " ", row(Zcon[,i]~)]));
    my(cj = 2*DB - thrunode[i]*EB - Zcon[,i]);
    if(cj~*GZ*cj != -2 || cj~*GZ*DB != 2 || cj~*GZ*EB != thrunode[i], error("conj ", i, " sanity"));
    filewrite(fh, concat(["X ", Str(i), " -1 ", Str(thrunode[i]), " ", Str(hq(i)), " ", row(cj~)])));
  filewrite(fh, concat("R 0 ", row(EB~)));
  for(k = 1, #ROOT1, filewrite(fh, concat("R 1 ", row(ROOT1[k]~))));
  for(k = 1, #ROOT2, filewrite(fh, concat("R 2 ", row(ROOT2[k]~))));
  fileclose(fh); }
print("BASEDONE ", BNAME);
\q

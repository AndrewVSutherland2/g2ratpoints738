/* finalize.gp : verification and output of a reduced model out/model_<NAME>_red.gp:
   - singular points from Magma (out/model_<NAME>_sing.gp, written by sing.m; the images of the contracted known curves
     must be among them), the exhaustive list of square lines of height <= H (out/model_<NAME>_sq<H>.txt from sqlines.sh),
     classified into tritangents (through no singular point) and node lines (through a singular point);
   - the known lines (images of the known curves of degree 1) of height <= H must all be found by the search;
   - on five random lines: all forced x-values (intersections with the tritangents and node lines) are squares, and
     the rational points of height <= 10^4 of the genus-2 curve are counted (Elkies' model on the same line for comparison);
   - a_p against the weight-3 CM form of level 163 for 47 <= p <= PMAX.
   Writes data/x163_<NAME>.txt (28 coefficients) and data/x163_<NAME>_lines.txt (tritangents, then node lines), and
   out/model_<NAME>.gp (c, S, sing, trit, nodelines, M, PAIR).
   usage: ./run.sh finalize.gp logs/finalize_<NAME>.log 'NAME="<NAME>"' [H=150] [PMAX=100] [WRITE=1] */
\r cplib.gp
if(type(NAME) != "t_STR", error("set NAME"));
H = if(type(H) == "t_INT", H, 150); PMAX = if(type(PMAX) == "t_INT", PMAX, 100); WRITE = if(type(WRITE) == "t_INT", WRITE, 1);
read(concat(["out/model_", NAME, "_red.gp"]));
S = sextic(c);
print("model ", NAME, " (pair ", PAIR, "): max |coeff| ", vecmax(abs(c)), ", content ", content(c));
\\ singular points (Magma)
singfile = concat(["out/model_", NAME, "_sing.gp"]);
if(type(fileopen(singfile, "r")) != "t_INT", error("run  magma -b file:=out/model_", NAME, "_red_coeffs.txt out:=", singfile, " sing.m  first"));
read(singfile);
singtypes = vector(#singpts, n, if(singmult[n] == 2, concat("A", Str(singtau[n])), concat(["mult", Str(singmult[n]), "tau", Str(singtau[n])])));
print("Magma: singular subscheme of degree ", singdeg, ", rational singular points ", singpts, " of types ", singtypes, if(#singnonrat > 0, concat(", NON-RATIONAL singular components [degree, length] ", Str(singnonrat)), ""));
for(n=1,#sing, if(!vecsearch(vecsort(singpts), sing[n]), error("the image ", sing[n], " of a contracted known curve is not among Magma's singular points: the witnesses and the model disagree (fatal, no data files are written)")));
for(n=1,#singpts, if(gradat(S, singpts[n]) != [0,0,0], error("Magma's singular point ", singpts[n], " is not singular?")));
if(#singpts == 0, error("no rational singular point"));
ALLSING = singpts;
\\ square lines of height <= H (exhaustive search)
found = readlines(concat(["out/model_", NAME, "_sq", Str(H), ".txt"]));
found = select(l -> issqline(S, l), found);      \\ re-verified in gp (the zero form is never a square here)
thru(l) = #ALLSING > 0 && vecmin(vector(#ALLSING, m, abs(l*ALLSING[m]~))) == 0;
ft = select(l -> !thru(l), found); fn = select(l -> thru(l), found);
print("exhaustive search up to height ", H, ": ", #found, " square lines = ", #ft, " tritangents + ", #fn, " node lines");
\\ the known lines: all of height <= H must have been found; those of larger height are added
known = Set(concat(trit, nodelines));
missed = select(l -> vecmax(abs(l)) <= H && !vecsearch(vecsort(found), l), known);
if(#missed > 0, error("known lines of height <= ", H, " missed by the search: ", missed));
extra = select(l -> vecmax(abs(l)) > H, known);
print(#known, " known lines (images of known curves): all of height <= ", H, " found; ", #extra, " of larger height added: ", extra);
alll = Set(concat(found, extra)); tl = select(l -> !thru(l), alll); nl = select(l -> thru(l), alll);
print("TOTAL: ", #tl, " tritangents + ", #nl, " node lines = ", #alll, " forced lines -> ", 2*#alll, " forced points on a generic line; max |coeff| of the lines ", vecmax(vector(#alll, n, vecmax(abs(alll[n])))));
print("node lines per singular point: ", vector(#ALLSING, m, #select(l -> l*ALLSING[m]~ == 0, nl)));
\\ random lines: forced x-values are squares; point counts to 10^4
SEl = SE; BASELBL = if(BASE == "", "Elkies' model", concat(["base ", BASE]));
countpts(Sp, l, Hb) = { my(F = restrict(Sp, l), fu = subst(subst(F, v, 1), u, x)); if(!issquarefree(fu), return(-1)); my(pts = hyperellratpoints(fu, Hb)); sum(n=1,#pts, if(pts[n][2] == 0, 1, 2)) }
setrand(2026);
{
for(k=1,5, my(l = vector(3, m, random(21) - 10)); if(l == [0,0,0] || l[3] == 0 || thru(l), k--; next);
  my(F = restrict(S, l), fu = subst(subst(F, v, 1), u, x)); if(!issquarefree(fu), print("line ", l, ": singular"); next);
  my(nsq = 0);
  for(n=1,#alll, my(m = alll[n], A = m[1]*l[3] - m[3]*l[1], Bc = m[2]*l[3] - m[3]*l[2], xv);
    xv = if(A == 0, oo, -Bc/A);
    nsq += if(xv == oo, issquare(pollead(fu)), issquare(subst(fu, x, xv))));
  print("line ", l, ": forced x-values that are squares: ", nsq, " of ", #alll, " (=> ", 2*nsq, " forced points); points of height <= 10^4: ", countpts(S, l, 10^4), " (", BASELBL, " on the same line: ", countpts(SEl, l, 10^4), ")"));
}
\\ a_p check: #X(F_p) + (number of exceptional curves) p = 1 + 20 p + p^2 + a_p (each rational A_n point adds n rational
\\ exceptional curves for odd n; for even n the middle curves are defined over F_p too since the point is rational and
\\ the tangent cone... -- only used with rational A1 points here: k nodes add k p)
nexc = sum(n=1,#singtau, singtau[n]);
if(#singnonrat > 0, print("WARNING: non-rational singular components: the a_p correction below assumes only the rational singular points"));
npts(Sp, p) = { my(cnt = 0, Spp = Sp*Mod(1,p), vv);
  for(a=0,p-1, for(b=0,p-1, vv = lift(substvec(Spp,[x0,x1,x2],[a,b,1])); cnt += 1 + kronecker(vv,p)));
  for(a=0,p-1, vv = lift(substvec(Spp,[x0,x1,x2],[a,1,0])); cnt += 1 + kronecker(vv,p));
  vv = lift(substvec(Spp,[x0,x1,x2],[1,0,0])); cnt += 1 + kronecker(vv,p);
  cnt }
cmap(p) = { if(kronecker(-163, p) != 1, return(0)); my(a); for(bb=1, sqrtint(4*p), my(tt = 4*p - 163*bb^2); if(tt >= 0 && issquare(tt, &a), return((a^2 - 163*bb^2)/2))); error("no rep") }
\\ number of singular points of S mod p in P^2(F_p) (a prime with more of them than the rational singular points is a prime
\\ of bad reduction of the model, where the point count is off by a multiple of p: not a failure)
nsingmodp(Sp, p) = { my(cnt = 0, Spp = Sp*Mod(1,p), g = vector(3, m, deriv(Spp, [x0,x1,x2][m])), pts = List());
  forvec(v = [[0,p-1],[0,p-1]], listput(pts, [v[1], v[2], 1])); for(a=0,p-1, listput(pts, [a, 1, 0])); listput(pts, [1, 0, 0]);
  for(n=1,#pts, my(P = pts[n]); if(substvec(Spp, [x0,x1,x2], P) == 0 && vector(3, m, substvec(g[m], [x0,x1,x2], P)) == [0,0,0], cnt++)); cnt }
{
my(bad = 0, badp = List()); forprime(p = 47, PMAX, my(n = npts(S, p) + nexc*p, ap = n - 1 - p^2 - 20*p, cm = cmap(p)); if(ap != cm,
  my(ns = nsingmodp(S, p)); if((ap - cm) % p == 0 && ns > #ALLSING, listput(badp, [p, (ap - cm)/p, ns]), bad++; print("p = ", p, ": a_p = ", ap, " but the CM form gives ", cm, " (", ns, " singular points mod p)"))));
print("a_p check against the CM form for 47 <= p <= ", PMAX, " (", nexc, " exceptional curves): ", if(bad == 0, "OK", concat(Str(bad), " mismatches")), if(#badp > 0, concat(["; primes of bad reduction (count off by k p, extra singular points mod p) [p, k, #sing mod p]: ", Str(Vec(badp))]), ""));
if(bad > 0, error("a_p verification failed at ", bad, " good prime(s): the model contradicts the CM form (fatal, nothing is written)"));
}
{
if(WRITE,
  my(F = fileopen(concat([REPO, "/data/x163_", NAME, ".txt"]), "w")); filewrite(F, concat(vector(28, n, concat(Str(c[n]), " ")))); fileclose(F);
  F = fileopen(concat([REPO, "/data/x163_", NAME, "_lines.txt"]), "w");
  filewrite(F, concat(["# lines on which the X_163 model data/x163_", NAME, ".txt is a perfect square: ", Str(#tl), " tritangents + ", Str(#nl), " node lines (singular points ", Str(ALLSING), "; exhaustive search to height ", Str(H), " plus the images of the known curves)"]));
  for(n=1,#tl, filewrite(F, concat([Str(tl[n][1]), " ", Str(tl[n][2]), " ", Str(tl[n][3])])));
  for(n=1,#nl, filewrite(F, concat([Str(nl[n][1]), " ", Str(nl[n][2]), " ", Str(nl[n][3])])));
  fileclose(F);
  writefile(concat(["out/model_", NAME, ".gp"]), Str("\\\\ double-plane model y^2 = S of X_163 for D' = C_i + sign C_j, PAIR = ", PAIR, " (models): ", #tl, " tritangents + ", #nl, " node lines, singular points ", ALLSING, " of types ", singtypes, "\nNAME = \"", NAME, "\";\nPAIR = ", PAIR, ";\nc = ", c, ";\nS = ", S, ";\nsing = ", ALLSING, ";\nsingtypes = ", singtypes, ";\ntrit = ", tl, ";\nnodelines = ", nl, ";\nM = ", M, ";"));
  print("written data/x163_", NAME, ".txt, data/x163_", NAME, "_lines.txt, out/model_", NAME, ".gp"),
  writefile(concat(["out/model_", NAME, ".gp"]), Str("\\\\ double-plane model y^2 = S of X_163 for D' = C_i + sign C_j, PAIR = ", PAIR, " (models; validation model, no data file): ", #tl, " tritangents + ", #nl, " node lines, singular points ", ALLSING, " of types ", singtypes, "\nNAME = \"", NAME, "\";\nPAIR = ", PAIR, ";\nc = ", c, ";\nS = ", S, ";\nsing = ", ALLSING, ";\nsingtypes = ", singtypes, ";\ntrit = ", tl, ";\nnodelines = ", nl, ";\nM = ", M, ";"));
  print("written out/model_", NAME, ".gp (WRITE = 0: no data files)"));
}
\q

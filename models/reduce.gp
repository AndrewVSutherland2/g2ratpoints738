/* reduce.gp : coefficient reduction of the raw model out/model_<NAME>_raw.gp by projective coordinate changes
   (reducelib.gp): (1) a frame made of three of the known new lines (all triples are
   tried, the best NKEEP by the size of the transformed coefficients are kept), (2) torus rescaling by small primes,
   (3) local search over shears and diagonal scalings (localsearch), and the transformed lines / singular points.
   Output: out/model_<NAME>_red.gp (c, S, M, sing, trit, nodelines, nodeimg; y = M x, S_red(y) = S_raw(M^-1 y) up to a
   square factor).  usage: ./run.sh reduce.gp logs/reduce_<NAME>.log 'NAME="<NAME>"' [NKEEP=6] */
\r cplib.gp
\r reducelib.gp
if(type(NAME) != "t_STR", error("set NAME"));
NKEEP = if(type(NKEEP) == "t_INT", NKEEP, 6);
read(concat(["out/model_", NAME, "_raw.gp"]));
Sraw = S; craw = c;
LNEW = Set(concat(trit_raw, nodelines_raw));
print(#LNEW, " known lines (", #trit_raw, " tritangents, ", #nodelines_raw, " through a singular point); raw model: ", sizev(craw), " (max digits, sum of digits)");
T0 = getabstime();
\\ 1. frames from triples of known lines
cands = List();
{
my(n = #LNEW, cnt = 0);
if(type(FRAMES) == "t_STR",
  \\ the program frames (frames.c) scores all triples and returns the best ones, which are recomputed here exactly
  my(fin = concat(["out/model_", NAME, "_frames.in"]), F = fileopen(fin, "w"));
  filewrite(F, concat(vector(28, k, Str(craw[k], " ")))); filewrite(F, Str(n));
  for(a=1,n, filewrite(F, Str(LNEW[a][1], " ", LNEW[a][2], " ", LNEW[a][3]))); fileclose(F);
  my(res = externstr(concat([FRAMES, " ", Str(NKEEP), " < ", fin])), seen = 0);
  for(k=1,#res, my(w = strsplit(res[k], " "));
    if(w[1] == "N", cnt = eval(w[2]); seen = 1; next);
    if(w[1] != "C" && w[1] != "F", error("frames: unexpected output ", res[k]));
    my(a = eval(w[2]), b = eval(w[3]), d = eval(w[4]), M = framemat(LNEW[a], LNEW[b], LNEW[d]));
    my(c1 = normvec(transformvec(Sraw, M)), s1 = scorev(c1));
    listput(cands, [s1, [a,b,d], M, c1]));
  if(!seen || #cands < min(NKEEP, cnt), error("frames: incomplete output")),
for(a=1,n, for(b=a+1,n, for(d=b+1,n, my(M = framemat(LNEW[a], LNEW[b], LNEW[d])); if(matdet(M) == 0, next); cnt++;
  my(c1 = normvec(transformvec(Sraw, M)), s1 = scorev(c1));
  listput(cands, [s1, [a,b,d], M, c1])))));
cands = vecsort(Vec(cands), (u, v) -> sign(u[1] - v[1]));
print(cnt, " frames tried (", (getabstime()-T0)/1000., " s); best scores: ", vector(min(NKEEP, #cands), k, cands[k][1]));
}
\\ 2. + 3. torus and shear reduction of the best frames
best = [10^9, 0, 0];
{
for(k=1,min(NKEEP, #cands), my(M0 = cands[k][3], c0 = cands[k][4], pr = relprimes(c0), T, c1, M1, LS, c2, M2);
  T = torussearch(c0, pr); c1 = T[1]; M1 = matdiagonal(vector(3, m, 1/T[2][m])) * M0;
  LS = localsearch(c1, matid(3), pr, 60, 0); c2 = LS[1]; M2 = LS[2]*M1;
  print("frame ", cands[k][2], ": ", sizev(c0), " -> torus ", sizev(c1), " -> shears ", sizev(c2), " max |coeff| ", vecmax(abs(c2)));
  if(scorev(c2) < best[1], best = [scorev(c2), M2, c2]));
}
M = best[2]; c = best[3];
c2chk = normvec(transformvec(Sraw, M)); if(c2chk != c, error("inconsistent transformation"));
S = polfromvec(c);
trit = translines(trit_raw, M); nodelines = translines(nodelines_raw, M);
sing = vector(#sing, n, prim(M*sing[n]~)~);
if(nodeimg != 0, nodeimg = prim(M*nodeimg~)~);     \\ nodeimg = 0: r maps to a conic (two conics through the node), nothing to transport
print("reduced model: max |coeff| ", vecmax(abs(c)), " (", #Str(vecmax(abs(c))), " digits), content ", content(c), "; singular points ", sing, " (image of r: ", nodeimg, ")");
print("squares: tritangents ", vector(#trit, n, issqline(S, trit[n])) == vector(#trit, n, 1), ", node lines ", vector(#nodelines, n, issqline(S, nodelines[n])) == vector(#nodelines, n, 1));
print("max |coefficient| of the known lines: ", vecmax(vector(#trit + #nodelines, n, vecmax(abs(concat(trit, nodelines)[n])))));
print("coefficients: ", c);
writefile(concat(["out/model_", NAME, "_red.gp"]), Str("\\\\ reduced model (models/reduce.gp) of ", NAME, ": y = M x, S(y) = S_raw(M^-1 y) up to a square\nNAME = \"", NAME, "\";\nPAIR = ", PAIR, ";\nc = ", c, ";\nS = ", S, ";\nM = ", M, ";\nsing = ", sing, ";\nnodeimg = ", nodeimg, ";\ntrit = ", trit, ";\nnodelines = ", nodelines, ";"));
writefile(concat(["out/model_", NAME, "_red_coeffs.txt"]), Str(c));
print("written out/model_", NAME, "_red.gp  (", (getabstime()-T0)/1000., " s)");
\q

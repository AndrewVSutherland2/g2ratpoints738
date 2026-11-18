/* minred_file.gp : the minimised and reduced version of an existing model file (see minred_file.sh).  Set SXF (sextic file),
   LNF (its file of square lines, "A B C" rows), NODES (its singular points, a vector of vectors), MRF (output of minred.m),
   OUTS, OUTL (output files) and NAME before reading this script. */
\r cplib.gp
\r reducelib.gp
read(MRF);
{
my(R = readstr(SXF), c0 = List(), L = List(), c1, S0, S1, k, nl, ns);
for(i = 1, #R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = strsplit(R[i], " ")); for(j = 1, #w, if(#w[j], listput(c0, eval(w[j])))));
c0 = Vec(c0); if(#c0 != 28, error("minred_file.gp: expected 28 coefficients"));
R = readstr(LNF); for(i = 1, #R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = select(t -> #t > 0, strsplit(R[i], " "))); listput(L, vector(3, j, eval(w[j]))));
L = Vec(L); S0 = polfromvec(c0);
c1 = normvec(transformvec(S0, MRN)); S1 = polfromvec(c1);
k = 0; for(n = 1, 28, if(MRG[n] != 0, k = c1[n]/MRG[n]; break));
if(k == 0 || c1 != k*MRG || !issquare(MRC*k), error("minred_file.gp: the transformed sextic is not S(N^-1 y) up to a square"));
nl = translines(L, MRN); ns = vector(#NODES, n, prim(MRN*NODES[n]~)~);
for(n = 1, #nl, if(!issqline(S1, nl[n]), error("minred_file.gp: a line of the file is not a square line of the minimised model")));
for(n = 1, #ns, if(vector(3, m, substvec(deriv(S1, [x0,x1,x2][m]), [x0,x1,x2], ns[n])) != [0,0,0], error("minred_file.gp: a transformed singular point is not singular")));
writefile(OUTS, concat(vector(28, n, Str(c1[n], " "))));
my(hdr = Str("# lines on which the minimised model of ", NAME, " (MinRedTernaryForm; y = N x with N = ", MRN, ") is a perfect square: images of the ", #L, " lines of the original file (singular points ", ns, ")\n"));
writefile(OUTL, concat(concat([hdr], vector(#nl, n, Str(nl[n][1], " ", nl[n][2], " ", nl[n][3], "\n")))));
print(NAME, " ", vecmax(abs(c0)), " ", vecmax(abs(c1)), " ", matdet(MRN), " ", vecmax(apply(l -> vecmax(abs(l)), L)), " ", vecmax(apply(l -> vecmax(abs(l)), nl)), " ", ns);
}
\q

/* lifts.gp : for every conic of data/elkies163_conics.txt with a parametrisation x(s,t): the binary sextic g(s,t) with
   S(x(s,t)) = g(s,t)^2 (sign: leading coefficient positive) and a cubic form G(x0,x1,x2) with G(x(s,t)) = g(s,t)
   (unique modulo q times a linear form; a small representative).  The lift C_k of conic k is {q_k = 0, w = G_k(x)},
   its conjugate C_k' is {q_k = 0, w = -G_k(x)}.  Likewise for the 52 lines of data/elkies163_lines.txt: S(tP + sQ) = c(t,s)^2
   (P, Q the primitive basis of src/linesearch.c) and a cubic form G_l with G_l(tP + sQ) = c(t,s).
   Output: out/lifts.gp (CON = [[q, m, hubs, g, Gvec], ...] indexed like the data file, LIN = [[l, PQ, c, Gvec], ...];
   Gvec = coefficients on the monomials MON3).
   usage: cd models && gp -q lifts.gp */
\r cplib.gp
CONICS = readconics();
print(#CONICS, " conics read");
CON = vector(#CONICS);
nnode = 0; nnopar = 0; t0 = getabstime();
{
for(k=1,#CONICS, my(rec = CONICS[k], q = rec[1], m = rec[2], xv, F, g, G);
  if(q[6] == 0, nnode++);
  if(m == 0, nnopar++; CON[k] = [q, 0, rec[3], 0, 0]; next);
  xv = cpar(m);
  if(qform(q, xv) != 0, error("conic ", k, ": parametrisation is not on the conic"));
  F = Sx(xv); g = bsqrt(F); if(g == 0, error("conic ", k, ": S is not a square on the conic"));
  g = possign(g);
  G = fitform(3, xv, g, 2);
  if(substvec(cubic(G[1]), [x0,x1,x2], xv) != g, error("conic ", k, ": lift check failed"));
  CON[k] = [q, m, rec[3], g, G[1]]);
}
print(nnode, " conics through the node (0:0:1), ", nnopar, " without parametrisation; lifts of the other ", #CONICS - nnopar, " computed in ", (getabstime()-t0)/1000., " s");
print("max |coefficient| of the cubic lifts G: ", vecmax(vector(#CON, k, if(CON[k][2] == 0, 0, vecmax(abs(CON[k][5]))))));
print("all cubic lifts integral: ", vector(#CON, k, if(CON[k][2] == 0, 1, denominator(CON[k][5]) == 1)) == vector(#CON, k, 1));
LINES = readlines(basefile("_lines.txt"));
LIN = vector(#LINES);
{
for(k=1,#LINES, my(l = LINES[k], PQ = linebasis(l), xv = lpar(l), F, c, G);
  F = Sx(xv); c = bsqrt(F); if(c == 0, error("line ", l, ": S is not a square"));
  c = possign(c);
  G = fitform(3, xv, c, 1);
  if(substvec(cubic(G[1]), [x0,x1,x2], xv) != c, error("line ", l, ": lift check failed"));
  LIN[k] = [l, PQ, c, G[1]]);
}
print(#LIN, " lines: ", #select(r -> r[1][3] == 0, LIN), " through the node, ", #select(r -> r[1][3] != 0, LIN), " tritangents; cubic lifts computed");
writefile(liftsfile(), Str("\\\\ lifts of the conics and lines of the base model (models/lifts.gp, BASE = \"", BASE, "\"); MON3 = ", MON3, "\nCON = ", CON, ";\nLIN = ", LIN, ";"));
print("written ", liftsfile());
\q

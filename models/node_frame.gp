/* node_frame.gp : a copy of a one-node model in coordinates in which the node is (0:0:1), for use as a base of the
   conic-pair construction (whose scripts assume that position).  x = V y with V unimodular and V (0,0,1)~ = node, so the
   sextic stays integral with the same content, and a line l (row) becomes l V.
   usage (from the root): gp -q with SRC = "<name>" (data/x163_<SRC>.txt, data/x163_<SRC>_lines.txt as written by
   finalize.gp, the singular point in the header of the lines file), DST = "<name>" (data/x163_<DST>.txt, _lines.txt). */
tok(str) = { my(t = strsplit(str, " "), v = List()); for(i = 1, #t, if(#t[i], listput(v, t[i]))); Vec(v) }
MON = []; for(i=0,6, my(ii=6-i); for(jj=0,6-ii, my(j=6-ii-jj, k=6-ii-j); MON = concat(MON, [[ii,j,k]])));
{ my(l = readstr(concat(["data/x163_", SRC, ".txt"])), c = apply(eval, tok(l[1])));
  if(#c != 28, error("28 coefficients expected"));
  S = sum(n=1,28, c[n]*x0^MON[n][1]*x1^MON[n][2]*x2^MON[n][3]); }
{ my(l = readstr(concat(["data/x163_", SRC, "_lines.txt"])), h = l[1], a = strsplit(h, "singular points [["), b);
  if(#a != 2, error("no singular points in the header of the lines file"));
  b = strsplit(a[2], "]]")[1]; if(#strsplit(b, "], [") != 1, error("more than one singular point"));
  node = eval(concat(["[", b, "]"]));
  LNS = List(); for(i = 2, #l, if(#l[i], listput(LNS, apply(eval, tok(l[i]))))); LNS = Vec(LNS); }
if(content(node) != 1, error("node not primitive"));
{ my(g = bezout(node[1], node[2]), h = bezout(g[3], node[3]), u = [g[1]*h[1], g[2]*h[1], h[2]]);     \\ u . node = 1
  if(u*node~ != 1, error("bezout"));
  my(K = matkerint(Mat(u)));                                                                      \\ basis of { x : u.x = 0 }
  K = K*qflll(K);
  V = matconcat([K, node~]); if(abs(matdet(V)) != 1, error("V not unimodular")); }
S2 = substvec(S, [x0,x1,x2], [V[1,1]*x0+V[1,2]*x1+V[1,3]*x2, V[2,1]*x0+V[2,2]*x1+V[2,3]*x2, V[3,1]*x0+V[3,2]*x1+V[3,3]*x2]);
c2 = vector(28, n, polcoef(polcoef(polcoef(S2, MON[n][1], x0), MON[n][2], x1), MON[n][3], x2));
if(substvec(S2, [x0,x1,x2], [0,0,1]) != 0 || vector(3, m, substvec(deriv(S2, [x0,x1,x2][m]), [x0,x1,x2], [0,0,1])) != [0,0,0], error("(0:0:1) is not singular"));
prim(v) = { v = v/content(v); if(v[1] < 0 || (v[1] == 0 && v[2] < 0) || (v[1] == 0 && v[2] == 0 && v[3] < 0), v = -v); v }
L2 = vecsort(vector(#LNS, i, prim(LNS[i]*V)));
{ my(f = fileopen(concat(["data/x163_", DST, ".txt"]), "w")); filewrite(f, concat(vector(28, n, Str(c2[n], " ")))); fileclose(f);
  f = fileopen(concat(["data/x163_", DST, "_lines.txt"]), "w");
  filewrite(f, concat(["# lines on which the sextic data/x163_", DST, ".txt is a perfect square (the model ", SRC, " in coordinates x = V y, V = ", Str(V), ", with the node at (0:0:1))"]));
  for(i = 1, #L2, filewrite(f, Str(L2[i][1], " ", L2[i][2], " ", L2[i][3]))); fileclose(f); }
print(DST, ": node ", node, " of ", SRC, " moved to (0:0:1); max |coeff| ", vecmax(abs(c2)), "; ", #L2, " lines, ", sum(i = 1, #L2, L2[i][3] == 0), " through the node");
\q

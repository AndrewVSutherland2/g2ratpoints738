/* config.gp : for each polarization D' of a dp_<B>.txt file, the configuration of its contracted curves and lines:
     one of the two roots +-r for each pair of roots of D'^perp (v.D' = 0, v^2 = -2), and
     the classes v with v.D' = 1, v^2 = -2 and v.r >= 0 for the roots chosen.
   When the root system is A1^k the chosen roots are, up to the reflections in them (isometries of NS that fix D'), the
   contracted curves, and the classes kept are then the components of the pullbacks of the tritangents and node lines;
   another choice of signs gives an isomorphic configuration.  All of these are among the vectors with
   Phi(v) = -v^2 + (v.D')^2 <= 3.  Printed on stdout, one row per D':
     "id n2 index nroots nlines | the vectors, 20 coordinates each, roots first"
   where index is the index in NS of the lattice they span (1 when the lines and the contracted curves generate NS;
   then the matrix of intersection numbers of these classes, up to a permutation, determines the model, see canon.c).
   usage: gp -q with LAT = "lat_<B>.txt", DP = "dp_<B>.txt" (or a piece of it). */
default(parisizemax, 8000000000);
tok(str) = { my(t = strsplit(str, " "), v = List()); for(i = 1, #t, if(#t[i], listput(v, t[i]))); Vec(v) }
{ my(l = readstr(LAT)); for(i = 1, #l, my(w = tok(l[i])); if(w[1] == "GZ", GZ = matrix(20, 20, a, b, eval(w[1 + 20*(a-1) + b])); break)); }
cfg(d) = { my(gd = GZ*d, Phi = -GZ + gd*gd~, V = qfminim(Phi, 3)[3], n = #V, tv = gd~*V, R = List(), L = List());
  for(k = 1, n, my(t = tv[k]); if(t == 0, listput(R, V[,k])));
  my(GR = if(#R, (GZ*matconcat(Vec(R)))~, matrix(0, 20)));
  for(k = 1, n, my(t = tv[k], v = V[,k]); if(abs(t) != 1, next); v = t*v;
    my(w = GR*v, ok = 1); for(j = 1, #w, if(w[j] < 0, ok = 0; break)); if(ok, listput(L, v)));
  my(M = matconcat(concat(Vec(R), Vec(L))), H = mathnf(M));
  [#R, #L, if(#H == 20, abs(matdet(H)), 0), M] }
{ my(l = readstr(DP));
  for(i = 1, #l, my(p = strsplit(l[i], " | "), h = tok(p[1]), d = apply(eval, tok(p[2]))~, c = cfg(d));
    if(c[1] != eval(h[2]), error("number of roots: ", c[1], " for ", h[1], ", expected ", h[2]));
    print(h[1], " ", h[2], " ", c[3], " ", c[1], " ", c[2], " | ", concat(vector(#c[4], k, concat(vector(20, a, Str(c[4][a, k], " "))))))); }
\q

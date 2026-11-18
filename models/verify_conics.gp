\\ verify_conics.gp -- exact verification and rational parametrisation of the conics found by conicmodp + lift_conics.py
\\
\\   echo 'read("models/verify_conics.gp"); main("data/elkies163.txt", "lifted.txt", "conics_param.txt")' | gp -q
\\
\\ Input lines: "c00 c01 c02 c11 c12 c22 | g0 .. g9 | h0 .. h14" (q integral, g cubic, h quartic with rational
\\ coefficients, monomials in the order i = d..0, j = d-i..0 for x0^i x1^j x2^k).  For each conic the script checks
\\ (1) S - g^2 - q*h == 0 exactly (so S is the square of g on the conic, whether or not the conic has rational points),
\\ (2) qfsolve: does the conic have a rational point?  If so, qfparam gives a parametrisation (x0,x1,x2) = M*[s^2,st,t^2]~
\\     which is made primitive and reduced (SL2(Z) reduction of the binary quadratic forms), and we check that
\\     q(M(s,t)) == 0 and that S(M(s,t)) is the square of a binary sextic (issquare).
\\ Output lines: "c00 c01 c02 c11 c12 c22  m00 m01 m02 m10 m11 m12 m20 m21 m22" (row i = coefficients of s^2, st, t^2
\\ of x_i) or "... none" when the conic has no rational point.

readsextic(f) = { my(c = apply(eval, strsplit(readstr(f)[1], " ")), idx = 1, S = 0); forstep(i=6,0,-1, forstep(j=6-i,0,-1, my(k=6-i-j); S += c[idx]*x0^i*x1^j*x2^k; idx++)); S; }
monos(d) = { my(v = []); forstep(i=d,0,-1, forstep(j=d-i,0,-1, v = concat(v, [x0^i*x1^j*x2^(d-i-j)]))); v; }
ev(F, v) = subst(subst(subst(F, x0, v[1]), x1, v[2]), x2, v[3]);

\\ reduce a parametrisation matrix M (3x3, rows = binary quadratic forms in (s,t)) by SL2(Z) acting on (s,t), to
\\ minimise the sum of squares of the entries (simple greedy descent with s -> s + k t, t -> t + k s, swap)
sizeM(M) = sum(i=1,3, sum(j=1,3, M[i,j]^2));
shiftS(M, k) = { my(N = M); for(i=1,3, my(a=M[i,1],b=M[i,2],c=M[i,3]); \\ form a s^2 + b s t + c t^2 with s -> s + k t
  N[i,1] = a; N[i,2] = 2*a*k + b; N[i,3] = a*k^2 + b*k + c); N; }
swapST(M) = { my(N = M); for(i=1,3, N[i,1] = M[i,3]; N[i,3] = M[i,1]); N; }
reduceM(M) = {
  my(best = M, improved = 1);
  while(improved, improved = 0;
    my(cur = sizeM(best));
    foreach([1,-1], k, my(N = shiftS(best, k)); if(sizeM(N) < cur, best = N; cur = sizeM(N); improved = 1));
    my(N = swapST(best)); foreach([1,-1], k, my(N2 = shiftS(N, k)); if(sizeM(N2) < cur, best = N2; cur = sizeM(N2); improved = 1)));
  best;
}

main(SEXTIC, IN, OUT) = {
  my(S = readsextic(SEXTIC), L = readstr(IN), M3 = monos(3), M4 = monos(4), M2 = monos(2));
  my(nid = 0, npt = 0, nnone = 0, nbad = 0, fo = fileopen(OUT, "w"));
  for(n=1,#L, my(parts = strsplit(L[n], "|")); if(#parts < 3, next);
    my(qc = apply(eval, select(w->#w, strsplit(parts[1], " "))));
    my(gc = apply(eval, select(w->#w, strsplit(parts[2], " "))));
    my(hc = apply(eval, select(w->#w, strsplit(parts[3], " "))));
    my(q = sum(i=1,6, qc[i]*M2[i]), g = sum(i=1,10, gc[i]*M3[i]), h = sum(i=1,15, hc[i]*M4[i]));
    if(S - g^2 - q*h != 0, print("IDENTITY FAILS for conic ", qc); nbad++; next);
    nid++;
    my(G = [2*qc[1], qc[2], qc[3]; qc[2], 2*qc[4], qc[5]; qc[3], qc[5], 2*qc[6]]);
    my(sol = qfsolve(G));
    if(type(sol) != "t_COL",
      nnone++; filewrite(fo, concat(concat(strjoin(apply(x->Str(x), qc), " "), " "), "none")); next);
    my(M = qfparam(G, sol));
    my(cont = 0); for(i=1,3, for(j=1,3, cont = gcd(cont, M[i,j]))); M = M/cont;
    M = reduceM(M);
    my(X = M*[s^2, s*t, t^2]~);
    if(ev(q, X) != 0, print("PARAM FAILS q(M) != 0 for ", qc); nbad++; next);
    my(F = ev(S, X));
    if(!issquare(F), print("PARAM FAILS S(M) not a square for ", qc); nbad++; next);
    \\ the parametrisation must be a bijection: the three forms must not have a common factor (content already removed)
    npt++;
    filewrite(fo, concat(concat(strjoin(apply(x->Str(x), qc), " "), " "), strjoin(apply(x->Str(x), [M[1,1],M[1,2],M[1,3],M[2,1],M[2,2],M[2,3],M[3,1],M[3,2],M[3,3]]), " "))));
  fileclose(fo);
  print("identity S = g^2 + q h verified for ", nid, " conics; with rational point (parametrised, S(M) square checked): ", npt, "; without rational points: ", nnone, "; failures: ", nbad);
}

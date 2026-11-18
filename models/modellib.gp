/* modellib.gp : from a basis f = [f0,f1,f2] of L(D') (each f_i = [A_i, B_i], function (A_i + y B_i)/x2^N) compute the
   branch sextic S_new with y_new = w := y * det(f, d_a f, d_b f) (affine chart a = x0/x2, b = x1/x2, y_aff = y/x2^3),
   w^2 = S_new(f0, f1, f2). */
SEaff = substvec(SE, [x0,x1,x2], [a, b, 1]);
SEa = deriv(SEaff, a); SEb = deriv(SEaff, b);
\\ affine polynomial pair [P, Q] for f = [A,B]: P + y Q with P = A(a,b,1), Q = B(a,b,1)
aff(f) = [substvec(f[1],[x0,x1,x2],[a,b,1]), substvec(f[2],[x0,x1,x2],[a,b,1])];
\\ derivative w.r.t. a of P + yQ on X: (P_a + Q S_a/(2S) * ... ) : d/da (P + yQ) = P_a + y Q_a + Q S_a/(2y) = P_a + y (Q_a + Q S_a/(2S))
\\ we keep functions as [P, Q] meaning P + y Q with P, Q rational functions in a, b
dA(g) = [deriv(g[1],a), deriv(g[2],a) + g[2]*SEa/(2*SEaff)];
dB(g) = [deriv(g[1],b), deriv(g[2],b) + g[2]*SEb/(2*SEaff)];
gmul(g,h) = [g[1]*h[1] + g[2]*h[2]*SEaff, g[1]*h[2] + g[2]*h[1]];
gsub(g,h) = [g[1]-h[1], g[2]-h[2]];
gadd(g,h) = [g[1]+h[1], g[2]+h[2]];
\\ evaluate [P,Q] at affine point (a0,b0,y0)
geval(g, pt) = substvec(g[1],[a,b],[pt[1],pt[2]]) + pt[3]*substvec(g[2],[a,b],[pt[1],pt[2]]);
\\ symbolic w = y * det(f, f_a, f_b) as [P,Q]
wsym(f) = { my(g = vector(3,i,aff(f[i])), ga = vector(3,i,dA(g[i])), gb = vector(3,i,dB(g[i])), det = [0,0]);
  \\ det of 3x3 with rows g, ga, gb (entries [P,Q])
  my(m = [g, ga, gb], perms = [[1,2,3],[2,3,1],[3,1,2]], nperms = [[1,3,2],[3,2,1],[2,1,3]]);
  for(t=1,3, my(p = perms[t]); det = gadd(det, gmul(gmul(m[1][p[1]], m[2][p[2]]), m[3][p[3]])));
  for(t=1,3, my(p = nperms[t]); det = gsub(det, gmul(gmul(m[1][p[1]], m[2][p[2]]), m[3][p[3]])));
  gmul([0,1], det) }
\\ affine point from a projective point [x, y] (weights 1,3): (x0/x2, x1/x2, y/x2^3)
affpt(P) = [P[1][1]/P[1][3], P[1][2]/P[1][3], P[2]/P[1][3]^3];
\\ fit S_new: w^2 = sum_m a_m m(f0,f1,f2) at points pts (affine), returns coefficient vector in MON6 order
fitsextic(f, w, pts) = { my(g = vector(3,i,aff(f[i])), rows = List(), rhs = List());
  for(t=1,#pts, my(pt = pts[t], fv = vector(3,i,geval(g[i],pt)), wv = geval(w, pt));
    listput(rows, vector(28, n, monval(MON6[n], fv))); listput(rhs, wv^2));
  my(M = matconcat(Vec(rows)~), r = Vec(rhs)~);
  if(matrank(M) < 28, error("fit: rank ", matrank(M)));
  my(sol = matsolve(M~*M, M~*r));   \\ least squares = exact solution if consistent
  if(M*sol != r, error("fit: inconsistent")); sol~ }
\\ normalise: multiply by a square to get integer coefficients, remove square content
normsq(c) = { my(d = 1, cf, m, q); for(n=1,#c, d = lcm(d, denominator(c[n]))); c = c*d^2; q = core(content(c),1)[2]; c/q^2 }
sextic(c) = sum(n=1,28, c[n]*x0^MON6[n][1]*x1^MON6[n][2]*x2^MON6[n][3]);
\\ collect evaluation points: on lifts k (avoiding lines x2 = 0) at parameters v; skip points with y = 0 or f0 = 0 or bad denominators
evalpoints(f, lifts, vs) = { my(pts = List(), g = aff(f[1]));
  for(t=1,#lifts, my(k = lifts[t]); for(s=1,#vs, my(P = ptonlift(k, vs[s])); if(P[1][3]==0 || P[2]==0, next);
    my(pt = affpt(P)); if(substvec(SEaff,[a,b],[pt[1],pt[2]]) == 0, next); if(geval(g, pt) == 0, next); listput(pts, pt)));
  Vec(pts) }
restrict(S, l) = if(l[3]!=0, substvec(S,[x0,x1,x2],[u,v,-(l[1]*u+l[2]*v)/l[3]]), if(l[2]!=0, substvec(S,[x0,x1,x2],[u,-l[1]*u/l[2],v]), substvec(S,[x0,x1,x2],[0,u,v])));

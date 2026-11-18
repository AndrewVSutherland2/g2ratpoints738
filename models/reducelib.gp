/* reducelib.gp : coefficient reduction of a ternary sextic by projective coordinate changes.
   c = coefficient vector (28, linesearch order); y = M x new coordinates; S_new(y) = S(M^-1 y). */
MON6 = []; for(i=0,6, my(ii=6-i); for(jj=0,6-ii, my(j=6-ii-jj, k=6-ii-j); MON6 = concat(MON6, [[ii,j,k]])));
coefvec(Sp) = vector(28, n, polcoeff(polcoeff(polcoeff(Sp, MON6[n][1], x0), MON6[n][2], x1), MON6[n][3], x2));
polfromvec(c) = sum(n=1,28, c[n]*x0^MON6[n][1]*x1^MON6[n][2]*x2^MON6[n][3]);
transformvec(Sp, M) = { my(Mi = M^-1, ff = vector(3, i, Mi[i,1]*x0 + Mi[i,2]*x1 + Mi[i,3]*x2)); coefvec(substvec(Sp, [x0,x1,x2], ff)) }
normvec(c) = { my(D = 1, m, ct, q); for(n=1,28, D = lcm(D, denominator(c[n]))); my(cf = core(D,1)); m = cf[1]*cf[2]; c = c*m^2; ct = content(c); q = core(ct,1)[2]; c/q^2 }
sizev(c) = { my(mx=0, sm=0); for(n=1,28, if(c[n]!=0, my(l = log(abs(c[n]))/log(10)); mx = max(mx,l); sm += l)); [mx, sm] }
scorev(c) = { my(s = sizev(c)); s[1] + s[2]/28/10 }
torus(c, lam) = vector(28, n, c[n]*lam[1]^MON6[n][1]*lam[2]^MON6[n][2]*lam[3]^MON6[n][3]);
PR = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47];
relprimes(c, bound=5000) = { my(res=List(PR)); forprime(p=53, bound, my(cnt=0); for(n=1,28, if(c[n]!=0 && valuation(c[n],p)>0, cnt++)); if(cnt>=6, listput(res,p))); Vec(res) }
torussearch(c, primes) = { my(lam=[1,1,1], best = normvec(c), bs = scorev(best), improved=1);
  while(improved, improved=0;
    for(i=1,3, for(pi=1,#primes, for(e=-1,1, if(e==0,next);
      my(l2=lam); l2[i] *= primes[pi]^e; my(c2 = normvec(torus(c,l2)), s2 = scorev(c2));
      if(s2 < bs - 1e-9, bs=s2; best=c2; lam=l2; improved=1)))));
  [best, lam, bs] }
QS = [1,-1,2,-2,3,-3,1/2,-1/2,4,-4,1/3,-1/3,5,-5,1/4,-1/4,3/2,-3/2,2/3,-2/3,6,-6,1/5,-1/5];
shearmat(i,j,q) = { my(E = matid(3)); E[i,j] = q; E }
diagmat(i,p) = { my(E = matid(3)); E[i,i] = p; E }
applyM(c, E) = normvec(transformvec(polfromvec(c), E));
\\ the two kinds of moves of localsearch on the coefficient vector itself (the same vectors as applyM gives, without a
\\ polynomial substitution): y = E x with E = diagmat(i, p) divides x_i by p; E = shearmat(i, j, q) replaces x_i by x_i - q x_j
MONIDX = Map(); for(n=1,28, mapput(MONIDX, MON6[n], n));
SHT = vector(3, i, vector(3, j, if(i == j, 0, vector(28, n, my(a = MON6[n][i]); vector(a+1, t, my(e = MON6[n]); e[i] -= t-1; e[j] += t-1; [mapget(MONIDX, e), binomial(a, t-1)])))));
shearvec(c, i, j, q) = { my(r = vector(28), pw = vector(7, t, (-q)^(t-1)), T = SHT[i][j]);
  for(n=1,28, if(c[n] == 0, next); my(L = T[n], cn = c[n]); for(t=1,#L, r[L[t][1]] += cn*L[t][2]*pw[t])); r }
diagvec(c, i, p) = vector(28, n, c[n]/p^MON6[n][i]);
localsearch(c, M, primes, maxit=60, verbose=0) = { my(bs, improved=1, it=0);
  c = normvec(c); bs = scorev(c);
  while(improved && it<maxit, improved=0; it++;
    my(bestmove=0, bestc, bestsc=bs);
    for(i=1,3, for(pi=1,#primes, for(e=-1,1, if(e==0,next);
      my(E=diagmat(i,primes[pi]^e), c2=normvec(diagvec(c,i,primes[pi]^e)), s2=scorev(c2));
      if(s2 < bestsc-1e-9, bestsc=s2; bestc=c2; bestmove=E))));
    for(i=1,3, for(j=1,2, if(i==j,next); for(qi=1,#QS,
      my(E=shearmat(i,j,QS[qi]), c2=normvec(shearvec(c,i,j,QS[qi])), s2=scorev(c2));
      if(s2 < bestsc-1e-9, bestsc=s2; bestc=c2; bestmove=E))));
    if(bestmove!=0, c=bestc; bs=bestsc; M=bestmove*M; improved=1; if(verbose, print("it ",it," size ", sizev(c)))));
  [c, M, bs] }
\\ frame from three lines (rows = the three linear forms): y = M x with y_i = l_i . x
framemat(l1, l2, l3) = matconcat([l1; l2; l3]);
\\ transform lines: a line l (row) in x-coordinates becomes l*M^-1 in y-coordinates
prim(v) = { my(d=1); for(j=1,#v, d=lcm(d,denominator(v[j]))); v=v*d; v=v/content(v); if(v[1]<0 || (v[1]==0 && v[2]<0) || (v[1]==0&&v[2]==0&&v[3]<0), v=-v); v };
translines(L, M) = vector(#L, i, prim(L[i]*M^-1));

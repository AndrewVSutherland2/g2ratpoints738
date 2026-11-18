/* cplib.gp : common definitions for models (Elkies' double plane w^2 = S(x0,x1,x2), its conics and lines).
   Coefficient order of a ternary sextic: the one of src/linesearch.c (for i = 6..0, j = 6-i..0, k = 6-i-j: x0^i x1^j x2^k).
   Variables: x0,x1,x2 (plane), s,t (parameters of a conic or a line), u,v (used by modellib/linesq), a,b,y (affine chart).
   All scripts are run from this directory (REPO = "..": the directory above holds data/ and src/). */
\\ (parisizemax is set by run.sh: a default(parisizemax) inside a file loaded with \r aborts the rest of that file)
REPO =if(type(REPO) == "t_STR", REPO, "..");            \\ path of the repository root (override with REPO = "..." before \r)
BASE =if(type(BASE) == "t_STR", BASE, "");                 \\ "" = Elkies' model as the base surface; "<M>" = data/x163_<M>.txt with its _conics/_lines tables (frontier B)
basefile(suf) = if(BASE == "", concat(REPO, concat("/data/elkies163", suf)), concat([REPO, "/data/x163_", BASE, suf]));
liftsfile() = if(BASE == "", "out/lifts.gp", concat(["out/lifts_", BASE, ".gp"]));
MON6 = []; for(i=0,6, my(ii=6-i); for(jj=0,6-ii, my(j=6-ii-jj, k=6-ii-j); MON6 = concat(MON6, [[ii,j,k]])));
mons(d) = { my(L=List()); for(i=0,d, for(j=0,d-i, listput(L, [d-i-j,j,i]))); Vec(L) }   \\ any fixed order; MON6 is linesearch's
monval(m, x) = x[1]^m[1]*x[2]^m[2]*x[3]^m[3];
sextic(c) = sum(n=1,28, c[n]*x0^MON6[n][1]*x1^MON6[n][2]*x2^MON6[n][3]);
coefvec(Sp) = vector(28, n, polcoeff(polcoeff(polcoeff(Sp, MON6[n][1], x0), MON6[n][2], x1), MON6[n][3], x2));
\\ parse a line of blank-separated integers (the token 'none' is kept as a string)
parseline(str) = { my(tk = strsplit(str, " "), L = List()); for(i=1,#tk, if(#tk[i] == 0, next); my(w = tk[i]); if(w == "none", listput(L, w), listput(L, eval(w)))); Vec(L) }
readints(file) = { my(ls = readstr(file), L = List()); for(i=1,#ls, if(#ls[i] == 0 || Vec(ls[i])[1] == "#", next); my(v = parseline(ls[i])); for(j=1,#v, listput(L, v[j]))); Vec(L) }
cE = readints(basefile(".txt"));
if(#cE != 28, error(basefile(".txt"), ": expected 28 coefficients"));
SE = sextic(cE);
Sx(xv) = substvec(SE, [x0,x1,x2], xv);
NODE = [0,0,1];
\\ conics: records [q (6 coeffs c00 c01 c02 c11 c12 c22), m (3x3 matrix, rows x_i = m_i0 s^2 + m_i1 s t + m_i2 t^2) or 0, hubs]
readconics() = { my(ls = readstr(basefile("_conics.txt")), L = List());
  for(i=1,#ls, if(#ls[i] == 0 || Vec(ls[i])[1] == "#", next); my(v = parseline(ls[i]), q = v[1..6], m, h);
    if(v[7] == "none", m = 0; h = v[8], m = matrix(3,3,r,c, v[6 + 3*(r-1) + c]); h = v[16]);
    listput(L, [q, m, h]));
  Vec(L) }
qform(q, xv) = q[1]*xv[1]^2 + q[2]*xv[1]*xv[2] + q[3]*xv[1]*xv[3] + q[4]*xv[2]^2 + q[5]*xv[2]*xv[3] + q[6]*xv[3]^2;
qpol(q) = qform(q, [x0,x1,x2]);
\\ parametrisation as a vector of binary quadratic forms in s,t
cpar(m) = vector(3, i, m[i,1]*s^2 + m[i,2]*s*t + m[i,3]*t^2);
\\ lines: "A B C" records
readlines(file) = { my(ls = readstr(file), L = List()); for(i=1,#ls, if(#ls[i] == 0 || Vec(ls[i])[1] == "#", next); listput(L, parseline(ls[i]))); Vec(L) }
\\ primitive basis P, Q of the line A x0 + B x1 + C x2 = 0, exactly as src/linesearch.c line_basis (x = t P + s Q)
linebasis(l) = { my(A = l[1], B = l[2], C = l[3], g, e, uu, vv, P, Q);
  if(B == 0 && C == 0, return([[0,1,0], [0,0,1]]));
  g = gcd(B, C); Q = [0, C/g, -B/g];
  e = bezout(abs(B/g), abs(C/g)); uu = e[1]; vv = e[2];
  if(B < 0, uu = -uu); if(C < 0, vv = -vv);
  P = [g, -A*uu, -A*vv];
  if(A*P[1]+B*P[2]+C*P[3] != 0 || A*Q[1]+B*Q[2]+C*Q[3] != 0, error("linebasis"));
  [P, Q] }
lpar(l) = { my(PQ = linebasis(l)); vector(3, i, t*PQ[1][i] + s*PQ[2][i]) }
\\ binary forms in s,t: exact square root or 0
bsqrt(F) = { my(r); if(F == 0, return(0)); if(!issquare(F, &r), return(0)); r }
\\ normalise the sign of a binary form: leading coefficient in s (then in t) positive
possign(F) = { my(c = pollead(F, s)); if(type(c) != "t_INT" && type(c) != "t_FRAC", c = pollead(c, t)); if(c < 0, -F, F) }
\\ degree of the gcd of two binary forms A, B (in s,t) of degrees dA, dB: common roots, (1:0) included
bgcddeg(A, B, dA, dB) = { my(a = subst(A, t, 1), bb = subst(B, t, 1), d0);
  if(a == 0 || bb == 0, error("bgcddeg: zero form"));
  d0 = poldegree(gcd(a, bb), s); if(d0 < 0, d0 = 0);
  d0 + min(dA - poldegree(a, s), dB - poldegree(bb, s)) }
\\ cubic forms: 10 monomials
MON3 = mons(3);
cubic(v) = sum(n=1,10, v[n]*monval(MON3[n], [x0,x1,x2]));
\\ solve for a form of degree d in x0,x1,x2 with prescribed restriction F (binary form in s,t of degree d*e) to the
\\ parametrised curve xv (vector of 3 binary forms of degree e); returns a small solution (Babai rounding against the kernel)
fitform(d, xv, F, e) = { my(M = mons(d), n = #M, de = d*e, rows, A, rhs, sol, K);
  rows = vector(n, j, my(m = monval(M[j], xv)); vector(de+1, i, polcoeff(polcoeff(m, de+1-i, s), i-1, t)));
  A = matconcat(rows~)~;                      \\ (de+1) x n: column j = coefficients of the restricted monomial j
  rhs = vector(de+1, i, polcoeff(polcoeff(F, de+1-i, s), i-1, t))~;
  sol = matinverseimage(A, rhs); if(sol == 0 && rhs != 0, error("fitform: no solution"));
  K = matker(A, 1);                           \\ integral kernel basis
  if(#K > 0, K = K*qflll(K); my(cf = round(matsolve(K~*K, K~*sol))); sol = sol - K*cf);
  [sol~, M] }
formfromvec(sol, M) = sum(j=1,#M, sol[j]*monval(M[j], [x0,x1,x2]));
\\ write a string to a file, truncating it (gp's write() appends, and a second run would then produce a double file)
writefile(name, str) = { my(F = fileopen(name, "w")); filewrite(F, str); fileclose(F) }
\\ restriction of a ternary form to the line l = [A,B,C] (binary form in u,v; as in linesq.gp)
restrict(Sp, l) = if(l[3]!=0, substvec(Sp,[x0,x1,x2],[u,v,-(l[1]*u+l[2]*v)/l[3]]), if(l[2]!=0, substvec(Sp,[x0,x1,x2],[u,-l[1]*u/l[2],v]), substvec(Sp,[x0,x1,x2],[0,u,v])));
issqline(Sp, l) = { my(F = numerator(restrict(Sp, l))); F != 0 && issquare(F) }
gradat(Sp, P) = vector(3, m, substvec(deriv(Sp, [x0,x1,x2][m]), [x0,x1,x2], P));
prim(v) = { my(d = 1); for(n=1,#v, d = lcm(d, denominator(v[n]))); v = v*d; v = v/content(v); if(v[1] < 0 || (v[1] == 0 && v[2] < 0) || (v[1] == 0 && v[2] == 0 && v[3] < 0), v = -v); v }
cross(l, m) = [l[2]*m[3]-l[3]*m[2], l[3]*m[1]-l[1]*m[3], l[1]*m[2]-l[2]*m[1]];

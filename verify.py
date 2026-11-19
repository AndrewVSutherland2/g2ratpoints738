#!/usr/bin/env python3
"""verify.py -- machine check of the formulas displayed in curve738.tex against the data files in data/.
Extracts the displayed polynomials from the LaTeX source, converts them to PARI/GP syntax and checks every identity
stated in the paper (python3 verify.py from anywhere; needs gp).  The check runs on this directory alone: it reads curve738.tex
and data/ next to this script."""
import re, subprocess, sys, os, tempfile
N_CHECKS = 46     # the number of chk(...) calls in the GP script below: the verdict requires every one of them to have run and passed
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, 'data')
tex = open(os.path.join(HERE, 'curve738.tex')).read()

def tex2gp(e):
    e = e.replace('\\\\', ' ').replace('&', ' ').replace('{}', ' ')
    for t in ('\\bigl', '\\bigr', '\\,', '\\;', '\\quad', '\\qquad', '\\colon', '\\Bigl', '\\Bigr', '\\ '): e = e.replace(t, ' ')
    e = re.sub(r'\\tfrac\{(-?\d+)\}\{(\d+)\}', r'(\1/\2)', e)
    e = re.sub(r'\\frac\{(-?\d+)\}\{(\d+)\}', r'(\1/\2)', e)
    e = re.sub(r'([xz])_(\d)', r'\1\2', e)
    e = re.sub(r'\^\{(\d+)\}', r'^\1', e)
    e = e.strip().rstrip('.,').strip()
    toks = re.findall(r'\d+|[xz]\d|[xwyp]|[-+*/^()]', e)
    assert ''.join(toks) == re.sub(r'\s+', '', e), ('untokenised input', e, ''.join(toks))
    out = []
    for t in toks:
        if out:
            a = out[-1]
            left = a[-1].isalnum() or a == ')'
            right = t[0].isalnum() or t == '('
            if left and right and not (out[-1] == '^'): 
                # exponent digits directly after ^ are part of the power: "x^2x" -> x^2*x is still a product
                out.append('*')
        out.append(t)
    return ''.join(out)

def grab(start, end, which=0):
    i = -1
    for _ in range(which + 1):
        i = tex.index(start, i + 1)
    j = tex.index(end, i)
    return tex[i + len(start):j]

F = {}
F['f738'] = tex2gp(grab('C_{738}\\colon\\quad y^2 ={}&', '\\end{align*}'))
F['f642'] = tex2gp(grab('C_{642}\\colon\\quad y^2 ={}&', '\\end{align*}'))
F['f654'] = tex2gp(grab('C_{654}\\colon\\quad y^2={}&', '\\end{align*}'))
F['SE'] = tex2gp(grab('S_E={}&', '\\end{align*}'))
F['f648'] = tex2gp(grab('y^2 = 98285087781', '\\]').join(['98285087781', '']) if False else '98285087781' + grab('y^2 = 98285087781', '\\]'))
blk = grab('q_1&=', '\\end{align*}'); q1t, q2t = blk.split('q_2&=')
F['q1'] = tex2gp(q1t); F['q2'] = tex2gp(q2t)
blk = grab('G_1&=', '\\end{align*}'); g1t, g2t = blk.split('G_2&=')
F['G1'] = tex2gp(g1t); F['G2'] = tex2gp(g2t)
blk = grab('B_0=', '\\]'); b = re.split(r'B_[12]=', blk)
for k in range(3): F['B%d' % k] = tex2gp(b[k].replace('\\quad', ' '))
blk = grab('A_0={}&', '\\end{align*}'); a = re.split(r'A_[12]=\{\}&', blk)
for k in range(3): F['A%d' % k] = tex2gp(a[k])
F['S3'] = tex2gp(grab('S_{3}={}&', '\\end{align*}'))
F['par'] = tex2gp(grab('S_3(6x-2,\\,4x+16,\\,10x+7)=64\\,\\bigl(&', '\\bigr)'))
F['p'] = tex2gp(grab('p(x)=', '\\]'))
F['res'] = tex2gp(grab('p(x)^2-f(x)=', '\\]', 1))
gp = ['default(parisizemax, 2000000000);', 'w; x0; x1; x2; z0; z1; z2;   \\\\ fixes the order of the variables: w is the main variable of the identity for the map', 'ok = 1; chk(name, b) = { print(if(b, "ok    ", "FAIL  "), name); if(!b, ok = 0); }']
for k, v in F.items(): gp.append('%s = %s;' % (k, v))
gp.append(r'''
MON6 = []; for(i=0,6, my(ii=6-i); for(jj=0,6-ii, my(j=6-ii-jj, k=6-ii-j); MON6 = concat(MON6, [[ii,j,k]])));
nums(file) = { my(R = readstr(file), L = List()); for(i=1,#R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = strsplit(R[i], " ")); for(j=1,#w, if(#w[j], listput(L, eval(w[j]))))); Vec(L) }
rows(file) = { my(R = readstr(file), L = List()); for(i=1,#R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = select(t -> #t > 0, strsplit(R[i], " "))); listput(L, apply(eval, w))); Vec(L) }
sx(c) = sum(n=1,28, c[n]*x0^MON6[n][1]*x1^MON6[n][2]*x2^MON6[n][3]);
ev(P, v) = substvec(P, [x0,x1,x2], v);
\\ 1. the equations
chk("f738 in the abstract = f738 in the parametrisation display", f738 == par);
chk("S_E = data/elkies163.txt, the input of the computation", SE == sx(nums(Str(ROOT, "/elkies163.txt"))));
chk("a = 53360 = 640320/12, and the coefficients of the Inose equation: 3 a^2 = 8541868800, 2 a^3 + 2 = 303862746112002", 12*53360 == 640320 && 3*53360^2 == 8541868800 && 2*53360^3 + 2 == 303862746112002);
chk("S_3 = data/S3_coefficients.txt", subst(subst(subst(S3, z0, x0), z1, x1), z2, x2) == sx(nums(Str(ROOT, "/S3_coefficients.txt"))));
S3x = subst(subst(subst(S3, z0, x0), z1, x1), z2, x2);
CON = rows(Str(ROOT, "/elkies163_conics.txt"));
chk("the table of conics of Elkies' sextic has 1241 rows, and the table of lines 52: 43 tritangents and 9 lines through the node", #CON == 1241 && #rows(Str(ROOT, "/elkies163_lines.txt")) == 52 && #select(l -> l[3] == 0, rows(Str(ROOT, "/elkies163_lines.txt"))) == 9);
qc(v) = v[1]*x0^2 + v[2]*x0*x1 + v[3]*x0*x2 + v[4]*x1^2 + v[5]*x1*x2 + v[6]*x2^2;
par(v) = vector(3, i, v[6+3*(i-1)+1]*s^2 + v[6+3*(i-1)+2]*s + v[6+3*(i-1)+3]);
chk("q1, q2 = conics 410 and 424 of the table", q1 == qc(CON[410]) && q2 == qc(CON[424]));
X1 = par(CON[410]); X2 = par(CON[424]);
chk("S_E = G1^2 on q1 and G2^2 on q2", ev(q1, X1) == 0 && ev(q2, X2) == 0 && ev(SE, X1) == ev(G1, X1)^2 && ev(SE, X2) == ev(G2, X2)^2);
\\ 2. the lifts agree at three of the four intersection points
{ my(h = ev(q2, X1), d = ev(G1, X1) - ev(G2, X1), a = ev(G1, X1) + ev(G2, X1), fa = factor(h));
  chk("q1, q2 meet in a rational point and a cubic point", matsize(fa)[1] == 2 && vecsort(vector(2, i, poldegree(fa[i,1]))) == [1, 3]);
  chk("Gamma1.Gamma2 = 3, Gamma1.Gamma2' = 1", poldegree(gcd(h, d)) == 3 && poldegree(gcd(h, a)) == 1); }
\\ 3. the map to the new model: S_3(A + B w) is a rational square at rational points of X
A = [A0, A1, A2]; B = [B0, B1, B2];
zmap(P, wv) = vector(3, m, ev(A[m], P) + wv*ev(B[m], P));
{ my(n = 0, bad = 0, sq(v) = issquare(v));
  forstep(k = 3, 1200, 7, my(c = CON[k]); if(#c < 15 || type(c[7]) != "t_INT", next); my(X = par(c));
    foreach([2, -3, 1/2, 5, -7/3], sv, my(P = subst(X, s, sv), v = ev(SE, P), wv); if(!issquare(v, &wv) || wv == 0, next);
      foreach([1, -1], e, my(z = zmap(P, e*wv)); if(z == [0,0,0], next); n++; if(!sq(ev(S3x, z)), bad++))));
  chk(Str("S_3(A + Bw) is a square at ", n, " rational points of X on conics of Elkies' model"), n > 300 && bad == 0); }
\\ the map as an exact identity, in the chart x0 = 1: with z_i = A_i + B_i w and T = w det(z, dz/dx1, dz/dx2), where the
\\ derivatives are taken on the surface (dw/dx_j = (dS_E/dx_j)/(2w)), S_3(z) = T^2/16 modulo w^2 - S_E.  Below c_D = w T.
c_s1(P) = subst(P, x0, 1);
c_SE = c_s1(SE); c_A = apply(c_s1, A); c_B = apply(c_s1, B); c_Z = vector(3, i, c_A[i] + c_B[i]*w);
c_N(v) = vector(3, i, w*deriv(c_A[i], v) + c_SE*deriv(c_B[i], v) + c_B[i]*deriv(c_SE, v)/2);
c_D = matdet(matrix(3, 3, r, c, if(r == 1, c_Z[c], if(r == 2, c_N(x1)[c], c_N(x2)[c]))));
chk("S_3(z) = T^2/16 modulo w^2 - S_E, an identity of polynomials (T = w det(z, dz/dx1, dz/dx2), chart x0 = 1)", variables(w*x1*x2)[1] == w && (16*c_SE*substvec(S3, [z0, z1, z2], c_Z) - c_D^2) % (w^2 - c_SE) == 0);
chk("A_j = B_j G_i modulo q_i for both conics and all three coordinates", vector(2, i, my(q = c_s1([q1, q2][i]), G = c_s1([G1, G2][i])); vector(3, j, (c_A[j] - c_B[j]*G) % q == 0)) == [[1,1,1],[1,1,1]]);
chk("x2 = 0 is a tritangent: S_E(x0,x1,0) = x0^2 (1920 x0^2 - 2265 x1^2)^2", ev(SE, [x0, x1, 0]) == x0^2*(1920*x0^2 - 2265*x1^2)^2);
{ my(good = 1); foreach([[1,1,0],[1,2,0],[3,-1,0],[2,5,0]], P, my(wv = P[1]*(1920*P[1]^2 - 2265*P[2]^2), z = zmap(P, wv)); z = z/content(z); if(z != [0,2,1] && z != [0,-2,-1], good = 0));
  chk("the lift w = x0 (1920 x0^2 - 2265 x1^2) of x2 = 0 is contracted to (0:2:1)", good); }
{ my(z = zmap([1,1,0], 345)); z = z/content(z); if(z[1] < 0, z = -z); chk("(1:1:0; w = 345) maps to (118:-118:-75) with S_3 = 7570880^2", z == [118,-118,-75] && ev(S3x, z) == 7570880^2); }
\\ 4. the new model
grad(Sp, P) = vector(3, m, ev(deriv(Sp, [x0,x1,x2][m]), P));
chk("S_3 is singular at (0:0:1), (0:2:1), (2:-2:-1)", grad(S3x, [0,0,1]) == [0,0,0] && grad(S3x, [0,2,1]) == [0,0,0] && grad(S3x, [2,-2,-1]) == [0,0,0]);
LN = rows(Str(ROOT, "/S3_lines.txt"));
restr(Sp, l) = if(l[3] != 0, ev(Sp, [u, v, -(l[1]*u + l[2]*v)/l[3]]), if(l[2] != 0, ev(Sp, [u, -l[1]*u/l[2], v]), ev(Sp, [0, u, v])));
{ my(nd = [[0,0,1],[0,2,1],[2,-2,-1]], nt = 0, nn = 0, allsq = 1);
  for(i = 1, #LN, my(l = LN[i], Fr = numerator(restr(S3x, l))); if(!issquare(Fr), allsq = 0); if(vecmin(vector(3, k, abs(l*nd[k]~))) == 0, nn++, nt++));
  chk("S_3 is a square on the 21 tritangents and 18 node lines of the lines file", allsq && nt == 21 && nn == 18);
  my(n1 = 0, n2 = 0); for(i = 1, #LN, my(k = sum(j = 1, 3, LN[i]*nd[j]~ == 0)); if(k == 1, n1++); if(k == 2, n2++));
  chk("15 node lines pass through one node and 3 through two", n1 == 15 && n2 == 3); }
chk("S_3(6x-2, 4x+16, 10x+7) = 64 f738", subst(subst(subst(S3, z0, 6*x-2), z1, 4*x+16), z2, 10*x+7) == 64*f738);
chk("(6x-2, 4x+16, 10x+7) lies on the line (66, 31, -52)", 66*(6*x-2) + 31*(4*x+16) - 52*(10*x+7) == 0);
\\ 5. the curve
chk("f738 irreducible with Galois group S6", polisirreducible(f738) && polgalois(f738)[1] == 720);
{ my(d = poldisc(f738), c = d / (2^9*3^6*5^10*7^9*11^8*13^2*29^2*37^3*61^2*137^3)); chk("disc f738 = 2^9 3^6 5^10 7^9 11^8 13^2 29^2 37^3 61^2 137^3 * (28-digit prime)", denominator(c) == 1 && isprime(abs(c)) && #Str(abs(c)) == 28); }
chk("disc f738 is negative", sign(poldisc(f738)) == -1);
c_ft(p) = { my(fa = factormod(f738, p)); if(vecmax(fa[,2]) > 1, return(0)); vecsort(apply(poldegree, fa[,1]~)) };
chk("f738 has irreducible factors of degrees [6], [1,5], [1,2,3] modulo 67, 47, 59", c_ft(67) == [6] && c_ft(47) == [1, 5] && c_ft(59) == [1, 2, 3]);
{ my(n43 = subst(hyperellcharpoly(Mod(f738, 43)), x, 1), n47 = subst(hyperellcharpoly(Mod(f738, 47)), x, 1));
  chk("#J(F_43) = 3022 = 2*1511 and #J(F_47) = 3537 = 3^3*131 are coprime", n43 == 3022 && n47 == 3537 && 2*1511 == 3022 && 27*131 == 3537 && gcd(n43, n47) == 1 && poldisc(f738) % 43 != 0 && poldisc(f738) % 47 != 0); }
chk("Frobenius at 43", hyperellcharpoly(Mod(f738, 43)) == x^4 + 22*x^3 + 204*x^2 + 946*x + 1849);
PT = rows(Str(ROOT, "/points_738.txt"));
{ my(good = 1, cnt = vector(9), X = Set(vector(#PT, i, PT[i][1]/PT[i][3])));
  for(i = 1, #PT, my(a = PT[i][1], b = PT[i][2], c = PT[i][3]); if(gcd(a, c) != 1 || (b/c^3)^2 != subst(f738, x, a/c) || b == 0, good = 0); for(e = 1, 9, if(max(abs(a), abs(c)) <= 10^e, cnt[e] += 2)));
  chk("369 distinct x-values, all on the curve, none with y = 0", good && #X == 369 && #PT == 369);
  chk("points by height 10..10^9: 36 188 390 554 666 712 732 736 738", cnt == [36, 188, 390, 554, 666, 712, 732, 736, 738]);
  chk("largest x-value 185880347/11934737", vecmax(vector(#PT, i, max(abs(PT[i][1]), abs(PT[i][3])))) == 185880347 && setsearch(X, 185880347/11934737));
  my(big = vecsort(select(h -> h > 10^5, vector(#PT, i, max(abs(PT[i][1]), abs(PT[i][3])))))); chk("333 x-values of height <= 10^5 and 36 above, from 1.06e5 to 1.86e8", #big == 36 && #PT - #big == 333 && big[1] == 105599 && big[36] == 185880347);
  chk("-102/109 and 1474359/1473364 are x-values of the list", setsearch(X, -102/109) && setsearch(X, 1474359/1473364)); }
chk("f6 = 743425 is not a square", !issquare(743425));
\\ 6. the closure example
{ my(xs = [0, -3/2, 1/3, -5/3], ys = [2615, 36025/8, 46475/27, 248815/27]);
  chk("P1..P4 lie on the curve", vector(4, i, subst(f738, x, xs[i]) == ys[i]^2) == [1,1,1,1]);
  chk("p interpolates P1..P4", vector(4, i, subst(p, x, xs[i]) == ys[i]) == [1,1,1,1]);
  chk("p^2 - f = (2/81) x (2x+3)(3x-1)(3x+5)(109x+102)(1473364x-1474359)", p^2 - f738 == res);
  chk("y at x = 1474359/1473364", subst(p, x, 1474359/1473364) == -47232542363775629256305/1473364^3); }
\\ 7. the other curves, and the sextics of the search
\\ the row of data/search_curves.txt for a model and a line: the reduced sextic g0..g6 used in the search, as a polynomial
currow(model, l) = { my(R = readstr(Str(ROOT, "/search_curves.txt")), g = 0);
  for(i = 1, #R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = select(t -> #t > 0, strsplit(R[i], " ")));
    if(w[3] == model && [eval(w[4]), eval(w[5]), eval(w[6])] == l, g = Pol(Vecrev(vector(7, k, eval(w[6 + k]))))));
  g }
chk("the search's sextic on the line (66, 31, -52) of t0453 is 64 f738(-x/2): the equation of the abstract with x -> -2x", currow("t0453", [66, 31, -52]) == 64*subst(f738, x, -x/2));
chk("the search's sextic on the line (95, -9, -82) of b4 is f642", currow("b4", [95, -9, -82]) == f642);
\\ the sextic of a row of data/models.txt
modelrow(name) = { my(R = readstr(Str(ROOT, "/models.txt")), S = 0);
  for(i = 1, #R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = strsplit(R[i], " | "));
    if(w[1] == name, S = sx(apply(eval, select(t -> #t > 0, strsplit(w[3], " "))))));
  S }
maxc(S) = vecmax(abs(Vec(apply(t -> polcoeff(polcoeff(polcoeff(S, t[1], x0), t[2], x1), t[3], x2), MON6))));
hl(l) = vecmax(abs(l));
{ my(E = modelrow("E"), R = readstr(Str(ROOT, "/models.txt")), lE = 0);
  for(i = 1, #R, if(#R[i] == 0 || Vecsmall(R[i])[1] == 35, next); my(w = strsplit(R[i], " | ")); if(w[1] == "E", lE = apply(t -> eval(Str("[", t, "]")), strsplit(w[5], ";"))));
  chk("Elkies' sextic: largest coefficient 8697600, lines of height at most 148; its minimal form (row E): 35356 and 37, 52 lines",
    maxc(SE) == 8697600 && vecmax(apply(hl, rows(Str(ROOT, "/elkies163_lines.txt")))) == 148 && maxc(E) == 35356 && #lE == 52 && vecmax(apply(hl, lE)) == 37); }
{ my(S5 = modelrow("b5"), T = substvec(S5, [x0, x1, x2], [x0 - x2, x0 + x2, x1 - x2]), M = [1,0,-1;1,0,1;0,1,-1], l = [21, -256, -50]*M^-1);
  chk("the model b5: largest coefficient 6984; S(z0 - z2, z0 + z2, z1 - z2)/4 is integral with largest coefficient 5134; the determinant is 2; its line (21, -256, -50) is the line (327, -285, -512) of b5",
    maxc(S5) == 6984 && content(T) == 4 && maxc(T/4) == 5134 && abs(matdet(M)) == 2 && l/content(l) == [327, -285, -512]); }
{ my(R = rows(Str(ROOT, "/minimisation.txt")), n = 0); for(i = 2, #R, if(abs(R[i][4]) != 1, n++));
  chk("minimisation.txt: 2380 sextics constructed, 397 of them not minimal as they come out; Elkies' sextic 8697600 -> 35356 with determinant 1/4", #R == 2381 && n == 397 && R[1][2] == 8697600 && R[1][3] == 35356 && abs(R[1][4]) == 1/4); }
chk("the row t0453 of data/models.txt is S_3", modelrow("t0453") == S3x);
chk("the search's sextic on the line (83, 6, 125) of t1272 is f654(-2x)/16", currow("t1272", [83, 6, 125]) == subst(f654, x, -2*x)/16);
chk("the 648-point curve is the even sextic [83877673125, 7793336764011, -1685253906933, 98285087781] (x -> -x symmetry, so a split Jacobian)", f648 == 98285087781*x^6 - 1685253906933*x^4 + 7793336764011*x^2 + 83877673125 && subst(f648, x, -x) == f648 && issquarefree(f648));
chk("C654 irreducible with Galois group S6", polisirreducible(f654) && polgalois(f654)[1] == 720);
chk("C642 irreducible with Galois group S6", polisirreducible(f642) && polgalois(f642)[1] == 720);
print("GP SCRIPT COMPLETE");
quit
''')
script = 'ROOT = "%s";\n' % ROOT + '\n'.join(gp)
assert len(re.findall(r'\bchk\((?!name, b\))', script)) == N_CHECKS, 'N_CHECKS does not match the number of chk calls in the GP script'
tmpf = tempfile.NamedTemporaryFile('w', suffix='.gp', delete=False); tmpf.write(script); tmpf.close()
r = subprocess.run(['gp', '-q', tmpf.name], stdin=subprocess.DEVNULL, capture_output=True, text=True); os.unlink(tmpf.name)
out = r.stdout.splitlines(); print('\n'.join(l for l in out if l != 'GP SCRIPT COMPLETE'))
err = [l for l in r.stderr.splitlines() if 'Warning' not in l]
if err: print('\n'.join(err[:20]))
# The verdict is not gp's: after an error gp skips the statement and carries on, so a run can reach its last line with checks missing.
n_ok = sum(l.startswith('ok    ') for l in out); n_fail = sum(l.startswith('FAIL  ') for l in out); why = []
if r.returncode != 0: why.append('gp exited with status %d' % r.returncode)
if any('***' in l or 'skipping file' in l for l in err + [l for l in out if 'Warning' not in l]): why.append('gp reported an error')
if 'GP SCRIPT COMPLETE' not in out: why.append('the GP script did not reach its end')
if n_fail: why.append('%d checks failed' % n_fail)
if n_ok != N_CHECKS: why.append('%d of %d checks ran and passed' % (n_ok, N_CHECKS))
print('VERIFICATION FAILED: ' + '; '.join(why) if why else 'ALL CHECKS PASSED (%d of %d)' % (n_ok, N_CHECKS))
sys.exit(1 if why else 0)

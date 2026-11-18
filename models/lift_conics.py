#!/usr/bin/env python3
"""lift_conics.py -- Hensel-lift conics found mod p by conicmodp to characteristic 0 and reconstruct them over Q.

usage: lift_conics.py sextic.txt p conics_mod_p.txt [-K steps] [-o out.txt]

For each smooth conic q mod p on which S = c*square, we build the identity S = g^2 + q*h over F_p (g a cubic form
in the 7 standard monomials not divisible by the leading monomial of q, h a quartic), Hensel-lift the 27 unknowns
(5 free coefficients of q, 7 of g, 15 of h) against the 28 coefficient equations of the sextic identity (the 28th
equation is a consistency check at every step), and rationally reconstruct q, g, h from Z/p^K.  Output: one line per
conic "c00 c01 c02 c11 c12 c22 | g (10 rationals) | h (15 rationals)" with q primitive integral, plus the exact
identity S - g^2 - q*h == 0 verified over Q with Fraction arithmetic before printing.
"""
import sys, itertools
from fractions import Fraction
from math import gcd

def readsextic(fn):
    c = [int(x) for x in open(fn).read().split()]
    assert len(c) == 28
    S = {}
    idx = 0
    for i in range(6, -1, -1):
        for j in range(6 - i, -1, -1):
            S[(i, j, 6 - i - j)] = c[idx]; idx += 1
    return S

def monos(d):
    return [(i, j, d - i - j) for i in range(d, -1, -1) for j in range(d - i, -1, -1)]

MON2, MON3, MON4, MON6 = monos(2), monos(3), monos(4), monos(6)
IDX6 = {m: n for n, m in enumerate(MON6)}

def pmul(A, B):
    R = {}
    for ma, ca in A.items():
        if not ca: continue
        for mb, cb in B.items():
            if not cb: continue
            m = (ma[0] + mb[0], ma[1] + mb[1], ma[2] + mb[2])
            R[m] = R.get(m, 0) + ca * cb
    return R

def padd(A, B, s=1):
    R = dict(A)
    for m, c in B.items(): R[m] = R.get(m, 0) + s * c
    return R

def tovec(A):
    return [A.get(m, 0) for m in MON6]

def modinv(a, p): return pow(a % p, -1, p)

# ---- linear algebra mod p ----
def solve_mod(M, b, p):
    """solve M x = b mod p (M: list of rows, may be overdetermined but consistent); returns x and pivot rows"""
    n = len(M); m = len(M[0])
    A = [[M[i][j] % p for j in range(m)] + [b[i] % p] for i in range(n)]
    piv = []; r = 0
    for c in range(m):
        pr = next((i for i in range(r, n) if A[i][c]), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        inv = modinv(A[r][c], p)
        A[r] = [(v * inv) % p for v in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(m + 1)]
        piv.append(c); r += 1
    for i in range(r, n):
        if A[i][m] % p: raise ValueError("inconsistent system")
    x = [0] * m
    for i, c in enumerate(piv): x[c] = A[i][m]
    return x, r

def inv_mod(M, p):
    n = len(M)
    A = [[M[i][j] % p for j in range(n)] + [1 if i == j else 0 for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c]), None)
        if pr is None: raise ValueError("singular")
        A[c], A[pr] = A[pr], A[c]
        inv = modinv(A[c][c], p)
        A[c] = [(v * inv) % p for v in A[c]]
        for i in range(n):
            if i != c and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[c][j]) % p for j in range(2 * n)]
    return [row[n:] for row in A]

def independent_rows(M, p):
    """indices of a maximal set of independent rows of M mod p (greedy)"""
    n = len(M); m = len(M[0])
    basis = []  # reduced rows
    rows = []
    red = []
    for i in range(n):
        v = [x % p for x in M[i]]
        for (pc, brow) in red:
            if v[pc]:
                f = v[pc]
                v = [(v[j] - f * brow[j]) % p for j in range(m)]
        pc = next((j for j in range(m) if v[j]), None)
        if pc is None: continue
        inv = modinv(v[pc], p)
        v = [(x * inv) % p for x in v]
        red.append((pc, v)); rows.append(i)
        if len(rows) == m: break
    return rows

# ---- binary forms mod p: list of coefficients of s^(d-i) t^i ----
def bmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    return r

def compose(S, Q, p):
    """S(Q0,Q1,Q2) mod p, Q_i binary quadratic forms; returns degree-12 form"""
    pw = []
    for i in range(3):
        P = [[1]]
        for e in range(1, 7): P.append(bmul(P[-1], Q[i], p))
        pw.append(P)
    F = [0] * 13
    for (i, j, k), c in S.items():
        if c % p == 0: continue
        m = bmul(bmul(pw[0][i], pw[1][j], p), pw[2][k], p)
        for n, v in enumerate(m): F[n] = (F[n] + c * v) % p
    return F

def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    # Tonelli-Shanks
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

def form_sqrt(F, p):
    """G with G^2 = F (binary forms, deg F = 12, deg G = 6) mod p, or None"""
    # dehomogenise at a shift a with F(1,a) != 0 : f(t) = F(1, t+a)
    d = 12
    for a in range(p):
        f = [0] * 13
        # F(1,t) = sum F[i] t^i ; substitute t -> t + a
        for i, c in enumerate(F):
            if not c: continue
            # (t+a)^i
            bc = 1
            for k in range(i + 1):
                f[k] = (f[k] + c * bc * pow(a, i - k, p)) % p
                bc = bc * (i - k) // (k + 1)
        if f[0]: break
    else:
        return None
    c0 = f[0]; r0 = sqrt_mod(c0, p)
    if r0 is None: return None
    inv = modinv(c0, p)
    H = [(x * inv) % p for x in f]     # H = 1 + ...
    # power series sqrt to degree 6: G = 1 + sum g_k t^k, 2 g_k = H_k - sum_{i=1}^{k-1} g_i g_{k-i}
    g = [1] + [0] * 6
    i2 = modinv(2, p)
    for k in range(1, 7):
        acc = H[k]
        for i in range(1, k): acc -= g[i] * g[k - i]
        g[k] = (acc * i2) % p
    g = [(x * r0) % p for x in g]
    # shift back: G(1,t) = g(t - a)
    G = [0] * 7
    for i, c in enumerate(g):
        if not c: continue
        bc = 1
        for k in range(i + 1):
            G[k] = (G[k] + c * bc * pow(-a, i - k, p)) % p
            bc = bc * (i - k) // (k + 1)
    # verify G^2 = F
    G2 = bmul(G, G, p)
    if [x % p for x in G2] != [x % p for x in F]: return None
    return G

def parametrize_modp(q, p):
    """point and parametrisation Q(s,t) (3 binary quadratic forms) of the conic q (6 coeffs) mod p"""
    c = [x % p for x in q]
    def qf(x): return (c[0]*x[0]*x[0] + c[1]*x[0]*x[1] + c[2]*x[0]*x[2] + c[3]*x[1]*x[1] + c[4]*x[1]*x[2] + c[5]*x[2]*x[2]) % p
    def twoB(x, y): return (qf([x[i] + y[i] for i in range(3)]) - qf(x) - qf(y)) % p
    P = None
    if c[0] == 0: P = [1, 0, 0]
    else:
        for x1, x2 in itertools.chain(((1, k) for k in range(p)), [(0, 1)]):
            b = (c[1]*x1 + c[2]*x2) % p; cc = (c[3]*x1*x1 + c[4]*x1*x2 + c[5]*x2*x2) % p
            disc = (b*b - 4*c[0]*cc) % p
            r = sqrt_mod(disc, p)
            if r is not None:
                P = [((-b + r) * modinv(2 * c[0], p)) % p, x1, x2]; break
    if P is None: return None
    if P[0]: M1, M2 = [0,1,0], [0,0,1]
    elif P[1]: M1, M2 = [1,0,0], [0,0,1]
    else: M1, M2 = [1,0,0], [0,1,0]
    qa, qb, qc = qf(M1), twoB(M1, M2), qf(M2)
    b1, b2 = twoB(P, M1), twoB(P, M2)
    Q = []
    for i in range(3):
        Q.append([(qa*P[i] - b1*M1[i]) % p, (qb*P[i] - b1*M2[i] - b2*M1[i]) % p, (qc*P[i] - b2*M2[i]) % p])
    return P, Q

def ratrecon(a, m):
    """rational reconstruction of a mod m: n/d with |n|,|d| <= sqrt(m/2)"""
    a %= m
    r0, r1 = m, a
    s0, s1 = 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound: raise ValueError("ratrecon failed")
    if s1 < 0: r1, s1 = -r1, -s1
    if gcd(r1, s1) != 1: raise ValueError("ratrecon failed (gcd)")
    return Fraction(r1, s1)

def lift_one(S, qmod, p, K):
    """returns (q_int (6 ints, primitive), g (10 Fractions), h (15 Fractions)) or raises"""
    par = parametrize_modp(qmod, p)
    if par is None: raise ValueError("no point mod p")
    P, Q = par
    F = compose(S, Q, p)
    G = form_sqrt(F, p)
    if G is None: raise ValueError("not a square mod p")
    # leading monomial m0 = first nonzero coefficient of q; standard cubic monomials = those not divisible by m0
    i0 = next(i for i in range(6) if qmod[i] % p)
    m0 = MON2[i0]
    std3 = [m for m in MON3 if not (m[0] >= m0[0] and m[1] >= m0[1] and m[2] >= m0[2])]
    assert len(std3) == 7
    free_q = [i for i in range(6) if i != i0]
    # solve g(Q(s,t)) = G(s,t): for each standard monomial, its composition
    cols = []
    pw = []
    for i in range(3):
        Pw = [[1]]
        for e in range(1, 4): Pw.append(bmul(Pw[-1], Q[i], p))
        pw.append(Pw)
    for (i, j, k) in std3:
        cols.append(bmul(bmul(pw[0][i], pw[1][j], p), pw[2][k], p))
    Mg = [[cols[c][r] for c in range(7)] for r in range(7)]
    gcoef, rk = solve_mod(Mg, G, p)
    if rk < 7: raise ValueError("g system deficient")
    g = {m: gcoef[n] for n, m in enumerate(std3)}
    q = {MON2[i]: qmod[i] % p for i in range(6)}
    # h from q*h = S - g^2 mod p
    rhs = tovec(padd(S, pmul(g, g), -1))
    Mh = [[0] * 15 for _ in range(28)]
    for c, m in enumerate(MON4):
        col = tovec(pmul(q, {m: 1}))
        for r in range(28): Mh[r][c] = col[r]
    hcoef, rk = solve_mod(Mh, rhs, p)
    if rk < 15: raise ValueError("h system deficient")
    h = {m: hcoef[n] for n, m in enumerate(MON4)}
    # unknown vector X = (q free (5), g (7), h (15)); integers (lifted representatives)
    X = [q[MON2[i]] for i in free_q] + [g[m] for m in std3] + [h[m] for m in MON4]
    def unpack(X):
        qq = {MON2[i]: 0 for i in range(6)}; qq[m0] = 1
        for n, i in enumerate(free_q): qq[MON2[i]] = X[n]
        gg = {m: X[5 + n] for n, m in enumerate(std3)}
        hh = {m: X[12 + n] for n, m in enumerate(MON4)}
        return qq, gg, hh
    def Phi(X):
        qq, gg, hh = unpack(X)
        return tovec(padd(padd(pmul(gg, gg), pmul(qq, hh)), S, -1))
    # Jacobian mod p (28 x 27)
    qq, gg, hh = unpack(X)
    J = [[0] * 27 for _ in range(28)]
    col = 0
    for i in free_q:
        v = tovec(pmul({MON2[i]: 1}, hh))
        for r in range(28): J[r][col] = v[r] % p
        col += 1
    for m in std3:
        v = tovec(pmul({m: 2}, gg))
        for r in range(28): J[r][col] = v[r] % p
        col += 1
    for m in MON4:
        v = tovec(pmul({m: 1}, qq))
        for r in range(28): J[r][col] = v[r] % p
        col += 1
    rows = independent_rows(J, p)
    if len(rows) < 27: raise ValueError("Jacobian rank %d < 27" % len(rows))
    A = [[J[r][c] for c in range(27)] for r in rows]
    Ainv = inv_mod(A, p)
    pk = p
    for k in range(1, K):
        R = Phi(X)
        if any(v % pk for v in R): raise ValueError("consistency failure at step %d" % k)
        rhs = [(-(R[r] // pk)) % p for r in rows]
        delta = [sum(Ainv[i][j] * rhs[j] for j in range(27)) % p for i in range(27)]
        X = [X[i] + pk * delta[i] for i in range(27)]
        pk *= p
    R = Phi(X)
    if any(v % pk for v in R): raise ValueError("consistency failure at final step")
    # rational reconstruction
    Xr = [ratrecon(x, pk) for x in X]
    qq, gg, hh = unpack(Xr)
    qq[m0] = Fraction(1)
    # exact verification over Q
    E = padd(padd(pmul(gg, gg), pmul(qq, hh)), {m: Fraction(c) for m, c in S.items()}, -1)
    if any(v != 0 for v in E.values()): raise ValueError("exact identity fails after reconstruction")
    # make q primitive integral: q_int = q * L, then g^2 + q h = S with h -> h / L
    den = 1
    for m in MON2: den = den * qq[m].denominator // gcd(den, qq[m].denominator)
    qi = [int(qq[m] * den) for m in MON2]
    gq = 0
    for v in qi: gq = gcd(gq, abs(v))
    qi = [v // gq for v in qi]
    sgn = 1 if next(v for v in qi if v) > 0 else -1
    qi = [sgn * v for v in qi]
    scale = Fraction(sgn * den, gq)      # q_int = scale * qq
    hh = {m: v / scale for m, v in hh.items()}
    gl = [gg.get(m, Fraction(0)) for m in MON3]
    hl = [hh.get(m, Fraction(0)) for m in MON4]
    return qi, gl, hl

def main():
    args = sys.argv[1:]
    K = 40; out = None
    pos = []
    i = 0
    while i < len(args):
        if args[i] == '-K': K = int(args[i + 1]); i += 2
        elif args[i] == '-o': out = args[i + 1]; i += 2
        else: pos.append(args[i]); i += 1
    S = readsextic(pos[0]); p = int(pos[1]); fn = pos[2]
    fo = open(out, 'w') if out else sys.stdout
    nok = nfail = 0
    for line in open(fn):
        w = line.split()
        if len(w) < 6 or w[0].startswith('#'): continue
        qmod = [int(x) for x in w[:6]]
        try:
            qi, gl, hl = lift_one(S, qmod, p, K)
        except ValueError as e:
            nfail += 1
            print("# FAIL", ' '.join(w[:6]), str(e), file=sys.stderr)
            continue
        nok += 1
        print(' '.join(str(v) for v in qi), '|', ' '.join(str(v) for v in gl), '|', ' '.join(str(v) for v in hl), file=fo)
    print("p=%d: lifted %d, failed %d" % (p, nok, nfail), file=sys.stderr)

if __name__ == '__main__':
    main()

# CRT and rational reconstruction of the echelon basis of the inverse-map solutions
import sys, glob, re, math
from fractions import Fraction
files = sorted(glob.glob('k34_p*.txt'), key=lambda f: int(re.findall(r'p(\d+)', f)[0]))
data = []
for f in files:
    L = open(f).read().split('\n')
    if 'DONE' not in L: continue
    p = int(L[0]); nr = int(L[1]); rows = [[int(x) for x in L[2+i].strip('[] ').split(',')] for i in range(nr)]
    data.append((p, rows))
print(len(data), 'primes; rows', {len(d[1]) for d in data})
nuse = int(sys.argv[1]) if len(sys.argv) > 1 else len(data) - 2
use, chk = data[:nuse], data[nuse:]
M = 1
for p, _ in use: M *= p
def crt(vals):
    x = 0
    for (p, _), v in zip(use, vals):
        Mi = M // p; x += v * Mi * pow(Mi, -1, p)
    return x % M
def ratrec(a, m):
    # rational reconstruction with balanced bounds
    bound = math.isqrt(m // 2)
    r0, r1 = m, a; s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q*r1; s0, s1 = s1, s0 - q*s1
    if s1 == 0 or abs(s1) > bound: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return Fraction(r1, s1)
nr = len(use[0][1]); nc = len(use[0][1][0])
out = []; fail = 0; maxd = 0
for i in range(nr):
    row = []
    for j in range(nc):
        v = crt([d[1][i][j] for d in use]); q = ratrec(v, M)
        if q is None: fail += 1; q = Fraction(0)
        row.append(q); maxd = max(maxd, len(str(abs(q.numerator))) + len(str(q.denominator)))
    out.append(row)
print('modulus digits', len(str(M)), ' failures', fail, ' largest numerator+denominator digits', maxd)
# check against the remaining primes
bad = 0
for p, rows in chk:
    for i in range(nr):
        for j in range(nc):
            q = out[i][j]
            if q.denominator % p == 0 or (q.numerator * pow(q.denominator, -1, p) - rows[i][j]) % p: bad += 1
print('check on', len(chk), 'further primes: mismatches', bad)
if fail == 0 and bad == 0:
    with open('inv_int.txt', 'w') as f:                  # the six solutions as primitive integer vectors, for lll.m
        for row in out:
            den = 1
            for q in row: den = den * q.denominator // math.gcd(den, q.denominator)
            v = [int(q * den) for q in row]; g = 0
            for x in v: g = math.gcd(g, x)
            f.write(' '.join(str(x // g) for x in v) + '\n')
    for i in range(nr):
        den = 1
        for q in out[i]: den = den * q.denominator // math.gcd(den, q.denominator)
        print('row', i, 'common denominator digits', len(str(den)), ' max |numerator| digits', max(len(str(abs(q.numerator * (den // q.denominator)))) for q in out[i]))

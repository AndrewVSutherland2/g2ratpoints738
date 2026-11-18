#!/usr/bin/env python3
# pairs.py LAT OUT [NMAX=3] : the polarizations D' = C_a + C_b of a one-node model B, from lat_<B>.txt (base.gp).
#   C_a, C_b run over the lifts of the conics of B with C_a.C_b = 3; then D'^2 = 2.
#   The roots of D'^perp are the classes v with v^2 = -2, v.D' = 0, and such a class has degree |v.D_B| <= 2 (no class E
#   with E^2 = 0, E.D_B = 1 exists, which excludes degree 3), so they are counted from the table R of base.gp:
#   n2(D') = #{ v in R : v.C_a + v.C_b = 0 }.
#   Output: one row per distinct D' with n2 <= NMAX whose roots are pairwise orthogonal (root system A1^n2):
#     "<B>_<k> n2 | d1..d20 | i,si,j,sj,nodeconics,height;..."   (all witness pairs, best first: fewest conics through
#     the node, then smallest height max(h_i, h_j), then (i, j))
import sys, numpy as np
lat, out = sys.argv[1], sys.argv[2]; NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
B = lat.split('lat_')[-1].split('.txt')[0]
X = []; lab = []; R = []
for line in open(lat):
    w = line.split()
    if w[0] == 'GZ': GZ = np.array(w[1:], dtype=np.int64).reshape(20, 20)
    elif w[0] == 'X': lab.append(tuple(int(t) for t in w[1:5])); X.append([int(t) for t in w[5:25]])
    elif w[0] == 'R': R.append([int(t) for t in w[2:22]])
X = np.array(X, dtype=np.int64); R = np.array(R, dtype=np.int64); n = len(X)
GC = X @ GZ @ X.T
A = R @ GZ @ X.T                       # A[r, a] = v_r . C_a
vals = np.unique(A)
N = np.zeros((n, n), dtype=np.float32)   # N[a, b] = #{ r : A[r, a] + A[r, b] = 0 }
for v in vals:
    if -v in vals:
        N += (A == v).astype(np.float32).T @ (A == -v).astype(np.float32)
N = np.rint(N).astype(np.int32)
ia, ib = np.nonzero(np.triu(GC == 3, 1))
n2 = N[ia, ib]
hist = np.bincount(n2, minlength=12)
RR = R @ GZ @ R.T
dp = {}                                 # D' -> [n2, witnesses]
for a, b, k in zip(ia.tolist(), ib.tolist(), n2.tolist()):
    if k > NMAX: continue
    d = tuple((X[a] + X[b]).tolist())
    e = dp.get(d)
    if e is None:
        rts = np.nonzero(A[:, a] + A[:, b] == 0)[0]
        g = RR[np.ix_(rts, rts)]
        orth = bool(np.all(g == np.diag(np.diag(g))))
        e = dp[d] = [k, orth, []]
    (i, si, ti, hi), (j, sj, tj, hj) = lab[a], lab[b]
    e[2].append((ti + tj, max(hi, hj), i, j, si, sj))
nd = np.zeros(NMAX + 1, dtype=int); nk = np.zeros(NMAX + 1, dtype=int)
with open(out, 'w') as fo:
    for k, (d, (m, orth, wit)) in enumerate(sorted(dp.items(), key=lambda t: (t[1][0], min(t[1][2])))):
        nd[m] += 1
        if not orth: continue
        nk[m] += 1
        wit.sort()
        fo.write('%s_%d %d | %s | %s\n' % (B, k + 1, m, ' '.join(map(str, d)), ';'.join('%d,%d,%d,%d,%d,%d' % (i, si, j, sj, tn, h) for tn, h, i, j, si, sj in wit)))
print('%s: %d conic lifts, %d pairs with C_a.C_b = 3; pairs by number of roots of D\'^perp (0..11+): %s' % (B, n, len(ia), hist[:12].tolist()))
print('%s: distinct D\' with n2 = 1..%d: %s; with pairwise orthogonal roots: %s' % (B, NMAX, nd[1:].tolist(), nk[1:].tolist()))

#!/usr/bin/env python3
# choose_pairs.py OUT B1 [B2 ...] [-skip HASHFILE] : one conic pair for every class of polarizations found on the bases B1, B2, ...
# Reads lat/canon_<B>.txt (canonical forms) and lat/dp_<B>.txt (witness pairs).  A class is a canonical form; its pair is
# the best witness on the first base (in the order given) on which the class occurs: fewest conics through the node,
# then smallest height max(h_i, h_j) of the two conics, then smallest (i, j).
# Rows of OUT: "nodes hash base i j sign_i sign_j nodeconics height npolarizations nlines", sorted by nodes, base order, pair;
# nlines = number of classes of line components (twice the number of lines on which the sextic is a square).
import sys, os
args = sys.argv[1:]; D = os.path.dirname(os.path.abspath(__file__))
out = args[0]; bases = args[1:]
best = {}
for bi, B in enumerate(bases):
    canon = {}
    for line in open('%s/lat/canon_%s.txt' % (D, B)):
        w = line.split(); canon[w[0]] = (int(w[1]), w[5], int(w[2]), int(w[4]))
    for line in open('%s/lat/dp_%s.txt' % (D, B)):
        h, d, wit = line.rstrip('\n').split(' | ')
        ident = h.split()[0]; n2, hsh, idx, nl = canon[ident]
        assert idx == 1, (ident, idx)
        i, si, j, sj, tn, ht = (int(t) for t in wit.split(';')[0].split(','))
        key = (bi, tn, ht, i, j, si, sj)
        e = best.get(hsh)
        if e is None: best[hsh] = [n2, key, 1, nl]
        else:
            assert e[0] == n2 and e[3] == nl
            e[2] += 1
            if key < e[1]: e[1] = key
rows = sorted((n2, key, hsh, cnt, nl) for hsh, (n2, key, cnt, nl) in best.items())
with open(out, 'w') as fo:
    for n2, (bi, tn, ht, i, j, si, sj), hsh, cnt, nl in rows:
        fo.write('%d %s %s %d %d %d %d %d %d %d %d\n' % (n2, hsh, bases[bi], i, j, si, sj, tn, ht, cnt, nl))
import collections
c = collections.Counter((r[0], bases[r[1][0]]) for r in rows)
print('classes by nodes and base of the chosen pair:', ' '.join('%d/%s: %d' % (k[0], k[1], v) for k, v in sorted(c.items())))

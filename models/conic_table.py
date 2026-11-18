#!/usr/bin/env python3
# conic_table.py PARAM OUT : the conic table of a model from the output of verify_conics.gp: one row per conic,
# "c00 c01 c02 c11 c12 c22  m00..m22  0" (q, a rational parametrisation, and a column that is not used), sorted by
# (largest |coefficient| of q, q).  This is the format and the order of data/elkies163_conics.txt.
import sys
rows = []
for line in open(sys.argv[1]):
    w = line.split()
    if len(w) < 15 or w[0].startswith('#'): continue
    q = tuple(int(x) for x in w[:6]); par = tuple(int(x) for x in w[6:15])
    rows.append((q, par))
rows.sort(key=lambda c: (max(abs(v) for v in c[0]), c[0]))
with open(sys.argv[2], 'w') as fo:
    for q, par in rows:
        fo.write(' '.join(map(str, q)) + '  ' + ' '.join(map(str, par)) + '  0\n')
print("%d conics, %d through (0:0:1)" % (len(rows), sum(1 for q, p in rows if q[5] == 0)))

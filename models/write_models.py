#!/usr/bin/env python3
# write_models.py OUT NAME=DIR ... | -l LIST DIR : the list of models in the format of data/models.txt,
#   name | number of nodes | 28 coefficients | singular points | lines on which S is a square
# from the files DIR/x163_<NAME>.txt (coefficients) and DIR/x163_<NAME>_lines.txt (the singular points in its header,
# as "singular points [[a, b, c], ...]", then one line "A B C" per row).
import sys, re
def row(name, d, label=None):
    c = open('%s/x163_%s.txt' % (d, name)).read().split()
    assert len(c) == 28, name
    L = open('%s/x163_%s_lines.txt' % (d, name)).read().split('\n')
    m = re.search(r'singular points \[(\[.*?\])\]', L[0]); assert m, name
    sing = [tuple(int(t) for t in s.split(',')) for s in re.findall(r'\[([^\[\]]*)\]', m.group(1))]
    lines = [tuple(int(t) for t in l.split()) for l in L[1:] if l.strip() and not l.startswith('#')]
    return '%s | %d | %s | %s | %s' % (label or name, len(sing), ' '.join(c), ';'.join('%d,%d,%d' % p for p in sing), ';'.join('%d,%d,%d' % l for l in lines))
args = sys.argv[1:]; out = args[0]; rows = []
i = 1
while i < len(args):
    if args[i] == '-l':
        for l in open(args[i + 1]):
            w = l.split()
            if w: rows.append(row(w[0], args[i + 2]))
        i += 3
    else:
        n, d = args[i].split('='); rows.append(row(n, d)); i += 1
HEAD = """# Double-plane models w^2 = S(z0,z1,z2) of the K3 surface X_163 with one, two or three nodes: one row per model,
#   name | number of nodes | the 28 coefficients of S | the singular points of S = 0 | the lines on which S is a square
# The coefficients are listed for d = 0..6 and j = 0..d: the coefficient of z0^(6-d) z1^(d-j) z2^j.  A point is a,b,c
# for (a : b : c), a line a,b,c is a z0 + b z1 + c z2 = 0.  Every sextic is minimal and reduced (Magma, MinRedTernaryForm).
# E is Elkies' model, b1..b6 the other models with one node, p001.. those with two nodes, t0001.. those with three.
# Written by models/build_models.sh from data/elkies163.txt; models.m proves that each row is a model of X_163 with
# exactly these singular points and lines, and that the rows are pairwise non-isomorphic models.
"""
open(out, 'w').write(HEAD + '\n'.join(rows) + '\n')
print('%d models written to %s' % (len(rows), out))

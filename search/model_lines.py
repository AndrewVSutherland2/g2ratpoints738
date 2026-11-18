#!/usr/bin/env python3
"""model_lines.py OUTDIR [LABEL] : the scores and counts of all the lines kept on one model, from the output directory
of  HIST=1 search/search_model.sh NAME OUTDIR  (written to standard output; data/S3_search.txt is this table for the
model t0453, the model of C_738)."""
import sys, re
d = sys.argv[1]; s = open(d + '/summary.txt').read()
m = re.match(r'model (\S+): (\d) nodes, (\d+) lines on which S is a square; L=(\d+) HP=(\d+) lines with score >= (\d+) first count to (\d+), cuts (\d+)/(\d+)/(\d+)', s)
name, L, HP, SMIN, H1, T1, T2 = m[1], m[4], m[5], m[6], m[7], m[8], m[9]
label = sys.argv[2] if len(sys.argv) > 2 else 'model ' + name
hist = [l.split() for l in open(d + '/score_hist.txt')]
sc = [l.split() for l in open(d + '/scores.txt')]
cnt = lambda f: {tuple(l.split()[1:4]): l.split()[0] for l in open(d + '/' + f)}
n1, n4, n5 = cnt('counts_first.txt'), cnt('counts_1e4.txt'), cnt('counts_1e5.txt')
assert sum(int(n) for v, n in hist if int(v) >= int(SMIN)) == len(sc) == len(n1)
print('''# The search on the %s of models.txt (search/search_model.sh with HIST=1, search/model_lines.py).
# First part: the number of lines of height <= %s, not through a node, with a given score (the number of surface points
# of height <= %s on the line), for scores of 15 and more: rows 'score S N'.
# Second part: the %d lines kept (score >= %s), rows 'line a b c score n1 n4 n5' with the numbers of rational points of
# height <= %s, 10^4, 10^5 on the curve the line cuts out, in the reduced coordinates of the search (- : not counted,
# because the previous count was below %s, respectively %s).''' % (label, L, HP, len(sc), SMIN, H1, T1, T2))
for v, n in hist: print('score', v, n)
for a, b, c, v in sc: print('line', a, b, c, v, n1[(a, b, c)], n4.get((a, b, c), '-'), n5.get((a, b, c), '-'))

#!/usr/bin/env python3
"""sample_table.py OUTDIR DATA : the table of search/sample.sh (written to standard output; a summary goes to
OUTDIR/sample_summary.txt).  For each model of OUTDIR/models.txt the curves closed by the search with the lowered
settings (OUTDIR/NAME) are compared with the rows of DATA/search_curves.txt for that model."""
import sys, os, re, collections
O, D = sys.argv[1], sys.argv[2]
names = open(O + '/models.txt').read().split()
nodes = {l.split(' | ')[0]: int(l.split(' | ')[1]) for l in open(D + '/models.txt') if not l.startswith('#')}
std = collections.defaultdict(dict)
for l in open(D + '/search_curves.txt'):
    if l.startswith('#'): continue
    r = l.split(); std[r[2]][tuple(r[3:6])] = r
rows = []; settings = None
for m in names:
    d = O + '/' + m; s = open(d + '/summary.txt').read()
    st = re.search(r'lines with score >= (\d+) first count to (\d+), cuts (\d+)/(\d+)/(\d+) closure from (\d+)', s).groups()
    settings = settings or st; assert st == settings, m
    score = {tuple(l.split()[:3]): l.split()[3] for l in open(d + '/scores.txt')}
    n1 = {tuple(l.split()[1:4]): l.split()[0] for l in open(d + '/counts_first.txt')}; n4 = {tuple(l.split()[1:4]): l.split()[0] for l in open(d + '/counts_1e4.txt')}
    here = {}
    for l in open(d + '/closed.txt'):
        c = l.split(); key = tuple(c[2:5]); here[key] = c
        if key in std[m]:
            r = std[m][key]; assert r[:2] == c[:2] and r[13] == score[key] and r[15] == n1[key] and r[16] == n4[key], (m, key)
            status = 'found'
        else: status = 'score' if int(score[key]) < 40 else 'count2000' if int(n1[key]) < 250 else 'count1e4' if int(n4[key]) < 330 else None
        assert status, (m, key)
        rows.append([c[0], c[1], m, str(nodes[m])] + list(key) + [score[key], n1[key], n4[key], status])
    assert set(std[m]) <= set(here), m                 # nothing that the standard search finds is lost with the lower cuts
rows.sort(key=lambda r: (-int(r[0]), -int(r[1]), r[2], int(r[4]), int(r[5]), int(r[6])))
k = {i: sum(1 for m in names if nodes[m] == i) for i in (1, 2, 3)}
print('''# A sample of %d models of models.txt (all %d with one node, %d with two and %d with three, at regular intervals in the file),
# searched again with the score bound lowered from 40 to %s and the thresholds of the first two counts from 250 and 330
# to %s and %s (search/sample.sh).  Every curve closed (at least %s points of height <= 10^5):
# points after closure, points of height <= 10^5, model, nodes, line a b c, score, points of height <= %s and <= 10^4, and
# 'found' if the search with the standard settings has the curve, else the cut that removes it there: 'score' (below 40),
# 'count2000' (fewer than 250 points of height <= 2000) or 'count1e4' (fewer than 330 of height <= 10^4).
# models: %s''' % (len(names), k[1], k[2], k[3], settings[0], settings[2], settings[3], settings[5], settings[1], ' '.join(names)))
for r in rows: print(' '.join(r))
lost = [r for r in rows if r[10] != 'found']
with open(O + '/sample_summary.txt', 'w') as f:
    f.write('%d models (%d, %d, %d with one, two, three nodes): %d curves closed with the lowered settings, %d of them by the standard search\n' % (len(names), k[1], k[2], k[3], len(rows), len(rows) - len(lost)))
    for i in (1, 2, 3):
        a = [r for r in rows if r[3] == str(i)]; b = [r for r in a if r[10] != 'found']
        f.write('  %d node%s: %d closed, %d not by the standard search (%s); of these with at least 500 / 550 / 600 points: %d / %d / %d; most points %s\n' % (i, '' if i == 1 else 's', len(a), len(b),
                ', '.join('%s: %d' % (c, sum(1 for r in b if r[10] == c)) for c in ('score', 'count2000', 'count1e4')), sum(int(r[0]) >= 500 for r in b), sum(int(r[0]) >= 550 for r in b), sum(int(r[0]) >= 600 for r in b), max([int(r[0]) for r in b], default='-')))

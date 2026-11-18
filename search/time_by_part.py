#!/usr/bin/env python3
# time_by_part.py WORK : where the time of a run of run_all.sh went, by part of the computation.  (When the run searched
# a given list of models there is no part for the models.)
#   models: wall time and processor time of each stage of models/build_models.sh (WORK/times_cpu.txt)
#   search: the jobs of the search overlap (JOBS models at a time), so the wall time of the search is divided among its
#           parts in proportion to the wall seconds the jobs spent in each of them; processor time is added up directly
#           (user + system time of the programs of each part, from the summaries)
import sys, glob, re, os
W = sys.argv[1]
t = {}
for l in open(W + '/all.log'):
    w = l.split()
    if w and w[0] in ('START', 'MODELS', 'DONE'): t[w[0]] = float(w[2])
built = os.path.exists(W + '/times_cpu.txt')
rows = [l.split() for l in open(W + '/times_cpu.txt')] if built else []
st = {}
for (n0, w0, c0), (n1, w1, c1) in zip(rows, rows[1:]): st[n1] = (float(w1) - float(w0), float(c1) - float(c0))
names = ['points', 'scores', 'reduction', 'counts', 'closure', 'total']
wall = dict.fromkeys(names, 0.0); cpu = dict.fromkeys(names, 0.0); cw = [0.0] * 3; cc = [0.0] * 3; n = 0
for p in glob.glob(W + '/search/*/summary.txt'):
    s = open(p).read(); n += 1
    r = re.search(r'\nseconds: points ([0-9.]+), scores ([0-9.]+), reduction ([0-9.]+), counts ([0-9.]+), closure ([0-9.]+), total ([0-9.]+)', s)
    q = re.search(r'\nprocessor seconds: points ([0-9.]+), scores ([0-9.]+), reduction ([0-9.]+), counts ([0-9.]+), closure ([0-9.]+), total ([0-9.]+)', s)
    a = re.search(r'\nseconds in the counts: first ([0-9.]+), to 10\^4 ([0-9.]+), to 10\^5 ([0-9.]+)', s)
    b = re.search(r'\nprocessor seconds in the counts: first ([0-9.]+), to 10\^4 ([0-9.]+), to 10\^5 ([0-9.]+)', s)
    for i, k in enumerate(names): wall[k] += float(r[i + 1]); cpu[k] += float(q[i + 1])
    for i in range(3): cw[i] += float(a[i + 1]); cc[i] += float(b[i + 1])
SW = t['DONE'] - t['MODELS']; MW = t['MODELS'] - t['START']
def row(name, w, c, ind=0): print('%-58s %8.1f min %8.2f h' % ('  ' * ind + name, w / 60, c / 3600))
print('%-58s %12s %10s' % ('part', 'wall time', 'processor'))
mc = sum(c for w, c in st.values())
if built:
    row('building the models (from Elkies\' sextic to %d sextics)' % sum(1 for l in open(W + '/models.txt') if not l.startswith('#')), MW, mc)
    row('lines of Elkies\' sextic', *st['lines'], 1)
    row('Elkies\' model: conics, lattice, conic pairs, classes', *st['base_E'], 1)
    row('the other one-node models: built, then the same', st['bases_built'][0] + st['bases_classified'][0], st['bases_built'][1] + st['bases_classified'][1], 1)
    row('the sextics with two and three nodes', st['selected'][0] + st['models_built'][0] + st['done'][0], st['selected'][1] + st['models_built'][1] + st['done'][1], 1)
f = SW / wall['total']
row('enumerating points on the surfaces', f * wall['points'], cpu['points'])
row('ranking lines (points on each line, selection)', f * wall['scores'], cpu['scores'])
chk = wall['reduction'] + wall['counts'] + wall['closure']; chc = cpu['reduction'] + cpu['counts'] + cpu['closure']
row('checking the curves of the lines kept', f * chk, chc)
row('restriction to the line and reduction', f * wall['reduction'], cpu['reduction'], 1)
row('point count to height 2000', f * cw[0], cc[0], 1)
row('point count to height 10^4', f * cw[1], cc[1], 1)
row('point count to height 10^5', f * cw[2], cc[2], 1)
row('closure', f * wall['closure'], cpu['closure'], 1)
row('total', t['DONE'] - t['START'], mc + cpu['total'])
print('search: %d models, %.0f s wall; the jobs account for %.0f s x threads' % (n, SW, wall['total']))

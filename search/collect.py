#!/usr/bin/env python3
"""collect.py WORK OUT [-compare DATA] [-classes FILE] : the results of a run of run_all.sh as two tables.

  OUT/search_models.txt   one row per model: what the search did on it
  OUT/search_curves.txt   one row per curve closed; without the last two columns of data/search_curves.txt (Galois group
                          and isomorphism class) unless -classes FILE gives them, one row per curve, as written by
                          magma -b file:=OUT/search_curves.txt classes:=FILE curves.m

Reads WORK/models.txt and the directories WORK/search/NAME written by search/search_model.sh.  With -compare DATA the two
tables are compared with DATA/search_models.txt and DATA/search_curves.txt (times and the two Magma columns excepted),
and the files models.txt, minimisation.txt, elkies163_lines.txt and elkies163_conics.txt of WORK with those of DATA.  When
the run searched a given list of models (it has not built them: WORK has no minimisation.txt), the models must be rows of
DATA/models.txt and the comparison is with the stored rows of these models.  The exit status is 1 if anything differs."""
import sys, os, re

MODELS_HEAD = """# The search of all lines of height <= {L} on the models of models.txt, one row per model:
# name, nodes, forced points, surface points of height <= {HP} off the lines on which S is a square, lines of height <= {L}
# with at least 15 of them, lines kept (score >= {SMIN}, not through a node), the highest score, curves counted to {H1}
# (lines kept that cut out a curve of genus 2), curves with >= {T1} points of height <= {H1} (counted to 10^4), with >= {T2}
# of height <= 10^4 (counted to 10^5), with >= {T3} of height <= 10^5, curves closed (>= {CMIN}), the largest count at 10^5
# and the largest total after closure (- if none); then processor seconds for the surface points, the scores, the
# reduction, the counts and the closure, and wall seconds of the model on {P} threads.
"""
CURVES_HEAD = """# The {n} curves closed by the search (at least {CMIN} points of height <= 10^5), by number of points:
# points after closure, points of height <= 10^5, model, line a b c (a z0 + b z1 + c z2 = 0 on the model of models.txt),
# the reduced sextic g0 .. g6 of the search (y^2 = g0 + g1 x + .. + g6 x^6), the score of the line, its position among
# the lines of its model by score (1 + the number of lines with a higher score), the points of height <= {H1} and <= 10^4
# on the reduced sextic"""
CURVES_HEAD2 = """, the Galois group of the sextic (S6, or its number among the transitive groups of degree 6), and a class
# number: curves isomorphic over an algebraic closure of Q share it.  curves.m checks the last two columns and that each
# sextic is the section of the model by the line.
"""

def read_models(path):
    return [[t.strip() for t in l.split(' | ')] for l in open(path) if not l.startswith('#')]

def collect(work):
    models = read_models(work + '/models.txt'); mrows = []; crows = []; settings = None
    for name, k, coef, nodes, lines in models:
        d = '%s/search/%s' % (work, name); s = open(d + '/summary.txt').read()
        m = re.match(r'model (\S+): (\d) nodes, (\d+) lines on which S is a square; L=(\d+) HP=(\d+) lines with score >= (\d+) first count to (\d+), cuts (\d+)/(\d+)/(\d+) closure from (\d+); (\d+) threads\n', s)
        assert m and m[1] == name and m[2] == k and int(m[3]) == len(lines.split(';')), name
        st = tuple(m[i] for i in range(4, 13))
        if settings is None: settings = st
        assert st == settings, (name, st, settings)
        L, HP, SMIN, H1, T1, T2, T3, CMIN, P = st
        g = lambda pat: re.search(pat, s)
        npts = g(r'surface points of height <= \d+ off those lines: (\d+)')[1]; n15 = g(r'lines of height <= \d+ with at least 15 of them: (\d+)')[1]
        kept = int(g(r'lines kept: (\d+)')[1]); c1 = g(r'curves counted to \d+: (\d+); with at least \d+ points: (\d+)\n')
        c4 = g(r'curves counted to 10\^4: (\d+); with at least \d+ points: (\d+)\n'); c5 = g(r'curves counted to 10\^5: (\d+); with at least \d+ points: (\d+)\n')
        ncl = int(g(r'curves closed: (\d+)')[1]); cpu = g(r'\nprocessor seconds: points ([0-9.]+), scores ([0-9.]+), reduction ([0-9.]+), counts ([0-9.]+), closure ([0-9.]+), total')
        wall = g(r'\nseconds: .*total ([0-9.]+)')[1]
        sc = [l.split() for l in open(d + '/scores.txt')]; assert len(sc) == kept, name
        scof = {tuple(l[:3]): int(l[3]) for l in sc}; allsc = sorted(scof.values())
        n1 = {tuple(l.split()[1:4]): l.split()[0] for l in open(d + '/counts_first.txt')}; n4 = {tuple(l.split()[1:4]): l.split()[0] for l in open(d + '/counts_1e4.txt')}
        cand = [l.split() for l in open(d + '/candidates.txt')]; closed = [l.split() for l in open(d + '/closed.txt')]
        assert len(n1) == int(c1[1]) and len(n4) == int(c4[1]) == int(c1[2]) and int(c5[1]) == int(c4[2]) and len(cand) == int(c5[2]) and len(closed) == ncl, name
        assert ncl == sum(1 for c in cand if int(c[0]) >= int(CMIN)), name
        best5 = max((int(c[0]) for c in cand), default=0); best = max((int(c[0]) for c in closed), default=0)
        mrows.append([name, k, str(2 * len(lines.split(';'))), npts, n15, str(kept), str(allsc[-1]) if allsc else '-', c1[1], c1[2], c4[2], c5[2], str(ncl), str(best5 or '-'), str(best or '-')]
                     + ['%.0f' % float(cpu[i]) for i in range(1, 6)] + ['%.0f' % float(wall)])
        import bisect
        for c in closed:
            assert len(c) == 12, (name, c)
            key = tuple(c[2:5]); rank = 1 + len(allsc) - bisect.bisect_right(allsc, scof[key])
            crows.append(c[:2] + [name] + c[2:] + [str(scof[key]), str(rank), n1[key], n4[key]])
    crows.sort(key=lambda r: (-int(r[0]), -int(r[1]), r[2], int(r[3]), int(r[4]), int(r[5])))
    keys = dict(zip('L HP SMIN H1 T1 T2 T3 CMIN P'.split(), settings))
    return models, mrows, crows, keys

def main():
    a = sys.argv[1:]; work, out = a[0], a[1]; data = a[a.index('-compare') + 1] if '-compare' in a else None
    models, mrows, crows, keys = collect(work); os.makedirs(out, exist_ok=True)
    extra = [l.split() for l in open(a[a.index('-classes') + 1]) if l.strip()] if '-classes' in a else None
    assert extra is None or (len(extra) == len(crows) and all(len(e) == 2 for e in extra))
    with open(out + '/search_models.txt', 'w') as f:
        f.write(MODELS_HEAD.format(**keys)); f.write(''.join(' '.join(r) + '\n' for r in mrows))
    with open(out + '/search_curves.txt', 'w') as f:
        f.write(CURVES_HEAD.format(n=len(crows), **keys) + (CURVES_HEAD2 if extra else '.\n'))
        f.write(''.join(' '.join(r + (extra[i] if extra else [])) + '\n' for i, r in enumerate(crows)))
    tot = lambda i: sum(int(r[i]) for r in mrows)
    print('%d models: %d lines kept, %d curves counted to 10^4, %d to 10^5, %d with at least %s points of height <= 10^5, %d closed; most points: %s'
          % (len(mrows), tot(5), tot(8), tot(9), tot(10), keys['T3'], len(crows), crows[0][0] if crows else '-'))
    if data is None: return 0
    bad = 0
    need = ['models.txt', 'minimisation.txt', 'elkies163_lines.txt', 'elkies163_conics.txt', 'search_models.txt', 'search_curves.txt']
    missing = [f for f in need if not os.path.exists(data + '/' + f)]
    if missing: print('no stored results to compare with: ' + ', '.join(data + '/' + f for f in missing)); return 1
    built = os.path.exists(work + '/minimisation.txt')       # a run on a given list of models has only models.txt
    ref = [l.split() for l in open(data + '/search_models.txt') if not l.startswith('#')]
    refc = [l.split()[:17] for l in open(data + '/search_curves.txt') if not l.startswith('#')]
    if built:
        for f in need[:4]:
            if open(work + '/' + f).read() != open(data + '/' + f).read(): print(f + ' differs from ' + data + '/' + f); bad = 1
    else:                                                    # the given models must be models of the stored list
        stored = set(l for l in open(data + '/models.txt') if not l.startswith('#'))
        new = [l.split()[0] for l in open(work + '/models.txt') if not l.startswith('#') and l not in stored]
        if new: print('models.txt: %d models are not rows of %s/models.txt: %s' % (len(new), data, ' '.join(new[:10]))); bad = 1
        here = set(r[0] for r in mrows); ref = [r for r in ref if r[0] in here]; refc = [r for r in refc if r[2] in here]
    if sorted(r[:14] for r in ref) != sorted(r[:14] for r in mrows):
        rd = dict((r[0], r[:14]) for r in ref); d = [q[0] for q in mrows if rd.get(q[0]) != q[:14]]
        print('search_models.txt: %d rows here, %d stored; differing: %s' % (len(mrows), len(ref), ' '.join(d[:10]))); bad = 1
    if refc != crows:
        s1 = set(map(tuple, refc)); s2 = set(map(tuple, crows)); print('search_curves.txt: %d curves here, %d stored; %d only here, %d only stored' % (len(crows), len(refc), len(s2 - s1), len(s1 - s2))); bad = 1
    if bad: print('THE RUN DIFFERS FROM THE STORED RESULTS')
    elif built: print('the models, the numbers of the search and the curves closed are those of %s' % data)
    else: print('the numbers of the search and the curves closed on the %d models searched are those of %s (the models were given, not built)' % (len(mrows), data))
    return bad

if __name__ == '__main__': sys.exit(main())

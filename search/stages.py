#!/usr/bin/env python3
# stages.py DIR THREADS : where the time of a run went, by stratum and stage (thread-seconds = wall seconds of the stage x threads of the job)
import sys, glob, os, re, collections
D = sys.argv[1]; P = int(sys.argv[2]); MODELS = sys.argv[3]
nd = {l.split(" | ")[0]: int(l.split(" | ")[1]) for l in open(MODELS) if not l.startswith("#")}
names = ['points', 'scores', 'reduction', 'counts', 'closure', 'total']
t = collections.defaultdict(lambda: [0.0] * 6); n = collections.Counter(); kept = collections.Counter(); c4 = collections.Counter(); c5 = collections.Counter(); cl = collections.Counter(); cand = collections.Counter()
for p in glob.glob(D + '/*/summary.txt'):
    s = open(p).read(); k = nd[os.path.basename(os.path.dirname(p))]
    r = re.search(r'seconds: points ([0-9.]+), scores ([0-9.]+), reduction ([0-9.]+), counts ([0-9.]+), closure ([0-9.]+), total ([0-9.]+)', s)
    for i in range(6): t[k][i] += float(r[i + 1])
    n[k] += 1; kept[k] += int(re.search(r'lines kept: (\d+)', s)[1]); c4[k] += int(re.search(r'counted to 10\^4: (\d+)', s)[1]); c5[k] += int(re.search(r'counted to 10\^5: (\d+)', s)[1])
    cand[k] += int(re.search(r'with at least \d+ points: (\d+)\ncurves closed', s)[1]); cl[k] += int(re.search(r'curves closed: (\d+)', s)[1])
tot = sum(t[k][5] for k in t) * P
print('nodes models  lines kept  to 1e4   to 1e5  >=400  closed | thread-hours: ' + '  '.join(f'{x:>9}' for x in names))
for k in sorted(n):
    print(f'{k:5d} {n[k]:6d} {kept[k]:11d} {c4[k]:7d} {c5[k]:8d} {cand[k]:6d} {cl[k]:7d} |               ' + '  '.join(f'{t[k][i] * P / 3600:9.2f}' for i in range(6)))
print(f'  all {sum(n.values()):6d} {sum(kept.values()):11d} {sum(c4.values()):7d} {sum(c5.values()):8d} {sum(cand.values()):6d} {sum(cl.values()):7d} |               ' + '  '.join(f'{sum(t[k][i] for k in t) * P / 3600:9.2f}' for i in range(6)))
print('share of the total:                                                    ' + '  '.join(f'{100 * sum(t[k][i] for k in t) * P / tot:8.1f}%' for i in range(6)))

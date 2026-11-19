#!/usr/bin/env python3
"""check_census.py -- the statements of the note about the lattice census, checked against the census file of the project in
which the note was written (tools/k3lattice/x163pol/classes_ranked.txt there: one row per class of D^perp met by the random
walks; the key column holds the numbers of pairs of roots, of norm-4 vectors, and of vectors of norm 1/2, 5/2 and 9/2 in the
nontrivial coset).

  * the numbers of classes in Table 1 (root types A1, A1^2, A2, A1^3, A1+A2);
  * no recorded class has a vector of norm 1/2 in the coset, i.e. a class E with E^2 = 0, E.D = 1: every recorded class is a
    polarization (a double plane), none is of the excluded kind D^perp = <r> + M;
  * the visit counts quoted as evidence for completeness, their mean and variance, and the number of walks that found each
    class with one or two nodes (walk_stats.txt next to the census file gives the number of steps and these counts);
  * the figures quoted for the computation: from data/models.txt (the 2381 sextics), data/minimisation.txt,
    data/search_models.txt (one row per sextic), data/search_curves.txt (the curves closed, with the score of their line and
    its rank), data/S3_search.txt (the lines kept on the model of C_738), data/search_sample.txt (the sample searched with
    lower cuts) and data/search_summary.txt (times of the run); the number of lines of height at most 400 is recomputed;
  * the logs of the Magma certificates models.m, linesection.m and curves.m: every row of the data files appears in them with
    the stated verdict.

The census file and the walks are not in the repository of the computation, https://github.com/AndrewVSutherland2/g2ratpoints738,
which has all the other files read here.  There the 14 checks that read them (the first three items above, and the comparison of
data/census.txt with the census) are left out, the other 29 are run, and the last line says so; check_census.log is the output
of all 43.

usage: python3 check_census.py"""
import os, sys, re, collections, math, gzip
HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "..", "..", "tools", "k3lattice", "x163pol", "classes_ranked.txt")
TABLE = {"A1": 7, "A1^2": 186, "A2": 8, "A1^3": 2192, "A1+A2": 550}     # Table 1 of the note
N_RECORDED = 53361
CENSUS = os.path.exists(PATH)     # the census files of the project: where they are absent, the N_CENSUS checks that read them are left out
N_CENSUS = 14
ok = True; ran = 0
def chk(name, cond):
    global ok, ran; ok &= bool(cond); ran += 1; print("%-4s %s" % ("ok" if cond else "FAIL", name))
if CENSUS:
    rows = [ln.split("|") for ln in open(PATH) if ln.strip() and not ln.startswith("#")]
    count = {}; excluded = 0; rootfree = 0; visits = {}; v16 = []
    for r in rows:
        assert len(r) >= 12, r
        typ = r[7].strip(); key = [int(t) for t in r[10].strip().split(",")]
        assert len(key) == 5, key
        count[typ] = count.get(typ, 0) + 1
        excluded += key[2] != 0; rootfree += key[0] == 0
        visits.setdefault(typ, []).append(int(r[9]))
        if typ == "A1^3" and r[8].strip() == "16": v16.append(int(r[9]))
    chk("%d classes recorded" % N_RECORDED, len(rows) == N_RECORDED)
    for typ, n in TABLE.items(): chk("%d classes of root type %s" % (n, typ), count.get(typ, 0) == n)
    chk("no class without roots (no smooth sextic model was met)", rootfree == 0)
    chk("no recorded class has a class E with E^2 = 0, E.D = 1 (0 of %d)" % len(rows), excluded == 0)
    stats = open(os.path.join(os.path.dirname(PATH), "logs", "walk6", "walk_stats.txt")).readline()
    m = re.match(r"workers (\d+), total steps (\d+), mass ([0-9.]+), classes (\d+)", stats)
    steps = int(m.group(2)); mass = float(m.group(3))
    chk("the walks made 9.9 million steps; a class with |Aut| = a is expected 606/a times", m and 9.9e6 <= steps < 9.95e6 and round(steps / mass) == 606 and int(m.group(4)) == N_RECORDED)
    chk("every class of Table 1 was visited at least 4 times", min(min(visits[t]) for t in TABLE) == 4)
    m16 = sum(v16) / len(v16); var16 = sum((v - m16) ** 2 for v in v16) / (len(v16) - 1)
    chk("every three-node class was visited at least 10 times; the 1141 with |Aut| = 16 computed, 37.6 times on average (37.9 expected), with variance 42.5",
        min(visits["A1^3"]) == 10 and len(v16) == 1141 and "%.1f" % m16 == "37.6" and "%.1f" % (steps / mass / 16) == "37.9" and "%.1f" % var16 == "42.5")
    walks = [ln for ln in open(os.path.join(os.path.dirname(PATH), "logs", "walk6", "walk_stats.txt")) if " workers=" in ln]
    found = collections.Counter(re.search(r"workers=(\d+)/(\d+)", ln).groups() for ln in walks)
    chk("seven walks; of the 193 classes with one or two nodes, 185 were found by all seven and 8 by six", m.group(1) == "7" and len(walks) == 193 and found == {("7", "7"): 185, ("6", "7"): 8}
        and count["A1"] + count["A1^2"] == 193)
    chk("heuristic miss probabilities exp(-606/a): about 1e-4, 1e-2, 1e-1 for a = 64, 128, 256", [round(math.log10(math.exp(-steps / mass / a))) for a in (64, 128, 256)] == [-4, -2, -1])
D = os.path.join(HERE, "data")
def rd(name, sep=None): return [ln.split(sep) for ln in open(os.path.join(D, name)) if ln.strip() and not ln.startswith("#")]
def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2
        i = j + 1
    return r
def corr(a, b):
    ma = sum(a) / len(a); mb = sum(b) / len(b)
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
def med(v): v = sorted(v); return (v[(len(v) - 1) // 2] + v[len(v) // 2]) / 2
tex = open(os.path.join(HERE, "curve738.tex")).read()
def intex(*frags): return all(f in tex for f in frags)          # the note says so, in these words
# ---------------------------------------------------------------- the models
MOD = [[t.strip() for t in r] for r in rd("models.txt", "|")]
name = [m[0] for m in MOD]; nodes = {m[0]: int(m[1]) for m in MOD}; coef = {m[0]: [int(t) for t in m[2].split()] for m in MOD}
nlin = {m[0]: len(m[4].split(";")) for m in MOD}; byk = {k: [m[0] for m in MOD if int(m[1]) == k] for k in (1, 2, 3)}
chk("2381 sextics, 7 + 186 + 2188 with one, two, three nodes (E, b1..b6, p001.., t0001..); forced points 104-108, 86-96, 72-86",
    len(MOD) == len(set(name)) == 2381 and [len(byk[k]) for k in (1, 2, 3)] == [7, 186, 2188] and all(len(c) == 28 for c in coef.values())
    and byk[1] == ["E", "b1", "b2", "b3", "b4", "b5", "b6"] and all(m[0] == "p" for m in byk[2]) and all(m[0] == "t" for m in byk[3])
    and [(min(2 * nlin[m] for m in byk[k]), max(2 * nlin[m] for m in byk[k])) for k in (1, 2, 3)] == [(104, 108), (86, 96), (72, 86)])
s3 = [int(t) for t in open(os.path.join(D, "S3_coefficients.txt")).read().split()]
s3l = sorted(tuple(int(t) for t in r) for r in rd("S3_lines.txt")); rowl = sorted(tuple(int(t) for t in l.split(",")) for l in [m for m in MOD if m[0] == "t0453"][0][4].split(";"))
chk("the row t0453 is the sextic S_3 of the note with its 39 lines; its largest coefficient is 6326, and the median of the largest coefficient over the sextics is 6073",
    coef["t0453"] == s3 and s3l == rowl and len(rowl) == 39 and max(abs(c) for c in s3) == 6326 and med([max(abs(c) for c in coef[m]) for m in name]) == 6073 and intex("which is $6073$"))
MIN = rd("minimisation.txt")
chk("minimisation.txt has one row per sextic, in the order of models.txt; its last number of coefficients is the largest coefficient of the sextic listed",
    [r[0] for r in MIN] == name and all(int(r[2]) == max(abs(c) for c in coef[r[0]]) for r in MIN))
# ---------------------------------------------------------------- the search, model by model
SM = rd("search_models.txt"); SC = rd("search_curves.txt"); sm = {r[0]: r for r in SM}
I = lambda r, c: int(r[c])
tot = lambda c, k=0: sum(I(r, c) for r in SM if k in (0, int(r[1])))
chk("search_models.txt has one row per sextic, with its number of nodes and forced points", [r[0] for r in SM] == name and all(int(r[1]) == nodes[r[0]] and int(r[2]) == 2 * nlin[r[0]] for r in SM))
Lh = 400; mu = [1] * (Lh + 1); isp = [True] * (Lh + 1)
for p in range(2, Lh + 1):
    if isp[p]:
        for m in range(p, Lh + 1, p):
            mu[m] = -mu[m]
            if m > p: isp[m] = False
        for m in range(p * p, Lh + 1, p * p): mu[m] = 0
nlines = sum(mu[d] * ((2 * (Lh // d) + 1) ** 3 - 1) for d in range(1, Lh + 1)) // 2
usl = 1e6 * (tot(14) + tot(15)) / (2381 * nlines)
ptm = ["%.1f" % (tot(3, k) / len(byk[k]) / 1e6) for k in (1, 2, 3)]; kpm = ["%d" % round(tot(5, k) / len(byk[k]), -2) for k in (1, 2, 3)]
chk("2.1 x 10^8 lines of height <= 400 per model (computed: %d), 5.1 x 10^11 in all; %s, %s, %s million surface points per model on average; scoring costs %.2f microseconds of processor time per line"
    % ((nlines,) + tuple(ptm) + (usl,)), nlines == 213603553 and "%.1e" % (2381 * nlines) == "5.1e+11"
    and intex("(on average $%s$, $%s$ and $%s$ million" % tuple(ptm), "take $%.2f$ microseconds of processor" % usl, "between half a million and a little over two million such points") and 0.05 < usl < 0.2
    and 500000 <= min(I(r, 3) for r in SM) and max(I(r, 3) for r in SM) < 2200000)
chk("the lines kept have score >= 40: on average %s / %s / %s on a model with one / two / three nodes; 6.6 x 10^6 lines kept in all" % tuple(kpm),
    tot(5) == 6605608 and intex("on average $%s$ lines on a model with one node, $%s$ with two nodes and $%s$ with three" % tuple(kpm)))
closed = collections.Counter(r[2] for r in SC)
chk("119,276 curves with >= 250 points of height <= 2000; 17,338 with >= 330 of height <= 10^4; 6,118 with >= 400 of height <= 10^5; 1,061 closed (the rows of search_curves.txt, model by model), 1,035 up to isomorphism",
    tot(8) == 119276 and tot(9) == 17338 and tot(10) == 6118 and tot(11) == len(SC) == 1061 and all(I(r, 11) == closed.get(r[0], 0) for r in SM)
    and len({r[18] for r in SC}) == 1035 and all(len(r) == 19 and I(r, 1) >= 440 and I(r, 0) >= I(r, 1) and I(r, 13) >= 40 and I(r, 15) >= 250 and I(r, 16) >= 330 for r in SC)
    and all((r[13] == "-") == (I(r, 11) == 0) and (r[13] == "-" or I(r, 13) == max(I(c, 0) for c in SC if c[2] == r[0])) for r in SM))
best = {}
for r in SC: best[r[18]] = (max(best.get(r[18], (0, 0))[0], I(r, 0)), nodes[r[2]])
st = {k: [sum(1 for t, kk in best.values() if kk == k and t >= b) for b in (500, 550, 600)] + [max(t for t, kk in best.values() if kk == k)] for k in (1, 2, 3)}
chk("by number of nodes: 72, 12, 1.6 curves with >= 400 points of height <= 10^5 per model; different curves with >= 500 / 550 / 600 points: 56 / 11 / 2, 202 / 28 / 5, 271 / 47 / 9; most points 642, 638, 738",
    ["%.2g" % (tot(10, k) / len(byk[k])) for k in (1, 2, 3)] == ["72", "12", "1.6"] and st == {1: [56, 11, 2, 642], 2: [202, 28, 5, 638], 3: [271, 47, 9, 738]})
top = [r for r in SC if I(r, 0) >= 600]
chk("sixteen curves with at least 600 points, nine of them on three-node models; three with at least 642: 738 on t0453, 654 on t1272, 642 on b4",
    len(top) == 16 and sum(nodes[r[2]] == 3 for r in top) == 9 and [(r[0], r[2], r[3], r[4], r[5]) for r in SC if I(r, 0) >= 642] == [("738", "t0453", "66", "31", "-52"), ("654", "t1272", "83", "6", "125"), ("642", "b4", "95", "-9", "-82")])
chk("the table of the sixteen curves in the note is the top of search_curves.txt (points, points of height <= 10^5, nodes, forced points, model, line, score, rank)",
    all(("$%s%s$ & $%s$ & $%d$ & $%d$ & \\texttt{%s} & $(%s,\\,%s,\\,%s)$ & $%s$ & $%s$\\\\" % (r[0], "^{*}" if r[17] != "S6" else "", r[1], nodes[r[2]], 2 * nlin[r[2]], r[2], r[3], r[4], r[5], r[13], r[14])) in tex for r in top))
chk("the 654-point curve: 570 points of height <= 10^5 and 42 further x-values, on a model with 84 forced points where no other line has 440 points of height <= 10^5; the 642-point curve: on a model with 106 forced points",
    [r[:2] for r in SC if r[2] == "t1272"] == [["654", "570"]] and (654 - 570) // 2 == 42 and sm["t1272"][2] == "84" and sm["b4"][2] == "106" and sm["t0453"][2] == "78" and sm["t0453"][10:12] == ["1", "1"])
chk("1,060 of the 1,061 sextics have Galois group S6; the exception has 638 points", sum(r[17] == "S6" for r in SC) == 1060 and [(r[0], r[2]) for r in SC if r[17] != "S6"] == [("638", "p144")])
gain = sorted(I(r, 0) - I(r, 1) for r in SC)
chk("the closure adds between 6 and 116 points (median 42); the curves with >= 600 points have >= 512 points of height <= 10^5, those with >= 550 have >= 444",
    (gain[0], gain[-1], med(gain)) == (6, 116, 42) and min(I(r, 1) for r in top) == 512 and min(I(r, 1) for r in SC if I(r, 0) >= 550) == 444)
bysc = collections.Counter(I(r, 13) for r in SC)
chk("the sixteen curves have scores >= 46, >= 310 points of height <= 2000 and >= 400 of height <= 10^4; seven closed curves have score 40 and ten have 41, one of the former with 584 points",
    (min(I(r, 13) for r in top), min(I(r, 15) for r in top), min(I(r, 16) for r in top)) == (46, 310, 400) and (bysc[40], bysc[41]) == (7, 10) and max(I(r, 0) for r in SC if I(r, 13) == 40) == 584)
three = [r for r in SM if r[1] == "3"]; g = collections.defaultdict(list)
for r in three: g[I(r, 2)].append(I(r, 10))
chk("three-node models: curves with >= 400 points of height <= 10^5 per model fall from 2.4 (72 forced points) to 1.3 (84) and 0.8 (86); rank correlation -0.13",
    ["%.1f" % (sum(g[k]) / len(g[k])) for k in (72, 84, 86)] == ["2.4", "1.3", "0.8"] and [len(g[k]) for k in (72, 84, 86)] == [10, 66, 5] and "%.2f" % corr(ranks([I(r, 2) for r in three]), ranks([I(r, 10) for r in three])) == "-0.13")
# ---------------------------------------------------------------- the model of C_738
hist = {}; lines = []
for t in rd("S3_search.txt"):
    if t[0] == "score": hist[int(t[1])] = int(t[2])
    else: lines.append((tuple(int(x) for x in t[1:4]), int(t[4]), int(t[5]), None if t[6] == "-" else int(t[6]), None if t[7] == "-" else int(t[7])))
bys = sorted(lines, key=lambda l: -l[1]); rk = {l[0]: i + 1 for i, l in enumerate(bys)}; r0 = sm["t0453"]
chk("on the model of C_738: 1,328,965 surface points; 249,542 lines through no node with at least 15 of them; the 2,460 with at least 40 are kept",
    r0[3] == "1328965" and sum(hist.values()) == 249542 and sum(n for s_, n in hist.items() if s_ >= 40) == len(lines) == 2460 == I(r0, 5) and min(l[1] for l in lines) == 40 and I(r0, 6) == 112)
c4 = [l for l in lines if l[3] is not None and l[3] >= 330]
chk("... the line (66, 31, -52) has score 112 and 424, 536, 662 points; the best other line 95, 310, 346, 376; 27 curves reach 250 points, four reach 330, one reaches 400",
    bys[0] == ((66, 31, -52), 112, 424, 536, 662) and (max(l[1] for l in bys[1:]), max(l[2] for l in bys[1:]), max(l[3] for l in bys[1:] if l[3] is not None), max(l[4] for l in bys[1:] if l[4] is not None)) == (95, 310, 346, 376)
    and sum(l[2] >= 250 for l in lines) == 27 == I(r0, 8) and len(c4) == 4 == I(r0, 9) and sum(1 for l in lines if l[4] is not None and l[4] >= 400) == 1 == I(r0, 10) and r0[12:14] == ["662", "738"])
chk("... the four curves with >= 330 points of height <= 10^4 have ranks 1, 2, 34, 2105 by score; the last has score 41, 262 points of height <= 2000 and 372 of height <= 10^5; the model takes %s seconds" % r0[19],
    sorted(rk[l[0]] for l in c4) == [1, 2, 34, 2105] and max(c4, key=lambda l: rk[l[0]])[1:] == (41, 262, 330, 372) and intex("closure included, takes $%s$ seconds on two threads" % r0[19]) and 10 <= I(r0, 19) <= 60)
# ---------------------------------------------------------------- the sample searched with lower cuts
SAMP = rd("search_sample.txt"); smod = [ln for ln in open(os.path.join(D, "search_sample.txt")) if ln.startswith("# models:")][0].split()[2:]
lost = [r for r in SAMP if r[10] != "found"]
chk("a sample of 87 models (7, 20, 60 with one, two, three nodes) searched with score >= 30 and thresholds 200, 300: the curves marked 'found' are the rows of search_curves.txt for these models, with the same numbers",
    len(smod) == len(set(smod)) == 87 and [sum(nodes[m] == k for m in smod) for k in (1, 2, 3)] == [7, 20, 60] and all(r[2] in smod and nodes[r[2]] == int(r[3]) and int(r[1]) >= 440 and int(r[7]) >= 30 for r in SAMP)
    and sorted((r[2], r[4], r[5], r[6], r[0], r[1], r[7], r[8], r[9]) for r in SAMP if r[10] == "found") == sorted((c[2], c[3], c[4], c[5], c[0], c[1], c[13], c[15], c[16]) for c in SC if c[2] in smod)
    and all((r[10] == "score") == (int(r[7]) < 40) and (r[10] == "count2000") == (int(r[7]) >= 40 and int(r[8]) < 250) and (r[10] == "count1e4") == (int(r[7]) >= 40 and int(r[8]) >= 250 and int(r[9]) < 330) for r in lost))
SAMPLE = (152, 144, 7, 1, 0, 556, 8, 0)          # (curves closed in the sample, of them found by the search, lost with >= 500, >= 550, >= 600 points, most points among the lost, lost by the score, by the two counts)
chk("... %d curves closed in the sample, %d of them by the search; of the %d others, %d / %d / %d have at least 500 / 550 / 600 points (most: %d); %d are lost by the score bound and %d by the counts"
    % (SAMPLE[0], SAMPLE[1], SAMPLE[0] - SAMPLE[1], SAMPLE[2], SAMPLE[3], SAMPLE[4], SAMPLE[5], SAMPLE[6], SAMPLE[7]),
    (len(SAMP), len(SAMP) - len(lost), sum(int(r[0]) >= 500 for r in lost), sum(int(r[0]) >= 550 for r in lost), sum(int(r[0]) >= 600 for r in lost), max(int(r[0]) for r in lost),
     sum(r[10] == "score" for r in lost), sum(r[10] != "score" for r in lost)) == SAMPLE
    and (min(int(r[7]) for r in lost), max(int(r[7]) for r in lost), min(int(r[0]) for r in lost), sum(r[3] == "1" for r in lost)) == (30, 39, 494, 7)
    and intex("$152$ curves are closed, where the search above closes $144$", "scores between $30$ and $39$", "between $494$ and $556$ points"))
# ---------------------------------------------------------------- times
T = []
for ln in open(os.path.join(D, "search_summary.txt")):
    m = re.match(r"(\s*)(\S.*?)\s+([0-9.]+) min\s+([0-9.]+) h\s*$", ln)
    if m: T.append((len(m.group(1)) > 0, m.group(2), float(m.group(3)), float(m.group(4))))
tt = [w for sub, n_, w, c in T if n_ == "total"][0]; tm = [w for sub, n_, w, c in T if n_.startswith("building the models")][0]
LAB = {"building the models": "constructing the models", "enumerating points on the surfaces": "listing the points of the surfaces", "ranking lines (points on each line, selection)": "scoring the lines",
       "checking the curves of the lines kept": "the curves of the lines kept", "restriction to the line and reduction": "restriction and reduction", "point count to height 2000": "points of height at most $2000$",
       "point count to height 10^4": "points of height at most $10^4$", "point count to height 10^5": "points of height at most $10^5$", "closure": "closure", "total": "total"}
trows = ["%s%s & $%.1f$ & $%.1f$\\\\" % ("\\quad " if sub else "", LAB["building the models" if n_.startswith("building the models") else n_], w, c) for sub, n_, w, c in T if ("building the models" if n_.startswith("building the models") else n_) in LAB]
chk("the run took %.0f minutes, %.0f of them for the models; the table of times in the note is that of search_summary.txt (10 rows)" % (tt, tm),
    len(trows) == 10 and all(r in tex for r in trows) and intex("took $%.0f$ minutes on a $16$-core machine, $%.0f$ of them for the construction" % (tt, tm), "took $%.0f$ minutes on that machine" % tt, "construction takes $%.0f$ minutes" % tm)
    and 30 < tt < 120 and abs(sum(w for sub, n_, w, c in T if not sub and n_ != "total") - tt) < 0.3)
wm = [med([I(r, 19) for r in SM if int(r[1]) == k]) for k in (3, 2, 1)]
chk("a model takes about %d / %d / %d seconds on two threads (median, three / two / one nodes); the processor seconds of search_models.txt add up to the search rows of the table of times" % tuple(wm),
    intex("a three-node model takes about $%d$ seconds, a two-node model about\n$%d$ and a one-node model about $%d$" % tuple(wm))
    and all(abs(tot(c) / 3600 - h) < 0.06 for c, h in ((14, [c_ for sub, n_, w, c_ in T if n_.startswith("enumerating")][0]), (15, [c_ for sub, n_, w, c_ in T if n_.startswith("ranking")][0]), (18, [c_ for sub, n_, w, c_ in T if n_ == "closure"][0]))))
# ---------------------------------------------------------------- the logs of the Magma certificates
ml = open(os.path.join(HERE, "models.log")).read().split("\n")
mrow = {ln.split()[0]: ln for ln in ml if "NS of rank 20, discriminant -163" in ln}
chk("models.log: every sextic of models.txt is certified (NS of rank 20 and discriminant -163, the stated number of forced points), and the 2381 sextics are 2381 isomorphism classes",
    len(mrow) == 2381 and all(m in mrow and ("%d forced points" % (2 * nlin[m])) in mrow[m] and mrow[m].split()[1] == "%d:" % nodes[m] for m in name)
    and all(("%d nodes: %d sextics in %d isomorphism classes" % (k, len(byk[k]), len(byk[k]))) in ml for k in (1, 2, 3))
    and "MODELS CHECK PASSED: 2381 sextics, 2381 isomorphism classes" in ml and not any("the same model" in ln for ln in ml))
chk("... the model of C_654 has 26 tritangents and 16 node lines, the model of C_738 has 21 and 18", "t1272 3: 42 lines (26, " in mrow["t1272"] and "t0453 3: 39 lines (21, 15, 3, 0 through" in mrow["t0453"])
ll = [ln for ln in open(os.path.join(HERE, "linesection.log")).read().split("\n") if "lines over Qbar with a section isomorphic to C_738" in ln]
lrow = {ln.split()[0]: ln for ln in ll}
chk("linesection.log: no line on 2380 of the sextics, and exactly one, (66, 31, -52), on S_3",
    len(ll) == len(lrow) == 2381 and set(lrow) == set(name) and all(("0 + 0 + 0 in the three charts" in lrow[m]) == (m != "t0453") for m in name)
    and "1 + 0 + 0 in the three charts; rational: (66, 31, -52)" in lrow["t0453"] and "the model S_3: exactly one line over Qbar, the line (66, 31, -52), and S_3 on it is 64 f" in open(os.path.join(HERE, "linesection.log")).read())
bycls = collections.defaultdict(list)
for r in SC: bycls[r[18]].append(r[2])
chk("no curve of the search lies on two models; 26 curves occur twice on their sextic, 7 of them on Elkies' sextic, which is even in z0",
    all(len(set(fs)) == 1 for fs in bycls.values()) and sum(len(fs) > 1 for fs in bycls.values()) == 26 and max(len(fs) for fs in bycls.values()) == 2
    and sum(fs.count("E") > 1 for fs in bycls.values()) == 7
    and all(c == 0 for c, e in zip(coef["E"], [(6 - d) for d in range(7) for j in range(d + 1)]) if e % 2 == 1))
cl = open(os.path.join(HERE, "curves.log")).read()
chk("curves.log: 1061 rows verified, 1060 with Galois group S6, 1035 classes, none on two models", "1061 rows: each sextic is the section of its model by its line; 1060 have Galois group S6" in cl
    and "1035 classes of curves" in cl and "classes on more than one model: 0; classes that occur more than once on their model: 26 (7 of them on Elkies' sextic E)" in cl and "CURVES CHECK PASSED" in cl)
il = open(os.path.join(HERE, "inose", "verify.log")).read(); iv = [ln for ln in open(os.path.join(HERE, "inose", "verify_inverse.log")) if ln.startswith("fibres")]
chk("inose/verify.log and inose/verify_inverse.log: the map to Inose's surface and its inverse are verified (no failure, 100 fibres)",
    "INOSE MAP CHECK PASSED" in il and "FAILED" not in il and len(iv) == 20 and all(ln.rstrip().endswith("0 failures: ok") for ln in iv)
    and [ln.split()[1].rstrip(":") for ln in iv] == ["%d..%d" % (5 * i + 1, 5 * i + 5) for i in range(20)])
# data/census.txt is the list of the classes with one, two or three nodes of the census files of the project: rebuilt here and compared
if CENSUS:
    want = {}
    for r in rows:
        if r[7].strip() in ("A1", "A1^2", "A1^3"): want[r[11].strip()] = (int(r[0]), r[7].strip(), r[2].strip())
    cen = []
    for ln in gzip.open(os.path.join(os.path.dirname(PATH), "classes_merged.txt.gz"), "rt"):
        if ln.startswith("C ") and ln.split()[1] in want:
            cid = ln.split()[1]; g = ln[ln.index(" G [") + 4:].strip().rstrip("]").split(",")
            cen.append((["A1", "A1^2", "A1^3"].index(want[cid][1]), want[cid][0], "%s | %s | %s | %s" % (cid, want[cid][1], want[cid][2], " ".join(t.strip() for t in g))))
    cen.sort()
    chk("data/census.txt is the list of the %d classes with one, two or three nodes of the census (id, root type, forced points, Gram matrix), rebuilt from the census files and compared exactly" % len(cen),
        len(cen) == 2385 == len(want) and [c[2] for c in cen] == [ln.rstrip("\n") for ln in open(os.path.join(D, "census.txt")) if not ln.startswith("#")])
cs = open(os.path.join(HERE, "census.log")).read(); cm = re.search(r"A1\^3: 2192 classes in the census, 2188 of them are models of the list; not: \[ (\d+), (\d+), (\d+), (\d+) \]\n", cs)
crow = dict(ln.split() for ln in cs.split("\n") if len(ln.split()) == 2 and ln.split()[0] in nodes)
chk("census.log: the 2381 models are 2381 classes of the census: all 7 with one node, all 186 with two, and 2188 of the 2192 with three",
    "CENSUS CHECK PASSED: the 2381 models are 2381 classes of the census" in cs and "A1: 7 classes in the census, 7 of them are models of the list\n" in cs
    and "A1^2: 186 classes in the census, 186 of them are models of the list\n" in cs and cm is not None and set(crow) == set(name) and len(set(crow.values())) == 2381
    and not set(cm.groups()) & set(crow.values()))
N = 43 if CENSUS else 43 - N_CENSUS
ok = ok and ran == N
verdict = "ALL CHECKS PASSED (%d of %d)" % (ran, N) if CENSUS else "CHECKS PASSED (%d of %d; the other %d read the files of the census, which are not here)" % (ran, N, N_CENSUS)
print(verdict if ok else "SOME CHECKS FAILED (%d run, %d expected)" % (ran, N)); sys.exit(0 if ok else 1)

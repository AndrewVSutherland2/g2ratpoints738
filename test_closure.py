#!/usr/bin/env python3
"""test_closure.py -- regression test of closure.js, the program embedded in the note.

Runs closure.js in headless Firefox on a few curves and compares what it returns with
  * the 369 x-values of data/points_738.txt, from the five starting values of the note (the certificate of the note),
    together with the counts of sets tried and relations found that the note quotes;
  * an independent closure written here in Python, in homogeneous coordinates with exact rational arithmetic, on small
    curves that have rational points with y = 0, rational points at infinity (leading coefficient a square) or degree 5.
The reference treats every point (x : z) alike, so it does not share the special cases of closure.js.

usage: python3 test_closure.py        (needs firefox; about a minute)"""
import http.server, itertools, json, math, os, shutil, subprocess, sys, tempfile, threading
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
N_CHECKS = 33                           # the verdict requires exactly this many checks to have run and passed

# ---------------------------------------------------------------- the reference closure
def isqrt_exact(n):
    if n < 0: return None
    r = math.isqrt(n); return r if r * r == n else None

def F(f, a, c):                         # c^6 f(a/c); f = [f6, ..., f0]
    return sum(f[k] * a ** (6 - k) * c ** k for k in range(7))

def pmul(p, q):                         # polynomials in x (homogenised by z), ascending in x
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r

def solve(M, v):                        # Gaussian elimination over Q, square nonsingular system
    n = len(M); A = [[Fr(x) for x in row] + [Fr(y)] for row, y in zip(M, v)]
    for i in range(n):
        p = next(r for r in range(i, n) if A[r][i] != 0); A[i], A[p] = A[p], A[i]
        A[i] = [x / A[i][i] for x in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0: A[r] = [x - A[r][i] * y for x, y in zip(A[r], A[i])]
    return [A[i][n] for i in range(n)]

def norm(a, c):
    g = math.gcd(a, c); a //= g; c //= g
    if c < 0 or (c == 0 and a < 0): a, c = -a, -c
    return (a, c)

def reference_closure(f, seeds, limit=400):
    """all (a : c) reachable from the seeds (and the rational points at infinity) by the rule: the binary cubic through four
    known points, with any choice of the signs of y, meets the curve in two further points; add them when they are rational"""
    asc = [Fr(c) for c in reversed(f)]                                  # f ascending in x; F(x, z) = sum asc[k] x^k z^(6-k)
    pts = {}
    def add(a, c):
        a, c = norm(a, c); y = isqrt_exact(F(f, a, c))
        if y is None or (a, c) in pts: return False
        pts[(a, c)] = y; return True
    for a, c in seeds: assert add(a, c), ("seed not on the curve", a, c)
    add(1, 0)
    done = set(); grew = True
    while grew:
        grew = False
        for S in itertools.combinations(sorted(pts), 4):
            if S in done: continue
            done.add(S)
            rows = [[Fr(a) ** k * Fr(c) ** (3 - k) for k in range(4)] for a, c in S]      # p(a, c) = sum p_k a^k c^(3-k)
            for sg in itertools.product((1, -1), repeat=3):
                ys = [pts[S[0]]] + [s * pts[P] for s, P in zip(sg, S[1:])]
                p = solve(rows, ys)                                      # ascending in x
                G = [u - v for u, v in zip(pmul(p, p), asc)]             # p^2 - F, a binary sextic, ascending in x
                for a, c in S:                                           # divide by c x - a z (synthetic division, both cases)
                    if c == 0:                                           # root at infinity: the x^6 coefficient vanishes
                        assert G[-1] == 0; G = G[:-1]
                    else:
                        r = Fr(a, c); q = [Fr(0)] * (len(G) - 1); carry = Fr(0)
                        for k in range(len(G) - 1, 0, -1): carry = G[k] + carry * r; q[k - 1] = carry
                        assert G[0] + carry * r == 0; G = q
                while len(G) < 3: G.append(Fr(0))
                assert len(G) == 3
                q0, q1, q2 = G                                           # residual q2 x^2 + q1 x z + q0 z^2
                if q2 == 0 and q1 == 0: continue                         # both further points at infinity (or G = 0)
                den = math.lcm(q0.denominator, q1.denominator, q2.denominator)
                A, B, C = int(q2 * den), int(q1 * den), int(q0 * den)
                if A == 0: new = [(1, 0), (-C, B)]
                else:
                    d = isqrt_exact(B * B - 4 * A * C)
                    if d is None: continue
                    new = [(-B + d, 2 * A), (-B - d, 2 * A)]
                for a, c in new:
                    if add(a, c): grew = True
            assert len(pts) <= limit, "the reference closure grew past the limit; choose a smaller example"
    return pts

def key(a, c): return "∞" if c == 0 else (str(a) if c == 1 else "%d/%d" % (a, c))

# ---------------------------------------------------------------- the cases
def prod_roots(roots):
    p = [Fr(1)]
    for r in roots: p = pmul(p, [Fr(-r), Fr(1)])
    return p
def desc(asc, deg=6):
    asc = list(asc) + [Fr(0)] * (deg + 1 - len(asc)); assert all(x.denominator == 1 for x in asc)
    return [int(x) for x in reversed(asc)]
def sub(p, q):
    n = max(len(p), len(q)); p = list(p) + [Fr(0)] * (n - len(p)); q = list(q) + [Fr(0)] * (n - len(q))
    return [a - b for a, b in zip(p, q)]

C738 = [743425, -16144220, 60026500, 172162750, 40133975, -46160230, 6838225]
SMALL = [   # name, f (descending), starting x-values, what the case is for
    ("Weierstrass point", desc(sub([Fr(0), Fr(0), Fr(1)], prod_roots(range(6)))), [(0, 1), (1, 1), (2, 1), (3, 1)],
     "y^2 = x^2 - x(x-1)...(x-5): the line y = x through (0,0), (1,1), (2,2), (3,3) gives x = 4, 5"),
    ("square leading coefficient", desc(sub([Fr(0)] * 6 + [Fr(1)], prod_roots(range(1, 6)))), [(1, 1), (2, 1), (3, 1)],
     "y^2 = x^6 - (x-1)...(x-5): the cubic y = x^3 through a point at infinity and x = 1, 2, 3 gives x = 4, 5"),
    ("degree 5", desc(sub([Fr(0), Fr(0), Fr(1)], prod_roots(range(5)))), [(0, 1), (1, 1), (2, 1), (3, 1)],
     "y^2 = x^2 - x(x-1)...(x-4): the line y = x through four points gives x = 4 and the point at infinity"),
    ("two Weierstrass points", desc(sub([Fr(0), Fr(0), Fr(1), Fr(-2), Fr(1)], prod_roots(range(6)))), [(0, 1), (1, 1), (2, 1), (3, 1)],
     "y^2 = (x(x-1))^2 - x(x-1)...(x-5): y = x(x-1) through (0,0), (1,0), (2,2), (3,6) gives x = 4, 5"),
]

PAGE = """<!doctype html><meta charset="utf-8"><title>closure.js test</title><script>%s</script><script>
var CASES = %s, out = [];
CASES.forEach(function (c) {
  var r;
  try {
    var st = Closure.start(c.coef, c.seeds, { K: c.K }); while (!st.step(1000)) {} r = Closure.result(st); r.rejected = st.rejected;
    var S = function (q) { return q === null ? null : [q[0].toString(), q[1].toString()]; };
    r.rel = st.pts.map(function (pt, n) {     // for every point: null for a starting value, else the relation that produced it
      var x = st.explain(n); if (!x) return null;
      return { key: x.key, n: n, ok: x.ok, round: x.round, parents: x.parents.map(function (pp) { return { n: pp.n, key: pp.key, y: S(pp.y) }; }),
               cubic: x.cubic.map(S), others: x.others };
    });
    r.cert = Closure.certText(st);            // the same run as a certificate (the text of "Show the certificate")
  }
  catch (e) { r = { error: String(e) }; }
  out.push(r);
});
var x = new XMLHttpRequest(); x.open('POST', '/result', false); x.send(JSON.stringify(out));
</script>"""

def run_in_firefox(cases):
    js = open(os.path.join(HERE, "closure.js"), encoding="utf-8").read()
    return run_page(PAGE % (js, json.dumps(cases)))
def run_page(page):
    """Serves the page to headless Firefox and returns the JSON it posts back."""
    page = page.encode("utf-8"); got = {}
    class H(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a): pass
        def do_GET(self):
            self.send_response(200); self.send_header("Content-Type", "text/html; charset=utf-8"); self.end_headers(); self.wfile.write(page)
        def do_POST(self):
            got["body"] = self.rfile.read(int(self.headers["Content-Length"])); self.send_response(204); self.end_headers()
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H); threading.Thread(target=srv.serve_forever, daemon=True).start()
    prof = tempfile.mkdtemp(prefix="closure_test_")
    try:
        subprocess.run(["firefox", "--headless", "--no-remote", "--profile", prof, "--screenshot", os.path.join(prof, "shot.png"),
                        "http://127.0.0.1:%d/" % srv.server_address[1]], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)
    finally:
        srv.shutdown(); shutil.rmtree(prof, ignore_errors=True)
    if "body" not in got: sys.exit("FAILED: the page did not report (is firefox installed?)")
    return json.loads(got["body"].decode("utf-8"))

def parse_key(k):
    if k == "∞": return None
    a, _, c = k.partition("/"); return Fr(int(a), int(c or 1))
def check_relations(f, r, nseed_keys):
    """Independent exact check of every relation closure.js recorded (f = [f6..f0]): the cubic q passes through the finite
    parents with the stated y, q^2 - f equals its leading coefficient times the product of (x - x_i) over the finite parents
    and the finite new points, the missing degree is the number of points at infinity among them, the point itself is one
    of the new points, and every parent was found before the point.  Returns (number of relations, number that verify,
    largest number of steps back to a starting value)."""
    fa = [Fr(c) for c in reversed(f)]; ev = lambda q, x: sum(c * x ** k for k, c in enumerate(q))
    rel = r["rel"]; good = 0; depth = [0] * len(rel)
    for n, x in enumerate(rel):
        if x is None: continue
        q = [Fr(int(a), int(b)) for a, b in x["cubic"]]; ok = len(q) == 4 and x["key"] in x["others"] and x["ok"]
        fin = []
        for pp in x["parents"]:
            ok = ok and pp["n"] < n
            if pp["key"] == "∞": ok = ok and pp["y"] is None; continue
            xv = parse_key(pp["key"]); yv = Fr(int(pp["y"][0]), int(pp["y"][1]))
            ok = ok and yv * yv == ev(fa, xv) and ev(q, xv) == yv; fin.append(xv)
        newfin = [parse_key(k) for k in x["others"] if k != "∞"]
        ninf = sum(pp["key"] == "∞" for pp in x["parents"]) + sum(k == "∞" for k in x["others"])
        lhs = sub(pmul(q, q), fa)
        while len(lhs) > 1 and lhs[-1] == 0: lhs.pop()
        rhs = [lhs[-1]]
        for t in fin + newfin: rhs = pmul(rhs, [-t, Fr(1)])
        ok = ok and len(fin) == len(set(fin)) and lhs == rhs and len(lhs) - 1 == 6 - ninf and len(x["parents"]) == 4 and len(x["others"]) == 2
        good += bool(ok); depth[n] = 1 + max(depth[pp["n"]] for pp in x["parents"])
    return sum(x is not None for x in rel), good, max(depth)

def main():
    ok = True; ran = 0
    def chk(name, cond, detail=""):
        nonlocal ok, ran; ok &= bool(cond); ran += 1; print("%-4s %s%s" % ("ok" if cond else "FAIL", name, (" -- " + detail) if detail else ""))
    listed = set()
    for ln in open(os.path.join(HERE, "data", "points_738.txt")):          # one line "a y c" per x-value a/c, y > 0
        if ln.strip() and not ln.startswith("#"):
            a, y, c = map(int, ln.split()); assert y > 0 and y * y == F(C738, a, c); listed.add(key(a, c))
    cases = [{"coef": " ".join(map(str, C738)), "seeds": "0, -1/2, 1/2, -3/2, 1/3", "K": 16}]
    refs = []
    for name, f, seeds, why in SMALL:
        ref = reference_closure(f, seeds); refs.append(ref)
        cases.append({"coef": " ".join(map(str, f)), "seeds": ", ".join(key(a, c) for a, c in seeds), "K": 0})
    res = run_in_firefox(cases)
    r = res[0]
    chk("C_738: no error", "error" not in r, r.get("error", ""))
    if "error" not in r:
        chk("C_738: the five starting values give exactly the 369 listed x-values", len(listed) == 369 and set(r["xs"]) == listed and len(r["xs"]) == 369)
        chk("C_738: 738 rational points, none at infinity", r["points"] == 738 and not r["infinity"])
        chk("C_738: 7,805,120 sets of four tried and 84,317 relations found, as quoted in the note", r["tests"] == 7805120 and r["relations"] == 84317,
            "%d sets, %d relations" % (r["tests"], r["relations"]))
        nrel, good, depth = check_relations(C738, r, 5)
        chk("C_738: each of the 364 x-values that is not a starting value has a recorded relation", nrel == 364 and sum(x is None for x in r["rel"]) == 5, "%d relations" % nrel)
        chk("C_738: all 364 recorded relations verify exactly (cubic through the four parents, q^2 - f factors as stated)", good == 364, "%d of %d; at most %d steps back to a starting value" % (good, nrel, depth))
        import check_witnesses as CW, copy
        shipped = open(os.path.join(HERE, "data", "witnesses_738.json"), encoding="utf-8").read()
        chk("C_738: the certificate text is data/witnesses_738.json, byte for byte", r["cert"] == shipped)
        cert = json.loads(r["cert"]); known, npts, nw = CW.check(cert)
        chk("C_738: check_witnesses.py verifies the certificate: 364 witnesses, the 369 listed x-values, 738 points", (nw, npts) == (364, 738) and set(known) == listed)
        def why(c):                           # the reason check_witnesses.py gives for rejecting c ("" if it accepts)
            try: CW.check(c); return ""
            except AssertionError as e: return str(e)
        W = cert["witnesses"]
        bad1 = copy.deepcopy(cert); bad1["witnesses"][100]["x"] = "7"
        bad2 = copy.deepcopy(cert); y = W[200]["parents"][1][1]; bad2["witnesses"][200]["parents"][1][1] = y[1:] if y.startswith("-") else "-" + y
        late = max(range(len(W)), key=lambda n: sum(k not in cert["start"] for k, _ in W[n]["parents"]))
        bad3 = copy.deepcopy(cert); bad3["witnesses"].insert(0, bad3["witnesses"].pop(late))
        bad4 = copy.deepcopy(cert); bad4["start"][1] = "7"
        bad5 = copy.deepcopy(cert); bad5["start"].append("0/1")                      # the same x-value under another name
        bad6 = copy.deepcopy(cert); bad6["witnesses"].append(copy.deepcopy(W[0])); bad6["witnesses"][-1]["x"] = "80/482"; bad6["witnesses"][-1]["others"][0] = "80/482"
        resp = copy.deepcopy(cert); resp["start"][1] = "-2/4"                        # other spellings of known values are the same values
        for w in resp["witnesses"]: w["parents"] = [["-3/6" if k == "-1/2" else k, y] for k, y in w["parents"]]
        resp["witnesses"][0]["x"] = "80/482"
        chk("C_738: a certificate with a wrong x, a wrong sign of y, a parent that is not known yet, or a starting value off the curve is rejected",
            "x is not one of the further points" in why(bad1) and "is not on the cubic" in why(bad2) and "is not known yet" in why(bad3)
            and "starting value 7 is not on the curve" in why(bad4))
        assert W[0]["x"] == "40/241" and W[0]["others"][0] == "40/241"
        k2, n2, w2 = CW.check(resp)
        chk("C_738: x-values are compared as numbers: a value repeated under another spelling (0/1 as a start, 80/482 as a witness) is rejected, and respelling known values changes nothing",
            "repeated starting value 0" in why(bad5) and "x was known already" in why(bad6) and (set(k2), n2, w2) == (listed, 738, 364) and len(k2) == 369)
    for (name, f, seeds, why), ref, r in zip(SMALL, refs, res[1:]):
        if "error" in r: chk(name + ": no error", False, r["error"]); continue
        want = {key(a, c) for a, c in ref}; npts = sum(1 if y == 0 else 2 for y in ref.values())
        chk("%s: same %d x-values as the reference closure (%s)" % (name, len(want), why), set(r["xs"]) == want and len(r["xs"]) == len(want),
            "" if set(r["xs"]) == want else "closure.js only: %s; reference only: %s" % (sorted(set(r["xs"]) - want), sorted(want - set(r["xs"]))))
        chk("%s: %d rational points" % (name, npts), r["points"] == npts, "closure.js says %d" % r["points"])
        chk("%s: every x-value is on the curve" % name, all(isqrt_exact(F(f, a, c)) is not None for a, c in ref))
        nrel, good, depth = check_relations(f, r, len(seeds))
        chk("%s: every recorded relation verifies exactly" % name, nrel >= 1 and good == nrel, "%d of %d" % (good, nrel))
        import check_witnesses as CW
        try: known, cpts, nw = CW.check(json.loads(r["cert"])); cok = set(known) == want and cpts == npts and nw == nrel
        except AssertionError as e: cok = False
        chk("%s: check_witnesses.py verifies its certificate and finds the same points" % name, cok)
    chk("Weierstrass point: x = 4 and 5 are found", {"4", "5"} <= set(res[1].get("xs", [])))
    chk("square leading coefficient: x = 4, 5 and the points at infinity are found", {"4", "5", "∞"} <= set(res[2].get("xs", [])))
    chk("degree 5: x = 4 and the point at infinity are found", {"4", "∞"} <= set(res[3].get("xs", [])))
    ok = ok and ran == N_CHECKS            # a case that errors runs fewer checks: that is a failure too
    print("ALL CLOSURE TESTS PASSED (%d of %d)" % (ran, N_CHECKS) if ok else "SOME CLOSURE TESTS FAILED (%d checks run, %d expected)" % (ran, N_CHECKS)); sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()

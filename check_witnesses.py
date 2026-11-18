#!/usr/bin/env python3
"""check_witnesses.py -- verifies a certificate written by closure.js (Closure.certText, or the button "Show the certificate"
on the page): a list of rational points of a genus 2 curve y^2 = f(x), each one either a starting value or explained by one
relation, a cubic through four points known before.  Nothing is searched for; every step is an exact identity over Q.

The file is JSON: "curve" (the coefficients of f from x^6 down to the constant term), "start" (the starting x-values, "∞"
for the points at infinity) and "witnesses", one per further x-value, in an order in which every parent is a starting
value or the x of an earlier witness:
    "x"        the x-value explained,
    "parents"  four pairs [x, y] (y is null for a point at infinity),
    "cubic"    the coefficients q0, q1, q2, q3 of the cubic y = q(x),
    "others"   the two further points in which the cubic meets the curve.
For each witness the script checks that
  * the four parents are known, have distinct x, lie on the curve with the stated y, and lie on the cubic;
  * q(x)^2 - f(x) = c * prod (x - x_i), the product over the finite parents and the finite further points, as an identity
    of polynomials, and the missing degree is the number of points at infinity among the six;
  * x is one of the further points.
Then (x, q(x)) is a rational point of the curve.  A point at infinity is accepted as a starting value only if the curve has
rational points there (degree 5, or a square leading coefficient).

usage: python3 check_witnesses.py [FILE]        default: data/witnesses_738.json, for which the result is also compared
                                                with the list data/points_738.txt"""
import json, os, sys
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
INF = "∞"

def rat(s):                                               # "a/b" or "a"
    try: return Fr(s)
    except (ValueError, ZeroDivisionError, TypeError): raise AssertionError("%r is not a rational number" % (s,))
def canon(k):                                             # one name per x-value: "0/1", "2/4" and "1/2 " must not count as new values
    return INF if k == INF else str(rat(k))
def ev(p, x): return sum(c * x ** k for k, c in enumerate(p))
def pmul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r
def is_square(v):                                         # of a rational number
    return v >= 0 and isqrt(v.numerator) ** 2 == v.numerator and isqrt(v.denominator) ** 2 == v.denominator

def check(cert):
    """Returns (list of x-values proved to be x-coordinates of rational points, number of rational points, number of witnesses);
    raises AssertionError with a message at the first failure.  x-values are compared as rational numbers, whatever their
    spelling in the file, and are returned in lowest terms."""
    assert len(cert["curve"]) == 7, "f must be given by seven coefficients"
    f = [rat(c) for c in reversed(cert["curve"])]         # ascending
    assert all(c.denominator == 1 for c in f) and (f[6] != 0 or f[5] != 0), "f must have integer coefficients and degree 5 or 6"
    inf_ok = f[6] == 0 or is_square(f[6])
    known = []; seen = set()
    for k in cert["start"]:
        k = canon(k)
        assert k not in seen, "repeated starting value %s" % k
        if k == INF: assert inf_ok, "the curve has no rational point at infinity"
        else: assert is_square(ev(f, rat(k))), "starting value %s is not on the curve" % k
        seen.add(k); known.append(k)
    for n, w in enumerate(cert["witnesses"]):
        tag = "witness %d (x = %s): " % (n + 1, w["x"])
        assert len(w["cubic"]) == 4 and len(w["parents"]) == 4 and len(w["others"]) == 2, tag + "malformed"
        q = [rat(c) for c in w["cubic"]]; x = canon(w["x"]); others = [canon(k) for k in w["others"]]
        fin = []; ninf = 0
        for k, y in w["parents"]:
            k = canon(k)
            assert k in seen, tag + "parent %s is not known yet" % k
            if k == INF: assert y is None, tag + "a point at infinity has no y"; ninf += 1; continue
            xv = rat(k); yv = rat(y)
            assert yv * yv == ev(f, xv), tag + "parent %s is not on the curve with the stated y" % k
            assert ev(q, xv) == yv, tag + "parent %s is not on the cubic" % k
            fin.append(xv)
        assert len(set(fin)) == len(fin) and ninf <= 1 + (f[6] != 0), tag + "the parents are not distinct"
        new = [rat(k) for k in others if k != INF]; ninf += sum(k == INF for k in others)
        lhs = [a - b for a, b in zip(pmul(q, q), f)]
        while len(lhs) > 1 and lhs[-1] == 0: lhs.pop()
        assert lhs != [0], tag + "q^2 = f"
        rhs = [lhs[-1]]
        for t in fin + new: rhs = pmul(rhs, [-t, Fr(1)])
        assert lhs == rhs, tag + "q^2 - f is not the stated product"
        assert len(lhs) - 1 == 6 - ninf, tag + "wrong number of points at infinity"
        assert x in others, tag + "x is not one of the further points"
        if x == INF: assert inf_ok
        else: assert is_square(ev(f, rat(x))), tag + "x is not on the curve"           # implied by the identity; checked anyway
        assert x not in seen, tag + "x was known already"
        seen.add(x); known.append(x)
    npts = 0
    for k in known:
        if k == INF: npts += 1 if f[6] == 0 else 2
        else: npts += 1 if ev(f, rat(k)) == 0 else 2
    return known, npts, len(cert["witnesses"])

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "witnesses_738.json")
    cert = json.load(open(path, encoding="utf-8"))
    try: known, npts, nw = check(cert)
    except AssertionError as e: print("CERTIFICATE REJECTED: %s" % e); sys.exit(1)
    print("curve: y^2 = f(x), f = %s (from x^6 down)" % " ".join(cert["curve"]))
    print("%d starting values, %d witnesses, all verified exactly" % (len(cert["start"]), nw))
    if len(sys.argv) <= 1:
        listed = set()
        for ln in open(os.path.join(HERE, "data", "points_738.txt")):
            if ln.strip() and not ln.startswith("#"):
                a, y, c = map(int, ln.split()); listed.add(str(Fr(a, c)))
        same = [int(c) for c in cert["curve"]] == [743425, -16144220, 60026500, 172162750, 40133975, -46160230, 6838225] and set(known) == listed and len(known) == 369
        print("the x-values are exactly the 369 of data/points_738.txt: %s" % ("yes" if same else "NO"))
        if not same or [canon(k) for k in cert["start"]] != ["0", "-1/2", "1/2", "-3/2", "1/3"]: print("CERTIFICATE DOES NOT MATCH THE NOTE"); sys.exit(1)
    print("CERTIFICATE VERIFIED: %d x-values, %d rational points" % (len(known), npts))

if __name__ == "__main__": main()

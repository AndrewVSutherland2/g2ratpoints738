/* closure.js -- four-point cubic closure on a genus 2 curve y^2 = f(x), deg f = 5 or 6, in exact arithmetic.
   Given some rational points (by their x-coordinates), repeatedly takes four known points, passes the cubic y = p(x)
   through them and looks at the two further intersection points with the curve: their x-coordinates are the roots of
   the quadratic r in  p^2 - f = (p3^2 - f6) (x-x1)(x-x2)(x-x3)(x-x4) r(x).  When the discriminant of r is a rational
   square the two points are rational and are added.  Every point that is reported has been verified directly:
   c^6 f(a/c) is checked to be a perfect square with BigInt arithmetic.  The modular prefilter (the discriminant must be
   a square modulo a dozen primes) only decides which relations are worth checking exactly.
   Points with y = 0 are used like any others (they have no sign to vary).  Rational points at infinity exist when
   the leading coefficient is a square (two points) or f has degree 5 (one point); they are then known points from the
   start, shown as "∞", and the relations through them are found as well: the program moves them to a finite place with
   the substitution u = 1/(x - t), for an integer t at which f(t) is neither 0 nor a square, and runs on the model
   v^2 = u^6 f(t + 1/u), which has no rational points at infinity.  x-values and heights are always those of the given model.
   Option K: each relation must use at least two of the K known points of smallest height (K = 0: no restriction,
   every set of four is tried; that is the full closure and is much slower).
   Every point that is not a starting value remembers the relation that first produced it (its four parents and their
   signs); Closure.explain returns that relation in the coordinates of the given model, as the cubic y = q(x), and
   checks the identity q^2 - f = c * prod (x - x_i) over the four parents and the two new points exactly.
   Closure.certText returns all these relations as one text, a certificate that check_witnesses.py verifies. */
var Closure = (function () {
  'use strict';
  function babs(a) { return a < 0n ? -a : a; }
  function bgcd(a, b) { a = babs(a); b = babs(b); while (b) { var t = a % b; a = b; b = t; } return a; }
  function bsqrt(n) {                      // exact square root of a BigInt, or null
    if (n < 0n) return null; if (n < 2n) return n;
    var x = 1n << BigInt((n.toString(2).length + 1) >> 1);
    for (;;) { var y = (x + n / x) >> 1n; if (y >= x) break; x = y; }
    return x * x === n ? x : null;
  }
  function parseRat(s) {                   // "a/b" or "a" -> [num, den] in lowest terms, den > 0
    s = s.trim().replace(/−/g, '-'); var m = /^(-?\d+)(?:\s*\/\s*(\d+))?$/.exec(s);
    if (!m) throw new Error('cannot read the number "' + s + '"');
    var a = BigInt(m[1]), c = m[2] ? BigInt(m[2]) : 1n; if (c === 0n) throw new Error('zero denominator in "' + s + '"');
    var g = bgcd(a, c); if (g > 1n) { a /= g; c /= g; } return [a, c];
  }
  function ratStr(a, c) { return c === 1n ? a.toString() : a.toString() + '/' + c.toString(); }
  function rnorm(n, d) { if (d < 0n) { n = -n; d = -d; } var g = bgcd(n, d); if (g > 1n) { n /= g; d /= g; } return [n, d]; }
  function radd(p, q) { return rnorm(p[0] * q[1] + q[0] * p[1], p[1] * q[1]); }
  function rsub(p, q) { return rnorm(p[0] * q[1] - q[0] * p[1], p[1] * q[1]); }
  function rmul(p, q) { return rnorm(p[0] * q[0], p[1] * q[1]); }
  function rdiv(p, q) { return rnorm(p[0] * q[1], p[1] * q[0]); }
  var PRIMES = [46337, 46327, 46309, 46307, 46301, 46279, 46273, 46271, 46261, 46237, 46229, 46219, 46199, 46187];

  function evalF(f, a, c) {                // c^6 f(a/c) = sum f[k] a^(6-k) c^k, f[0] the coefficient of x^6
    var s = 0n, cp = 1n; for (var k = 0; k < 7; k++) { s = s * a + f[k] * cp; cp *= c; } return s;
  }
  function moved(f, t) {                   // coefficients (u^6 .. 1) of u^6 f(t + 1/u) = sum_k f_k (t u + 1)^k u^(6-k)
    var g = [0n, 0n, 0n, 0n, 0n, 0n, 0n], pw = [1n], k, d;      // g[d] is the coefficient of u^d; pw = (t u + 1)^k, ascending
    for (k = 0; k <= 6; k++) {
      for (d = 0; d <= k; d++) g[d + 6 - k] += f[6 - k] * pw[d];
      var nx = [0n]; for (d = 0; d <= k; d++) { nx[d] += pw[d]; nx.push(t * pw[d]); } pw = nx;
    }
    return g.reverse();
  }
  function signs(wi, wj, wk, wl) {         // sign patterns worth trying when some of the four points have y = 0 (argument: is y = 0?)
    var w = 255; if (wj) w &= 0x55; if (wk) w &= 0x33; if (wl) w &= 0x0F;          // such a point has no sign to vary,
    if (wi) { if (!wj) w &= 0x55; else if (!wk) w &= 0x33; else if (!wl) w &= 0x0F; }   // and the first point with y != 0 is taken with +
    return w;
  }
  function State(coefDesc, opts) {         // coefDesc: BigInt f6..f0
    this.given = coefDesc.slice(); this.K = (opts && opts.K) | 0; this.t = null; this.nweier = 0;
    if (bsqrt(coefDesc[0]) !== null) {     // leading coefficient 0 or a square: rational points at infinity, move them to u = 0
      for (var n = 0; this.t === null; n++) {
        if (n > 2000) throw new Error('f takes square values at too many integers: is it squarefree?');
        var t = BigInt(n & 1 ? (n + 1) >> 1 : -(n >> 1)), v = evalF(coefDesc, t, 1n);
        if (v !== 0n && bsqrt(v) === null) this.t = t;
      }
      coefDesc = moved(coefDesc, this.t);
    }
    this.f = coefDesc.slice();
    this.pts = []; this.index = new Map(); this.rejected = []; this.via = null;
    this.round = 0; this.tests = 0; this.exact = 0; this.relations = 0; this.done = false;
    this.mods = PRIMES.map(function (p) {
      var qr = new Uint8Array(p), inv = new Int32Array(p), i;
      for (i = 0; i < p; i++) qr[(i * i) % p] = 1;
      inv[1] = 1; for (i = 2; i < p; i++) inv[i] = (p - Math.floor(p / i) * inv[p % i] % p) % p;
      var fm = coefDesc.map(function (c) { return Number(((c % BigInt(p)) + BigInt(p)) % BigInt(p)); });
      return { p: p, qr: qr, inv: inv, f6: fm[0], f5: fm[1], f4: fm[2], x: [], y: [], bad: [] };
    });
    this.cursor = null;
  }
  /* add the point with x = a/c (given model; c = 0 is the point at infinity) if it is on the curve and new; returns 1 if added */
  State.prototype.addX = function (a, c, isSeed) {
    var key = c === 0n ? '∞' : ratStr(a, c); if (this.index.has(key)) return 0;
    var xa = a, h = c === 0n ? 1n : (babs(a) > c ? babs(a) : c);
    if (this.t !== null) {                 // the coordinate the closure runs in: u = c / (a - t c)
      var d = a - this.t * c; if (d === 0n) { if (isSeed) this.rejected.push(key); return 0; }    // x = t is not on the curve
      a = d < 0n ? -c : c; c = d < 0n ? -d : d;
    }
    var b = c === 0n ? null : bsqrt(evalF(this.f, a, c));
    if (b === null) { if (isSeed) this.rejected.push(key); return 0; }
    var i = this.pts.length;
    this.pts.push({ a: a, c: c, b: b, h: h, key: key, xa: xa, born: this.round, weier: b === 0n, from: this.via });
    this.index.set(key, i); if (b === 0n) this.nweier++;
    for (var m = 0; m < this.mods.length; m++) {
      var M = this.mods[m], P = BigInt(M.p), cm = Number(((c % P) + P) % P);
      if (cm === 0) { M.bad.push(1); M.x.push(0); M.y.push(0); continue; }
      var ci = M.inv[cm], am = Number(((a % P) + P) % P), bm = Number(b % P);
      M.bad.push(0); M.x.push(am * ci % M.p); M.y.push(bm * (ci * ci % M.p) % M.p * ci % M.p);
    }
    return 1;
  };
  /* modular test: bit e of the result is set when sign pattern e (signs of y_j, y_k, y_l; y_i is +) survives prime m */
  State.prototype.mask = function (m, i, j, k, l, want) {
    var M = this.mods[m], p = M.p, X = M.x, Y = M.y, INV = M.inv;
    if (M.bad[i] | M.bad[j] | M.bad[k] | M.bad[l]) return want;
    var xi = X[i], xj = X[j], xk = X[k], xl = X[l];
    var dij = INV[(xi - xj + p) % p], dik = INV[(xi - xk + p) % p], dil = INV[(xi - xl + p) % p];
    var djk = INV[(xj - xk + p) % p], djl = INV[(xj - xl + p) % p], dkl = INV[(xk - xl + p) % p];
    if (!(dij && dik && dil && djk && djl && dkl)) return want;
    var Ti = Y[i] * (dij * dik % p * dil % p) % p;
    var Tj = Y[j] * ((p - dij) * djk % p * djl % p) % p;
    var Tk = Y[k] * (dik * djk % p * dkl % p) % p;
    var Tl = Y[l] * ((p - dil) * djl % p * dkl % p) % p;
    var e1 = (xi + xj + xk + xl) % p, e2 = (xi * xj + xi * xk + xi * xl + xj * xk + xj * xl + xk * xl) % p;
    var ai = (e1 - xi + p) % p, aj = (e1 - xj + p) % p, ak = (e1 - xk + p) % p, al = (e1 - xl + p) % p;
    var Ui = Ti * ai % p, Uj = Tj * aj % p, Uk = Tk * ak % p, Ul = Tl * al % p;
    var Vi = Ti * ((e2 - xi * ai % p + p) % p) % p, Vj = Tj * ((e2 - xj * aj % p + p) % p) % p;
    var Vk = Tk * ((e2 - xk * ak % p + p) % p) % p, Vl = Tl * ((e2 - xl * al % p + p) % p) % p;
    var out = 0, QR = M.qr, f6 = M.f6, f5 = M.f5, f4 = M.f4;
    for (var e = 0; e < 8; e++) {
      if (!((want >> e) & 1)) continue;
      var p3 = Ti + ((e & 1) ? p - Tj : Tj) + ((e & 2) ? p - Tk : Tk) + ((e & 4) ? p - Tl : Tl);
      var q2 = Ui + ((e & 1) ? p - Uj : Uj) + ((e & 2) ? p - Uk : Uk) + ((e & 4) ? p - Ul : Ul);   // p2 = -q2
      var p1 = Vi + ((e & 1) ? p - Vj : Vj) + ((e & 2) ? p - Vk : Vk) + ((e & 4) ? p - Vl : Vl);
      p3 %= p; q2 %= p; p1 %= p;
      var A = (p3 * p3 - f6) % p; if (A < 0) A += p;
      if (A === 0) { out |= 1 << e; continue; }
      var S = (f5 + 2 * p3 * q2 - A * e1) % p; if (S < 0) S += p;            // A s = -A e1 - (2 p3 p2 - f5), p2 = -q2
      var T = (2 * p3 * p1 + q2 * q2 - f4 - A * e2 - e1 * S) % p; if (T < 0) T += p;
      var D = (S * S - 4 * (A * T % p)) % p; if (D < 0) D += p;
      if (QR[D]) out |= 1 << e;
    }
    return out;
  };
  /* exact check of one relation; adds the new points it gives */
  State.prototype.exactCheck = function (i, j, k, l, e) {
    this.exact++;
    var P = this.pts, idx = [i, j, k, l], sg = [1n, (e & 1) ? -1n : 1n, (e & 2) ? -1n : 1n, (e & 4) ? -1n : 1n];
    var x = idx.map(function (t) { return [P[t].a, P[t].c]; });
    var y = idx.map(function (t, n) { return rnorm(sg[n] * P[t].b, P[t].c ** 3n); });
    var e1 = [0n, 1n], e2 = [0n, 1n], m, t;
    for (m = 0; m < 4; m++) { e1 = radd(e1, x[m]); for (t = m + 1; t < 4; t++) e2 = radd(e2, rmul(x[m], x[t])); }
    var p3 = [0n, 1n], p2 = [0n, 1n], p1 = [0n, 1n];
    for (m = 0; m < 4; m++) {
      var w = [1n, 1n]; for (t = 0; t < 4; t++) if (t !== m) w = rmul(w, rsub(x[m], x[t]));
      var T = rdiv(y[m], w), a = rsub(e1, x[m]);
      p3 = radd(p3, T); p2 = rsub(p2, rmul(T, a)); p1 = radd(p1, rmul(T, rsub(e2, rmul(x[m], a))));
    }
    var f = this.f, A = rsub(rmul(p3, p3), [f[0], 1n]); if (A[0] === 0n) return 0;
    var S = rsub(rsub([0n, 1n], rmul(A, e1)), rsub(rmul([2n, 1n], rmul(p3, p2)), [f[1], 1n]));
    var Tt = rsub(rsub(radd(rmul([2n, 1n], rmul(p3, p1)), rmul(p2, p2)), [f[2], 1n]), radd(rmul(A, e2), rmul(e1, S)));
    var D = rsub(rmul(S, S), rmul([4n, 1n], rmul(A, Tt)));
    var rn = bsqrt(D[0]), rd = bsqrt(D[1]); if (rn === null || rd === null) return 0;
    this.relations++;
    var r = [rn, rd], twoA = rmul([2n, 1n], A), added = 0;
    var r1 = rdiv(radd(S, r), twoA), r2 = rdiv(rsub(S, r), twoA);
    this.via = { p: idx, e: e };          // a point added now is remembered as coming from this relation
    added += this.addU(r1[0], r1[1]); added += this.addU(r2[0], r2[1]);
    this.via = null;
    return added;
  };
  State.prototype.addU = function (n, d) {  // a root u = n/d in the running coordinate: back to x = t + 1/u = (t n + d)/n
    if (this.t === null) return this.addX(n, d, false);
    if (n === 0n) return this.addX(1n, 0n, false);
    var a = this.t * n + d; return n < 0n ? this.addX(-a, -n, false) : this.addX(a, n, false);
  };
  /* one slice of work (about `ms` milliseconds); returns true when the closure is complete */
  State.prototype.step = function (ms) {
    if (this.done) return true;
    var t0 = Date.now(), P = this.pts, n, c = this.cursor;
    for (;;) {
      if (!c) {                          // start a round on the points known now
        n = P.length; if (n < 4) { this.done = true; return true; }
        var order = []; for (var q = 0; q < n; q++) order.push(q);
        order.sort(function (u, v) { return P[u].h < P[v].h ? -1 : P[u].h > P[v].h ? 1 : u - v; });
        var isA = new Uint8Array(n), K = this.K > 0 ? Math.min(this.K, n) : n;
        for (q = 0; q < K; q++) isA[order[q]] = 1;
        var anchors = []; for (q = 0; q < n; q++) if (isA[q]) anchors.push(q);
        var isNew = new Uint8Array(n), anyNew = 0; for (q = 0; q < n; q++) if (P[q].born === this.round) { isNew[q] = 1; anyNew = 1; }
        if (!anyNew) { this.done = true; return true; }
        this.round++;
        c = this.cursor = { n: n, isA: isA, anchors: anchors, isNew: isNew, ai: 0, aj: 1, added: 0 };
      }
      var A = c.anchors, na = A.length, isA2 = c.isA, isNew2 = c.isNew; n = c.n;
      while (c.ai < na - 1) {
        var i = A[c.ai], j = A[c.aj], newij = isNew2[i] | isNew2[j];
        for (var k = 0; k < n; k++) {
          if (k === i || k === j || (isA2[k] && k < j)) continue;
          var newk = newij | isNew2[k];
          for (var l = k + 1; l < n; l++) {
            if (l === i || l === j || (isA2[l] && l < j)) continue;
            if (!(newk | isNew2[l])) continue;
            this.tests++;
            var w = this.mask(0, i, j, k, l, this.nweier ? signs(P[i].weier, P[j].weier, P[k].weier, P[l].weier) : 255);
            for (var m = 1; w && m < this.mods.length; m++) w = this.mask(m, i, j, k, l, w);
            if (w) for (var e = 0; e < 8; e++) if ((w >> e) & 1) c.added += this.exactCheck(i, j, k, l, e);
          }
        }
        if (++c.aj >= na) { c.ai++; c.aj = c.ai + 1; }
        if (Date.now() - t0 > ms) return false;
      }
      c = this.cursor = null;            // round finished; the next loop iteration starts a new one or stops
    }
  };
  /* polynomials with rational coefficients, ascending: [[n, d], ...] */
  function pmul(u, v) {
    var r = [], i, j; for (i = 0; i < u.length + v.length - 1; i++) r.push([0n, 1n]);
    for (i = 0; i < u.length; i++) for (j = 0; j < v.length; j++) r[i + j] = radd(r[i + j], rmul(u[i], v[j]));
    return r;
  }
  function padd(u, v) { var r = [], i; for (i = 0; i < Math.max(u.length, v.length); i++) r.push(radd(u[i] || [0n, 1n], v[i] || [0n, 1n])); return r; }
  function ptrim(u) { u = u.slice(); while (u.length > 1 && u[u.length - 1][0] === 0n) u.pop(); return u; }
  /* the relation that produced point number n: the cubic y = q(x) of the given model through its four parents, the two
     further intersection points, and the exact check of  q^2 - f = c * prod (x - x_i)  (the product runs over the finite
     ones among the four parents and the two new points; a point at infinity lowers the degree instead) */
  State.prototype.explain = function (n) {
    var P = this.pts, me = P[n], v = me.from; if (!v) return null;
    var self = this, sg = [1n, (v.e & 1) ? -1n : 1n, (v.e & 2) ? -1n : 1n, (v.e & 4) ? -1n : 1n];
    var x = v.p.map(function (t) { return [P[t].a, P[t].c]; });
    var y = v.p.map(function (t, m) { return rnorm(sg[m] * P[t].b, P[t].c ** 3n); });
    var p = [[0n, 1n]], m, t;                // the interpolating cubic in the coordinate the closure runs in
    for (m = 0; m < 4; m++) {
      var num = [[1n, 1n]], den = [1n, 1n];
      for (t = 0; t < 4; t++) if (t !== m) { num = pmul(num, [rsub([0n, 1n], x[t]), [1n, 1n]]); den = rmul(den, rsub(x[m], x[t])); }
      p = padd(p, pmul(num, [rdiv(y[m], den)]));
    }
    while (p.length < 4) p.push([0n, 1n]);
    var q = p;
    if (this.t !== null) {                   // back to the given model: y = v (x - t)^3, so q(x) = sum p_k (x - t)^(3-k)
      var lin = [[-this.t, 1n], [1n, 1n]], pw = [[1n, 1n]]; q = [[0n, 1n]];
      for (m = 3; m >= 0; m--) { q = padd(q, pmul(pw, [p[m]])); pw = pmul(pw, lin); }
    }
    while (q.length < 4) q.push([0n, 1n]);
    var parents = v.p.map(function (t, m) {  // each parent with its y in the given model
      var pt = P[t], yy = null;
      if (pt.key !== '∞') yy = self.t !== null ? rnorm(sg[m] * pt.b, pt.a ** 3n) : rnorm(sg[m] * pt.b, pt.c ** 3n);
      return { n: t, key: pt.key, y: yy, seed: !pt.from };
    });
    var lhs = padd(pmul(q, q), this.given.slice().reverse().map(function (c) { return [-c, 1n]; }));
    lhs = ptrim(lhs);
    var roots = [], others = [];
    parents.forEach(function (pp) { if (pp.key !== '∞') roots.push(parseRat(pp.key)); });
    var ndeg = lhs.length - 1 - roots.length;         // the number of finite new intersection points
    var rest = lhs, k;
    roots.forEach(function (r) {            // divide by (x - r), exactly
      var out = [], carry = [0n, 1n];
      for (k = rest.length - 1; k >= 1; k--) { carry = radd(rest[k], rmul(carry, r)); out.unshift(carry); }
      if (radd(rest[0], rmul(carry, r))[0] !== 0n) throw new Error('a parent is not on the cubic');
      rest = out;
    });
    var lead = rest[rest.length - 1], ok = rest.length - 1 === ndeg && ndeg <= 2;
    if (ok && ndeg === 2) {
      var D = rsub(rmul(rest[1], rest[1]), rmul([4n, 1n], rmul(rest[2], rest[0])));
      var dn = bsqrt(D[0]), dd = bsqrt(D[1]);
      if (dn === null || dd === null) ok = false;
      else { var two = rmul([2n, 1n], rest[2]), mb = rsub([0n, 1n], rest[1]); others = [rdiv(radd(mb, [dn, dd]), two), rdiv(rsub(mb, [dn, dd]), two)]; }
    } else if (ok && ndeg === 1) others = [rdiv(rsub([0n, 1n], rest[0]), rest[1])];
    var keys = others.map(function (r) { return ratStr(r[0], r[1]); });
    if (ok) {                                // rebuild the right-hand side and compare with q^2 - f
      var rhs = [lead]; roots.concat(others).forEach(function (r) { rhs = pmul(rhs, [rsub([0n, 1n], r), [1n, 1n]]); });
      rhs = ptrim(rhs); ok = rhs.length === lhs.length && rhs.every(function (c, i) { return c[0] === lhs[i][0] && c[1] === lhs[i][1]; });
    }
    var infNew = 6 - (lhs.length - 1) - parents.filter(function (pp) { return pp.key === '∞'; }).length;   // new points at infinity
    for (k = 0; k < infNew; k++) keys.push('∞');
    return { key: me.key, round: me.born, parents: parents, cubic: q, lead: lead, others: keys, ok: ok && keys.indexOf(me.key) >= 0 };
  };
  /* the whole run as a certificate that can be checked without searching: the curve, the starting values, and for every
     other point one relation (its four parents with their y, the cubic y = q(x) through them with coefficients q0..q3,
     and the two further intersection points), in the order of discovery, so that every parent is a starting value or
     the point of an earlier relation.  check_witnesses.py verifies such a file. */
  State.prototype.certificate = function () {
    var self = this, S = function (q) { return ratStr(q[0], q[1]); }, W = [];
    this.pts.forEach(function (pt, n) {
      var x = self.explain(n); if (!x) return;
      W.push({ x: x.key, parents: x.parents.map(function (pp) { return [pp.key, pp.y ? S(pp.y) : null]; }), cubic: x.cubic.map(S), others: x.others });
    });
    return { format: 'four-point closure certificate, version 1', curve: this.given.map(String),
             start: this.pts.filter(function (pt) { return !pt.from; }).map(function (pt) { return pt.key; }), witnesses: W };
  };
  function certText(st) {                  // the certificate as JSON text, one relation per line
    var c = st.certificate();
    return '{"format": ' + JSON.stringify(c.format) + ',\n "curve": ' + JSON.stringify(c.curve) + ',\n "start": ' + JSON.stringify(c.start) + ',\n "witnesses": [\n' +
      c.witnesses.map(function (w) { return JSON.stringify(w); }).join(',\n') + '\n]}\n';
  }
  function polyStr(q, name) {              // "a x^3 + b x^2 + c x + d" with rational coefficients (ascending input)
    var s = '', k;
    for (k = q.length - 1; k >= 0; k--) {
      var c = q[k]; if (c[0] === 0n) continue;
      var neg = c[0] < 0n, mag = ratStr(neg ? -c[0] : c[0], c[1]);
      s += (s ? (neg ? ' − ' : ' + ') : (neg ? '−' : '')) + ((mag === '1' && k > 0) ? '' : mag + (k > 0 ? '·' : '')) + (k > 0 ? name + (k > 1 ? '^' + k : '') : '');
    }
    return s || '0';
  }
  function start(coefText, seedText, opts) {
    var co = coefText.trim().split(/[\s,;]+/).filter(Boolean).map(function (s) { return BigInt(s.replace(/−/g, '-')); });
    if (co.length !== 7 || (co[0] === 0n && co[1] === 0n)) throw new Error('give the seven coefficients of f, from x^6 down to the constant term (f of degree 5 or 6)');
    var st = new State(co, opts), seeds = seedText.split(/[\s,;]+/).filter(Boolean);
    st.seedCount = 0;
    seeds.forEach(function (s) { var r = /^(∞|inf|infinity|oo)$/i.test(s.trim()) ? [1n, 0n] : parseRat(s); st.seedCount += st.addX(r[0], r[1], true); });
    st.addX(1n, 0n, false);                // rational points at infinity, when there are any, are known from the start
    return st;
  }
  function result(st) {
    var xs = st.pts.slice().sort(function (u, v) { return u.h < v.h ? -1 : u.h > v.h ? 1 : u.xa < v.xa ? -1 : u.xa > v.xa ? 1 : (u.key === '∞' ? 1 : -1); });
    var npts = 0; xs.forEach(function (p) { npts += p.weier ? 1 : 2; });
    return { xs: xs.map(function (p) { return p.key; }), points: npts, infinity: st.index.has('∞'), rounds: st.round, tests: st.tests, exact: st.exact, relations: st.relations };
  }
  return { start: start, result: result, certText: certText, polyStr: polyStr, ratStr: ratStr, PRIMES: PRIMES };
})();

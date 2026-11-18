/*
 * cubic - find new rational points on y^2 = f(x) (deg f = 6) from known ones:
 * for 4 known points (x_i, y_i) the cubic p through them satisfies
 * p^2 - f = (p3^2 - f6) prod (x - x_i) * q(x) with q quadratic; rational roots
 * of q give new points (the group-law relation P5 + P6 ~ 3(inf+ + inf-) - sum P_i).
 * A sieve on the classes of the points in the Jacobians J(F_p) for small primes p ("index sieve", see below) leaves the
 * 4-subsets that can have rational residual points; these are verified exactly (GMP), and the passes are iterated to
 * closure.  When the index sieve cannot be used (-K, a quintic, too few usable primes; or CUBIC_OLD set) the earlier
 * prefilter is used: the discriminant of q modulo several 16-bit primes with quadratic-residue tables, then modulo four
 * 31-bit primes.  Both filters only discard 4-subsets that have no rational residual point, so the closure is the same.
 *
 * usage: cubic 'f0 .. f6' points.txt [-t threads] [-maxpass k] [-notangent] [-noinf] [-six n] [-K k]
 *   points.txt: lines "(a : b : c)" (as printed by g2rp) or "a b c"
 *   -K k: restricted closure — the four-point pass tests only the 4-subsets whose two lowest indices are among the k
 *   lowest-height seed points (the seed is sorted by max(|a|, |c|) when -K is given; k is reduced to the number of seed
 *   x-values when the seed is smaller, so points found later never enter those two slots), so that pass costs
 *   O(k^2 n^2) subsets instead of O(n^4).  The tangent (Hermite) and infinity passes
 *   are NOT restricted: they are O(n^3) at most and run in full (use -notangent -noinf for the four-point relations
 *   alone).  A heuristic, not a proof of completeness: K = 60 recovers the full known count of every catalogue curve,
 *   with the tangent/infinity passes and without them (tools/cubic_k_check.sh), K = 30 loses two x-values of the 614
 *   curve (strategy report, 2.3).
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <gmp.h>
#include <omp.h>
#include <immintrin.h>
#include <sys/mman.h>

typedef unsigned __int128 u128;
typedef __int128 i128;

#define NP 8                 /* table primes */
static uint64_t PR[NP] = {65521, 65519, 65497, 65479, 65449, 65447, 65437, 65423};   /* primes just below 2^16 */
static uint8_t *QRB[NP];   /* byte tables of quadratic residuosity */
static uint8_t *QR[NP];      /* quadratic residue bit tables */
#define NBIG 4
static uint64_t BIGP[NBIG] = {2147483647ULL, 2147483629ULL, 2147483587ULL, 2147483579ULL};  /* primes < 2^31 */
static double PINV[NP + NBIG];

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (a * b) % p; }   /* a, b < p < 2^32 */
/* fast reduction of x < 2^62 modulo p (< 2^31) using a double reciprocal */
static inline uint64_t redp(uint64_t x, uint64_t p, double pinv) {
    uint64_t q = (uint64_t)((double)x * pinv);
    int64_t r = (int64_t)(x - q * p);
    while (r < 0) r += p;
    while (r >= (int64_t)p) r -= p;
    return (uint64_t)r;
}
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t p) { uint64_t r = 1; a %= p; while (e) { if (e & 1) r = mulmod(r, a, p); a = mulmod(a, a, p); e >>= 1; } return r; }
static uint64_t invmod(uint64_t a, uint64_t p) { return powmod(a, p - 2, p); }
static inline int legendre_big(uint64_t a, uint64_t p) { if (a == 0) return 1; return powmod(a, (p-1)/2, p) == 1; }

typedef struct { mpz_t a, b, c; } pt_t;          /* x = a/c, y = b/c^3 */
static pt_t *pts; static int npts, cappts;
static uint64_t *xr[NP+NBIG], *yr[NP+NBIG];            /* residues of x, y per prime */
static uint64_t *invd[NP+NBIG];   /* invd[k][i*ncap+j] = (x_i - x_j)^-1 mod p (0 if equal) */
static uint64_t *inv2y[NP+NBIG];  /* (2 y_i)^-1 mod p (0 if y_i = 0) */
static uint64_t *yd1[NP+NBIG];    /* y_i' = f'(x_i) / (2 y_i) mod p (meaningless if inv2y = 0) */
static uint64_t *sx3[NP+NBIG];    /* s x_i^3 mod p, s = sqrt(f6) (0 if none) */
static uint8_t *badr[NP+NBIG];    /* 1 if the point's denominator c is 0 mod p: its residues are meaningless, every modular test involving it is undecided for that prime (Codex round 3, PR #15) */
static mpz_t sqrt_f6z; static int has_sqrt_f6 = 0;   /* s with s^2 = f6 (kept in GMP: a leading coefficient >= 2^126 has a root beyond int64 -- Codex round 5, PR #15); has_sqrt_f6 = 0 if f6 = 0 or not a square (then no inf relations) */
static int f6_is_zero = 0;
static int ncap = 0;
static uint64_t fr[NP+NBIG][7];   /* f mod p */

static void add_pt(mpz_t a, mpz_t b, mpz_t c) {
    if (npts == cappts) {
        int newcap = cappts ? 2*cappts : 1024;
        pts = realloc(pts, newcap * sizeof(pt_t));
        for (int k = 0; k < NP+NBIG; k++) {
            xr[k] = realloc(xr[k], newcap * 8); yr[k] = realloc(yr[k], newcap * 8); inv2y[k] = realloc(inv2y[k], newcap * 8); yd1[k] = realloc(yd1[k], newcap * 8); sx3[k] = realloc(sx3[k], newcap * 8); badr[k] = realloc(badr[k], newcap);
            uint64_t *nd = calloc((size_t)newcap * newcap, 8);
            if (invd[k]) { for (int i = 0; i < npts; i++) memcpy(nd + (size_t)i * newcap, invd[k] + (size_t)i * cappts, (size_t)npts * 8); free(invd[k]); }
            invd[k] = nd;
        }
        cappts = newcap; ncap = newcap;
    }
    mpz_init_set(pts[npts].a, a); mpz_init_set(pts[npts].b, b); mpz_init_set(pts[npts].c, c);
    for (int k = 0; k < NP+NBIG; k++) {
        uint64_t p = k < NP ? PR[k] : BIGP[k-NP];
        uint64_t am = mpz_fdiv_ui(a, p), bm = mpz_fdiv_ui(b, p), cm = mpz_fdiv_ui(c, p);
        uint64_t ci = invmod(cm, p);
        uint64_t x = mulmod(am, ci, p), y = mulmod(bm, mulmod(ci, mulmod(ci, ci, p), p), p);
        badr[k][npts] = (cm == 0);   /* x = a/c has no residue mod p: the point is excluded from every test mod p (inverse tables below are 0 for it) */
        xr[k][npts] = x; yr[k][npts] = y;
        inv2y[k][npts] = y ? invmod(mulmod(2, y, p), p) : 0;
        { uint64_t f1 = 0, pw = 1; for (int d = 1; d <= 6; d++) { f1 = (f1 + mulmod(mulmod(d, fr[k][d], p), pw, p)) % p; pw = mulmod(pw, x, p); } yd1[k][npts] = mulmod(f1, inv2y[k][npts], p); }
        { uint64_t sk = has_sqrt_f6 ? mpz_fdiv_ui(sqrt_f6z, p) : 0; sx3[k][npts] = mulmod(sk, mulmod(x, mulmod(x, x, p), p), p); }
        for (int i = 0; i < npts; i++) {
            uint64_t d = (x + p - xr[k][i]) % p;
            uint64_t iv = (d && !badr[k][npts] && !badr[k][i]) ? invmod(d, p) : 0;
            invd[k][(size_t)npts * ncap + i] = iv;
            invd[k][(size_t)i * ncap + npts] = iv ? (p - iv) % p : 0;
        }
        invd[k][(size_t)npts * ncap + npts] = 0;
    }
    npts++;
}

static mpz_t f[7];
static int is_weier(int i);

/* per-subset, per-prime precomputation: Lagrange basis L_i (cubic) and quartic Q = prod (x - x_i) */
typedef struct { uint64_t L[4][4]; uint64_t Q[5]; uint64_t yL[4][4]; } sub_t;

static int prep_subset(int k, const int *idx, sub_t *S) {
    uint64_t p = k < NP ? PR[k] : BIGP[k-NP]; double pinv = PINV[k];
    uint64_t x[4], y[4];
    for (int i = 0; i < 4; i++) { if (badr[k][idx[i]]) return 0; x[i] = xr[k][idx[i]]; y[i] = yr[k][idx[i]]; }
    uint64_t Q[5] = {1, 0, 0, 0, 0}; int dq = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t nx = x[i] ? p - x[i] : 0;
        for (int d = dq + 1; d >= 1; d--) Q[d] = redp(Q[d-1] + Q[d] * nx, p, pinv);
        Q[0] = redp(Q[0] * nx, p, pinv); dq++;
    }
    memcpy(S->Q, Q, sizeof Q);
    for (int i = 0; i < 4; i++) {
        uint64_t di = 1;
        for (int j = 0; j < 4; j++) if (j != i) { uint64_t iv = invd[k][(size_t)idx[i] * ncap + idx[j]]; if (!iv) return 0; di = redp(di * iv, p, pinv); }
        uint64_t c[4]; uint64_t carry = Q[4];
        c[3] = carry;
        for (int d = 3; d >= 1; d--) { carry = redp(Q[d] + carry * x[i], p, pinv); c[d-1] = carry; }
        for (int d = 0; d < 4; d++) { S->L[i][d] = redp(c[d] * di, p, pinv); S->yL[i][d] = redp(S->L[i][d] * y[i], p, pinv); }
    }
    return 1;
}

/* for sign pattern eps (bit i = sign of y_i), compute discriminant of q; returns 1 if it is a QR mod p (k-th prime) */
static int test_pattern(int k, const sub_t *S, int eps, uint64_t *disc_out) {
    uint64_t p = k < NP ? PR[k] : BIGP[k-NP]; double pinv = PINV[k];
    uint64_t P[4] = {0,0,0,0};
    for (int i = 0; i < 4; i++) {
        if (eps & (1 << i)) for (int d = 0; d < 4; d++) { P[d] += p - S->yL[i][d]; }
        else for (int d = 0; d < 4; d++) { P[d] += S->yL[i][d]; }
    }
    for (int d = 0; d < 4; d++) P[d] %= p;
    uint64_t r[7] = {0,0,0,0,0,0,0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r[i+j] += redp(P[i] * P[j], p, pinv);
    for (int d = 0; d < 7; d++) { r[d] = (r[d] + 7 * p - fr[k][d]) % p; }
    uint64_t q[3];
    for (int d = 6; d >= 4; d--) {
        uint64_t t = r[d]; q[d-4] = t;
        if (t) for (int e = 0; e <= 4; e++) { r[d-4+e] = (r[d-4+e] + p - redp(t * S->Q[e], p, pinv)) % p; }
    }
    if (q[2] == 0) { *disc_out = 0; return 1; }
    uint64_t disc = (redp(q[1] * q[1], p, pinv) + p - redp(4 * redp(q[2] * q[0], p, pinv), p, pinv)) % p;
    *disc_out = disc;
    if (k < NP) return (QR[k][disc >> 3] >> (disc & 7)) & 1;
    return legendre_big(disc, p);
}

/* exact verification: returns number of new points appended to out (0..2) */
typedef struct { mpz_t a, b, c; int idx[4]; int eps; int ys; } newpt_t;   /* ys (set by the four-point exact check and by inf_exact): the sign of y at the residual point, relative to the stored b >= 0 (0 if b = 0) */
static int exact_check_mpq(const int *idx, int eps, newpt_t *out) {
    mpq_t x[4], y[4], L[4][4], P[4], t, u;
    for (int i = 0; i < 4; i++) { mpq_init(x[i]); mpq_init(y[i]); mpq_init(P[i]); for (int d = 0; d < 4; d++) mpq_init(L[i][d]); }
    mpq_init(t); mpq_init(u);
    for (int i = 0; i < 4; i++) {
        mpq_set_num(x[i], pts[idx[i]].a); mpq_set_den(x[i], pts[idx[i]].c); mpq_canonicalize(x[i]);
        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, pts[idx[i]].c, 3);
        mpq_set_num(y[i], pts[idx[i]].b); mpq_set_den(y[i], c3); mpq_canonicalize(y[i]); mpz_clear(c3);
        if (eps & (1 << i)) mpq_neg(y[i], y[i]);
    }
    /* Lagrange basis */
    for (int i = 0; i < 4; i++) {
        mpq_t poly[4]; for (int d = 0; d < 4; d++) mpq_init(poly[d]);
        mpq_set_ui(poly[0], 1, 1); int deg = 0;
        mpq_t den; mpq_init(den); mpq_set_ui(den, 1, 1);
        for (int j = 0; j < 4; j++) if (j != i) {
            for (int d = deg + 1; d >= 1; d--) { mpq_mul(t, poly[d], x[j]); mpq_sub(poly[d], poly[d-1], t); }
            mpq_mul(poly[0], poly[0], x[j]); mpq_neg(poly[0], poly[0]); deg++;
            mpq_sub(t, x[i], x[j]); mpq_mul(den, den, t);
        }
        for (int d = 0; d < 4; d++) { mpq_div(L[i][d], poly[d], den); mpq_mul(L[i][d], L[i][d], y[i]); mpq_clear(poly[d]); }
        mpq_clear(den);
    }
    for (int d = 0; d < 4; d++) { mpq_set_ui(P[d], 0, 1); for (int i = 0; i < 4; i++) mpq_add(P[d], P[d], L[i][d]); }
    /* r = P^2 - f */
    mpq_t r[7]; for (int d = 0; d < 7; d++) mpq_init(r[d]);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { mpq_mul(t, P[i], P[j]); mpq_add(r[i+j], r[i+j], t); }
    for (int d = 0; d < 7; d++) { mpq_set_z(t, f[d]); mpq_sub(r[d], r[d], t); }
    /* divide by prod (x - x_i) */
    for (int i = 0; i < 4; i++) { /* synthetic division by (x - x_i), from the top */
        int deg = 6;
        for (int d = deg; d >= 1; d--) { mpq_mul(t, r[d], x[i]); mpq_add(r[d-1], r[d-1], t); }
        /* now r[0] is the remainder (should be 0); shift down */
        for (int d = 0; d < 6; d++) mpq_set(r[d], r[d+1]);
        mpq_set_ui(r[6], 0, 1);
    }
    /* quotient is r[0..2]: q2 x^2 + q1 x + q0 */
    int nnew = 0;
    if (mpq_sgn(r[2]) != 0) {
        mpq_t disc; mpq_init(disc);
        mpq_mul(disc, r[1], r[1]); mpq_mul(t, r[2], r[0]); mpq_set_ui(u, 4, 1); mpq_mul(t, t, u); mpq_sub(disc, disc, t);
        if (mpq_sgn(disc) >= 0 && mpz_perfect_square_p(mpq_numref(disc)) && mpz_perfect_square_p(mpq_denref(disc))) {
            mpq_t sq; mpq_init(sq); mpz_sqrt(mpq_numref(sq), mpq_numref(disc)); mpz_sqrt(mpq_denref(sq), mpq_denref(disc)); mpq_canonicalize(sq);
            for (int s = 0; s < 2; s++) {
                mpq_t xn, yn; mpq_init(xn); mpq_init(yn);
                if (s) mpq_sub(xn, r[1], sq); else mpq_add(xn, r[1], sq);
                mpq_neg(xn, xn); mpq_set_ui(u, 2, 1); mpq_mul(t, r[2], u); mpq_div(xn, xn, t);
                /* y = P(xn) */
                mpq_set(yn, P[3]); for (int d = 2; d >= 0; d--) { mpq_mul(yn, yn, xn); mpq_add(yn, yn, P[d]); }
                /* verify y^2 = f(x) */
                mpq_t fx; mpq_init(fx); mpq_set_z(fx, f[6]); for (int d = 5; d >= 0; d--) { mpq_mul(fx, fx, xn); mpq_set_z(t, f[d]); mpq_add(fx, fx, t); }
                mpq_mul(t, yn, yn);
                if (mpq_equal(t, fx)) {
                    /* a = num(x), c = den(x), b = y * c^3 */
                    mpz_init_set(out[nnew].a, mpq_numref(xn)); mpz_init_set(out[nnew].c, mpq_denref(xn));
                    mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, mpq_denref(xn), 3);
                    mpq_set_z(u, c3); mpq_mul(u, yn, u);
                    mpz_init_set(out[nnew].b, mpq_numref(u));   /* denominator must be 1 */
                    if (mpz_cmp_ui(mpq_denref(u), 1) != 0) { fprintf(stderr, "warning: non-integral b\n"); }
                    if (mpz_sgn(out[nnew].b) < 0) mpz_neg(out[nnew].b, out[nnew].b);
                    memcpy(out[nnew].idx, idx, sizeof(int)*4); out[nnew].eps = eps;
                    mpz_clear(c3);
                    nnew++;
                }
                mpq_clear(fx); mpq_clear(xn); mpq_clear(yn);
            }
            mpq_clear(sq);
        }
        mpq_clear(disc);
    }
    for (int d = 0; d < 7; d++) mpq_clear(r[d]);
    for (int i = 0; i < 4; i++) { mpq_clear(x[i]); mpq_clear(y[i]); mpq_clear(P[i]); for (int d = 0; d < 4; d++) mpq_clear(L[i][d]); }
    mpq_clear(t); mpq_clear(u);
    return nnew;
}

/* exact verification, fraction-free (mpz only).  With x_i = a_i/c_i, y_i = b_i/c_i^3 (gcd(a_i,c_i) = 1, enforced when
 * the points are read), D_ij = a_i c_j - a_j c_i and N_i(x) = prod_{j != i} (c_j x - a_j):
 *     y_i L_i(x) = b_i N_i(x) / prod_{j != i} D_ij            (L_i = Lagrange basis; the c_i^3 cancel),
 * so with W = prod_{i<j} D_ij and W_i = W / prod_{j != i} D_ij (integers, sign included)
 *     G(x) := W p(x) = sum_i s_i b_i W_i N_i(x)             (integer cubic).
 * R = G^2 - W^2 f is divisible in Z[x] by the primitive quartic Q_Z = prod_i (c_i x - a_i) (Gauss's lemma); the
 * integer quotient q equals (W^2/C) times the residual quadratic (C = prod c_i), so it has the same roots and
 * disc(q) = (W^2/C)^2 disc(residual): the rational-square tests agree with the mpq version.  Only the three top
 * coefficients of R are needed for the quotient.  Returns the number of residual rational points appended (0..2). */
static int check_exact = 0;   /* CUBIC_CHECKEXACT: cross-validate the fraction-free exact check against the mpq version */
static int exact_check_ff(const int *idx, int eps, newpt_t *out) {
    static __thread int init = 0;
    static __thread mpz_t D[4][4], W, Wi, bw[4], N[4][4], G[4], R[3], Qz[5], q[3], disc, sq, t, u, num, den, F, cp;
    if (!init) {
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { mpz_init(D[i][j]); mpz_init(N[i][j]); }
        mpz_init(W); mpz_init(Wi); mpz_init(disc); mpz_init(sq); mpz_init(t); mpz_init(u); mpz_init(num); mpz_init(den); mpz_init(F); mpz_init(cp);
        for (int i = 0; i < 4; i++) { mpz_init(bw[i]); mpz_init(G[i]); }
        for (int i = 0; i < 3; i++) { mpz_init(R[i]); mpz_init(q[i]); }
        for (int i = 0; i < 5; i++) mpz_init(Qz[i]);
        init = 1;
    }
    const pt_t *P[4]; for (int i = 0; i < 4; i++) P[i] = &pts[idx[i]];
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) {
        mpz_mul(D[i][j], P[i]->a, P[j]->c); mpz_submul(D[i][j], P[j]->a, P[i]->c);
        if (mpz_sgn(D[i][j]) == 0) return 0;   /* equal x-values: not a valid 4-subset */
    }
    mpz_mul(W, D[0][1], D[0][2]); mpz_mul(W, W, D[0][3]); mpz_mul(W, W, D[1][2]); mpz_mul(W, W, D[1][3]); mpz_mul(W, W, D[2][3]);
    for (int i = 0; i < 4; i++) {
        mpz_set_ui(t, 1); int neg = 0;
        for (int j = 0; j < 4; j++) if (j != i) { if (j > i) mpz_mul(t, t, D[i][j]); else { mpz_mul(t, t, D[j][i]); neg ^= 1; } }
        if (neg) mpz_neg(t, t);
        mpz_divexact(Wi, W, t);
        mpz_mul(bw[i], P[i]->b, Wi); if (eps & (1 << i)) mpz_neg(bw[i], bw[i]);
        /* N_i = prod_{j != i} (c_j x - a_j) */
        mpz_set_ui(N[i][0], 1); int deg = 0;
        for (int j = 0; j < 4; j++) if (j != i) {
            mpz_mul(N[i][deg+1], P[j]->c, N[i][deg]);
            for (int d = deg; d >= 1; d--) { mpz_mul(N[i][d], N[i][d], P[j]->a); mpz_neg(N[i][d], N[i][d]); mpz_addmul(N[i][d], P[j]->c, N[i][d-1]); }
            mpz_mul(N[i][0], N[i][0], P[j]->a); mpz_neg(N[i][0], N[i][0]);
            deg++;
        }
    }
    for (int d = 0; d < 4; d++) { mpz_set_ui(G[d], 0); for (int i = 0; i < 4; i++) mpz_addmul(G[d], bw[i], N[i][d]); }
    /* R[0..2] = coefficients of x^4, x^5, x^6 of G^2 - W^2 f */
    mpz_mul(R[2], G[3], G[3]);
    mpz_mul(R[1], G[2], G[3]); mpz_mul_2exp(R[1], R[1], 1);
    mpz_mul(R[0], G[1], G[3]); mpz_mul_2exp(R[0], R[0], 1); mpz_addmul(R[0], G[2], G[2]);
    mpz_mul(t, W, W);
    mpz_submul(R[2], t, f[6]); mpz_submul(R[1], t, f[5]); mpz_submul(R[0], t, f[4]);
    /* Q_Z = prod_i (c_i x - a_i) */
    mpz_set_ui(Qz[0], 1); { int deg = 0; for (int i = 0; i < 4; i++) {
        mpz_mul(Qz[deg+1], P[i]->c, Qz[deg]);
        for (int d = deg; d >= 1; d--) { mpz_mul(Qz[d], Qz[d], P[i]->a); mpz_neg(Qz[d], Qz[d]); mpz_addmul(Qz[d], P[i]->c, Qz[d-1]); }
        mpz_mul(Qz[0], Qz[0], P[i]->a); mpz_neg(Qz[0], Qz[0]); deg++; } }
    /* long division by Q_Z, top three steps */
    mpz_divexact(q[2], R[2], Qz[4]);
    mpz_submul(R[1], q[2], Qz[3]); mpz_submul(R[0], q[2], Qz[2]);
    mpz_divexact(q[1], R[1], Qz[4]);
    mpz_submul(R[0], q[1], Qz[3]);
    mpz_divexact(q[0], R[0], Qz[4]);
    if (mpz_sgn(q[2]) == 0) return 0;
    mpz_mul(disc, q[1], q[1]); mpz_mul(t, q[2], q[0]); mpz_mul_2exp(t, t, 2); mpz_sub(disc, disc, t);
    if (mpz_sgn(disc) < 0 || !mpz_perfect_square_p(disc)) return 0;
    mpz_sqrt(sq, disc);
    int nnew = 0;
    for (int s = 0; s < 2; s++) {
        if (s) mpz_sub(num, q[1], sq); else mpz_add(num, q[1], sq);
        mpz_neg(num, num); mpz_mul_2exp(den, q[2], 1);
        mpz_gcd(t, num, den); mpz_divexact(num, num, t); mpz_divexact(den, den, t);
        if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }
        /* F = c^6 f(a/c) = sum_d f_d a^d c^(6-d) (Horner), b = sqrt(F) */
        mpz_set(F, f[6]); mpz_set_ui(cp, 1);
        for (int d = 5; d >= 0; d--) { mpz_mul(cp, cp, den); mpz_mul(F, F, num); mpz_addmul(F, f[d], cp); }
        if (mpz_sgn(F) < 0 || !mpz_perfect_square_p(F)) { fprintf(stderr, "warning: residual root is not on the curve\n"); continue; }
        mpz_init_set(out[nnew].a, num); mpz_init_set(out[nnew].c, den); mpz_init(out[nnew].b); mpz_sqrt(out[nnew].b, F);
        memcpy(out[nnew].idx, idx, sizeof(int)*4); out[nnew].eps = eps;
        /* the cubic is G / W, so y at the root is (G_3 a^3 + G_2 a^2 c + G_1 a c^2 + G_0 c^3) / (W c^3), with c > 0 */
        mpz_set(u, G[3]); mpz_set_ui(cp, 1);
        for (int d = 2; d >= 0; d--) { mpz_mul(cp, cp, den); mpz_mul(u, u, num); mpz_addmul(u, G[d], cp); }
        out[nnew].ys = mpz_sgn(u) * mpz_sgn(W);
        if (check_exact) { mpz_mul(t, W, out[nnew].b); if (mpz_cmpabs(u, t) != 0) { fprintf(stderr, "CUBIC_CHECKEXACT: the cubic does not take the value +-y at a residual root (subset %d %d %d %d eps %d)\n", idx[0], idx[1], idx[2], idx[3], eps); abort(); } }
        nnew++;
    }
    return nnew;
}
static int exact_check(const int *idx, int eps, newpt_t *out) {
    int n = exact_check_ff(idx, eps, out);
    if (check_exact) {
        newpt_t o2[2]; int n2 = exact_check_mpq(idx, eps, o2);
        int bad = (n != n2);
        for (int i = 0; i < n && !bad; i++) bad = mpz_cmp(out[i].a, o2[i].a) || mpz_cmp(out[i].b, o2[i].b) || mpz_cmp(out[i].c, o2[i].c);
        if (bad) { fprintf(stderr, "CUBIC_CHECKEXACT: mismatch for subset %d %d %d %d eps %d (%d vs %d points)\n", idx[0], idx[1], idx[2], idx[3], eps, n, n2); abort(); }
        for (int i = 0; i < n2; i++) { mpz_clear(o2[i].a); mpz_clear(o2[i].b); mpz_clear(o2[i].c); }
    }
    return n;
}

static int known_x(mpz_t a, mpz_t c) {   /* linear scan (fine for a few hundred points) */
    for (int i = 0; i < npts; i++) if (mpz_cmp(pts[i].a, a) == 0 && mpz_cmp(pts[i].c, c) == 0) return 1;
    return 0;
}


/* ---------------- Hermite (tangent) relations ----------------
 * conditions: p(x_i) = y_i, and for multiplicity m_i >= 2: p'(x_i) = y_i^(1) = f'(x_i)/(2 y_i),
 * m_i >= 3: p''(x_i) = y_i^(2), m_i = 4: p'''(x_i) = y_i^(3)  (derivatives of y = sqrt(f)).
 * modular version: solve the 4x4 system by Gaussian elimination mod p; residual quadratic
 * from dividing p^2 - f by prod (x - x_i)^{m_i}.
 */
static void y_derivs_mod2(uint64_t p, const uint64_t *fm, uint64_t x, uint64_t y, uint64_t *d, uint64_t inv2yv) {
    uint64_t pw[7]; pw[0] = 1; for (int i = 1; i < 7; i++) pw[i] = mulmod(pw[i-1], x, p);
    uint64_t f1 = 0, f2 = 0, f3 = 0;
    for (int i = 1; i <= 6; i++) f1 = (f1 + mulmod(mulmod(i, fm[i], p), pw[i-1], p)) % p;
    for (int i = 2; i <= 6; i++) f2 = (f2 + mulmod(mulmod((uint64_t)i*(i-1) % p, fm[i], p), pw[i-2], p)) % p;
    for (int i = 3; i <= 6; i++) f3 = (f3 + mulmod(mulmod((uint64_t)i*(i-1)*(i-2) % p, fm[i], p), pw[i-3], p)) % p;
    d[0] = y;
    d[1] = mulmod(f1, inv2yv, p);
    d[2] = mulmod((f2 + p - mulmod(2, mulmod(d[1], d[1], p), p)) % p, inv2yv, p);
    d[3] = mulmod((f3 + p - mulmod(6, mulmod(d[1], d[2], p), p)) % p, inv2yv, p);
}
static void y_derivs_mod(uint64_t p, const uint64_t *fm, uint64_t x, uint64_t y, uint64_t *d) {
    uint64_t pw[7]; pw[0] = 1; for (int i = 1; i < 7; i++) pw[i] = mulmod(pw[i-1], x, p);
    uint64_t f1 = 0, f2 = 0, f3 = 0;
    for (int i = 1; i <= 6; i++) f1 = (f1 + mulmod(mulmod(i, fm[i], p), pw[i-1], p)) % p;
    for (int i = 2; i <= 6; i++) f2 = (f2 + mulmod(mulmod((uint64_t)i*(i-1) % p, fm[i], p), pw[i-2], p)) % p;
    for (int i = 3; i <= 6; i++) f3 = (f3 + mulmod(mulmod((uint64_t)i*(i-1)*(i-2) % p, fm[i], p), pw[i-3], p)) % p;
    uint64_t inv2y = invmod(mulmod(2, y, p), p);
    d[0] = y;
    d[1] = mulmod(f1, inv2y, p);
    d[2] = mulmod((f2 + p - mulmod(2, mulmod(d[1], d[1], p), p)) % p, inv2y, p);
    d[3] = mulmod((f3 + p - mulmod(6, mulmod(d[1], d[2], p), p)) % p, inv2y, p);
}
/* fraction-free solve of the 4x4 system A c = rhs (A[i][0..3], rhs = A[i][4]) mod p:
   returns det != 0 and writes det*c into A[i][4] (Bareiss / Cramer via determinants), no inverses */
static uint64_t det4_mod(uint64_t p, double pinv, const uint64_t M[4][4]) {
    uint64_t d = 0;
    /* Laplace on first row with 3x3 minors */
    static const int m3[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for (int j = 0; j < 4; j++) {
        const int *c = m3[j];
        uint64_t a = M[1][c[0]], b = M[1][c[1]], cc = M[1][c[2]];
        uint64_t dd = M[2][c[0]], e = M[2][c[1]], ff = M[2][c[2]];
        uint64_t g = M[3][c[0]], h = M[3][c[1]], i = M[3][c[2]];
        uint64_t t1 = redp(e * i, p, pinv), t2 = redp(ff * h, p, pinv), t3 = redp(dd * i, p, pinv), t4 = redp(ff * g, p, pinv), t5 = redp(dd * h, p, pinv), t6 = redp(e * g, p, pinv);
        uint64_t minor = (redp(a * ((t1 + p - t2) % p), p, pinv) + p - redp(b * ((t3 + p - t4) % p), p, pinv) + redp(cc * ((t5 + p - t6) % p), p, pinv)) % p;
        uint64_t term = redp(M[0][j] * minor, p, pinv);
        d = (j & 1) ? (d + p - term) % p : (d + term) % p;
    }
    return d;
}
static int solve4_mod(uint64_t p, uint64_t A[4][5]) {
    double pinv = 1.0 / (double)p;
    uint64_t M[4][4]; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) M[i][j] = A[i][j];
    uint64_t det = det4_mod(p, pinv, M);
    if (det == 0) return 0;
    uint64_t sol[4];
    for (int c = 0; c < 4; c++) {   /* Cramer: replace column c by rhs */
        uint64_t N[4][4]; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) N[i][j] = (j == c) ? A[i][4] : A[i][j];
        sol[c] = det4_mod(p, pinv, N);   /* = det * c_c */
    }
    for (int c = 0; c < 4; c++) A[c][4] = sol[c];
    A[0][0] = det;   /* pass det back in A[0][0] */
    return 1;
}
/* returns 1 = candidate, 0 = rejected, -1 = undecided (y_i = 0 mod p or singular system mod p) */
static int hermite_test_mod(int kp, const int *idx, const int *m, int k, int eps, uint64_t *disc_out) {
    uint64_t p = kp < NP ? PR[kp] : BIGP[kp-NP];
    for (int i = 0; i < k; i++) if (badr[kp][idx[i]]) return -1;   /* a point without residues mod p: undecided */
    uint64_t A[4][5]; int row = 0;
    uint64_t xs[4];
    for (int i = 0; i < k; i++) {
        uint64_t x = xr[kp][idx[i]], y = yr[kp][idx[i]];
        if (eps & (1 << i)) y = (p - y) % p;
        if (m[i] >= 2 && y == 0) return -1;   /* only tangent conditions need 1/(2y); a simple point with y = 0 (Weierstrass point) is fine */
        uint64_t d[4]; y_derivs_mod2(p, fr[kp], x, y, d, (eps & (1 << i)) ? (p - inv2y[kp][idx[i]]) % p : inv2y[kp][idx[i]]);
        xs[i] = x;
        uint64_t pw[4]; pw[0] = 1; for (int j = 1; j < 4; j++) pw[j] = mulmod(pw[j-1], x, p);
        for (int t = 0; t < m[i]; t++) {
            for (int j = 0; j < 4; j++) {
                if (j < t) { A[row][j] = 0; continue; }
                uint64_t fac = 1; for (int u = 0; u < t; u++) fac *= (j - u);
                A[row][j] = mulmod(fac % p, pw[j-t], p);
            }
            A[row][4] = d[t]; row++;
        }
    }
    if (!solve4_mod(p, A)) return -1;
    uint64_t P[4]; for (int j = 0; j < 4; j++) P[j] = A[j][4];
    uint64_t det = A[0][0], det2 = mulmod(det, det, p);
    uint64_t r[7] = {0,0,0,0,0,0,0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { r[i+j] += mulmod(P[i], P[j], p); if (r[i+j] >= p) r[i+j] -= p; }
    for (int d = 0; d < 7; d++) r[d] = (r[d] + p - mulmod(det2, fr[kp][d], p)) % p;
    int deg = 6;
    for (int i = 0; i < k; i++) for (int t = 0; t < m[i]; t++) {
        for (int d = deg; d >= 1; d--) r[d-1] = (r[d-1] + mulmod(r[d], xs[i], p)) % p;
        for (int d = 0; d < deg; d++) r[d] = r[d+1];
        r[deg] = 0; deg--;
    }
    if (r[2] == 0) { *disc_out = 0; return 1; }
    uint64_t disc = (mulmod(r[1], r[1], p) + p - mulmod(4 % p, mulmod(r[2], r[0], p), p)) % p;
    *disc_out = disc;
    if (kp < NP) return (QR[kp][disc >> 3] >> (disc & 7)) & 1;
    return legendre_big(disc, p);
}
static int hermite_exact(const int *idx, const int *m, int k, int eps, newpt_t *out) {
    mpq_t A[4][5], t, u, x[4], y[4];
    for (int r = 0; r < 4; r++) for (int c = 0; c < 5; c++) mpq_init(A[r][c]);
    mpq_init(t); mpq_init(u);
    int row = 0;
    for (int i = 0; i < k; i++) {
        mpq_init(x[i]); mpq_init(y[i]);
        mpq_set_num(x[i], pts[idx[i]].a); mpq_set_den(x[i], pts[idx[i]].c); mpq_canonicalize(x[i]);
        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, pts[idx[i]].c, 3);
        mpq_set_num(y[i], pts[idx[i]].b); mpq_set_den(y[i], c3); mpq_canonicalize(y[i]); mpz_clear(c3);
        if (eps & (1 << i)) mpq_neg(y[i], y[i]);
        mpq_t d[4], f1, f2, f3, pw[7]; for (int j = 0; j < 4; j++) mpq_init(d[j]); mpq_init(f1); mpq_init(f2); mpq_init(f3);
        for (int j = 0; j < 7; j++) { mpq_init(pw[j]); if (j == 0) mpq_set_ui(pw[j], 1, 1); else mpq_mul(pw[j], pw[j-1], x[i]); }
        for (int j = 1; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j, 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-1]); mpq_add(f1, f1, t); }
        for (int j = 2; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j*(j-1), 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-2]); mpq_add(f2, f2, t); }
        for (int j = 3; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j*(j-1)*(j-2), 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-3]); mpq_add(f3, f3, t); }
        mpq_set(d[0], y[i]);
        if (m[i] >= 2) {   /* derivatives of y = sqrt(f) only for tangent conditions (y_i = 0 for a Weierstrass point used as a simple point) */
            mpq_set_ui(u, 2, 1); mpq_mul(u, u, y[i]);
            mpq_div(d[1], f1, u);
            mpq_mul(t, d[1], d[1]); mpq_set_ui(d[2], 2, 1); mpq_mul(t, t, d[2]); mpq_sub(t, f2, t); mpq_div(d[2], t, u);
            mpq_mul(t, d[1], d[2]); mpq_set_ui(d[3], 6, 1); mpq_mul(t, t, d[3]); mpq_sub(t, f3, t); mpq_div(d[3], t, u);
        }
        for (int tt = 0; tt < m[i]; tt++) {
            for (int j = 0; j < 4; j++) {
                if (j < tt) { mpq_set_ui(A[row][j], 0, 1); continue; }
                unsigned long fac = 1; for (int v = 0; v < tt; v++) fac *= (j - v);
                mpq_set_ui(A[row][j], fac, 1); mpq_mul(A[row][j], A[row][j], pw[j-tt]);
            }
            mpq_set(A[row][4], d[tt]); row++;
        }
        for (int j = 0; j < 4; j++) mpq_clear(d[j]); mpq_clear(f1); mpq_clear(f2); mpq_clear(f3); for (int j = 0; j < 7; j++) mpq_clear(pw[j]);
    }
    int ok = 1;
    for (int col = 0; col < 4 && ok; col++) {
        int piv = -1; for (int r = col; r < 4; r++) if (mpq_sgn(A[r][col])) { piv = r; break; }
        if (piv < 0) { ok = 0; break; }
        if (piv != col) for (int c = 0; c < 5; c++) mpq_swap(A[col][c], A[piv][c]);
        for (int c = 4; c >= col; c--) mpq_div(A[col][c], A[col][c], A[col][col]);
        for (int r = 0; r < 4; r++) if (r != col && mpq_sgn(A[r][col])) { for (int c = 4; c >= col; c--) { mpq_mul(t, A[r][col], A[col][c]); mpq_sub(A[r][c], A[r][c], t); } }
    }
    int nnew = 0;
    if (ok) {
        mpq_t P[4], r[7]; for (int j = 0; j < 4; j++) { mpq_init(P[j]); mpq_set(P[j], A[j][4]); } for (int d = 0; d < 7; d++) mpq_init(r[d]);
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { mpq_mul(t, P[i], P[j]); mpq_add(r[i+j], r[i+j], t); }
        for (int d = 0; d < 7; d++) { mpq_set_z(t, f[d]); mpq_sub(r[d], r[d], t); }
        int deg = 6;
        for (int i = 0; i < k; i++) for (int tt = 0; tt < m[i]; tt++) {
            for (int d = deg; d >= 1; d--) { mpq_mul(t, r[d], x[i]); mpq_add(r[d-1], r[d-1], t); }
            if (mpq_sgn(r[0]) != 0) ok = 0;
            for (int d = 0; d < deg; d++) mpq_set(r[d], r[d+1]);
            mpq_set_ui(r[deg], 0, 1); deg--;
        }
        if (ok && mpq_sgn(r[2]) != 0) {
            mpq_t disc; mpq_init(disc);
            mpq_mul(disc, r[1], r[1]); mpq_mul(t, r[2], r[0]); mpq_set_ui(u, 4, 1); mpq_mul(t, t, u); mpq_sub(disc, disc, t);
            if (mpq_sgn(disc) >= 0 && mpz_perfect_square_p(mpq_numref(disc)) && mpz_perfect_square_p(mpq_denref(disc))) {
                mpq_t sq; mpq_init(sq); mpz_sqrt(mpq_numref(sq), mpq_numref(disc)); mpz_sqrt(mpq_denref(sq), mpq_denref(disc)); mpq_canonicalize(sq);
                for (int s = 0; s < 2; s++) {
                    mpq_t xn, yn; mpq_init(xn); mpq_init(yn);
                    if (s) mpq_sub(xn, r[1], sq); else mpq_add(xn, r[1], sq);
                    mpq_neg(xn, xn); mpq_set_ui(u, 2, 1); mpq_mul(t, r[2], u); mpq_div(xn, xn, t);
                    mpq_set(yn, P[3]); for (int d = 2; d >= 0; d--) { mpq_mul(yn, yn, xn); mpq_add(yn, yn, P[d]); }
                    mpq_t fx; mpq_init(fx); mpq_set_z(fx, f[6]); for (int d = 5; d >= 0; d--) { mpq_mul(fx, fx, xn); mpq_set_z(t, f[d]); mpq_add(fx, fx, t); }
                    mpq_mul(t, yn, yn);
                    if (mpq_equal(t, fx)) {
                        mpz_init_set(out[nnew].a, mpq_numref(xn)); mpz_init_set(out[nnew].c, mpq_denref(xn));
                        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, mpq_denref(xn), 3);
                        mpq_set_z(u, c3); mpq_mul(u, yn, u);
                        mpz_init_set(out[nnew].b, mpq_numref(u));
                        if (mpz_sgn(out[nnew].b) < 0) mpz_neg(out[nnew].b, out[nnew].b);
                        for (int q = 0; q < 4; q++) out[nnew].idx[q] = q < k ? idx[q] : -1;
                        out[nnew].eps = eps | (m[0] << 8) | (k << 12);
                        mpz_clear(c3); nnew++;
                    }
                    mpq_clear(fx); mpq_clear(xn); mpq_clear(yn);
                }
                mpq_clear(sq);
            }
            mpq_clear(disc);
        }
        for (int j = 0; j < 4; j++) mpq_clear(P[j]); for (int d = 0; d < 7; d++) mpq_clear(r[d]);
    }
    for (int i = 0; i < k; i++) { mpq_clear(x[i]); mpq_clear(y[i]); }
    for (int r = 0; r < 4; r++) for (int c = 0; c < 5; c++) mpq_clear(A[r][c]);
    mpq_clear(t); mpq_clear(u);
    return nnew;
}
#define FLUSH_LOC() do { _Pragma("omp critical") { for (int q_ = 0; q_ < nloc; q_++) { if (*nfound == *cfound) { *cfound = *cfound ? 2 * *cfound : 64; *found = realloc(*found, *cfound * sizeof(newpt_t)); } (*found)[(*nfound)++] = loc[q_]; } nloc = 0; } } while (0)



/* ---------------- relations through the points at infinity ----------------
 * sextic with f6 = s^2 (s > 0): contact order k (1..4) at inf_sigma fixes p3 = sigma*s, p2, p1, p0 in turn
 * (the polynomial part of sqrt(f)); 4-k further finite conditions (values / tangents).
 * quintic (f6 = 0): k = 1 means p3 = 0 (a conic); 3 finite conditions.
 * p^2 - f then has degree 6-k; divide by the finite factors -> residual quadratic.
 */

/* returns 1 = candidate, 0 = rejected, -1 = undecided (degenerate mod p) */
static int inf_test_mod(int kp, const int *idx, const int *m, int nfin, int eps, int k, int sigma) {
    uint64_t p = kp < NP ? PR[kp] : BIGP[kp-NP];
    for (int i = 0; i < nfin; i++) if (badr[kp][idx[i]]) return -1;   /* a point without residues mod p: undecided */
    uint64_t A[4][5]; int row = 0;
    uint64_t xs[4];
    for (int i = 0; i < nfin; i++) {
        uint64_t x = xr[kp][idx[i]], y = yr[kp][idx[i]];
        if (eps & (1 << i)) y = (p - y) % p;
        uint64_t d[4] = {y, 0, 0, 0};
        if (m[i] >= 2) { if (y == 0) return -1; y_derivs_mod2(p, fr[kp], x, y, d, (eps & (1 << i)) ? (p - inv2y[kp][idx[i]]) % p : inv2y[kp][idx[i]]); }
        xs[i] = x;
        uint64_t pw[4]; pw[0] = 1; for (int j = 1; j < 4; j++) pw[j] = mulmod(pw[j-1], x, p);
        for (int t = 0; t < m[i]; t++) {
            for (int j = 0; j < 4; j++) {
                if (j < t) { A[row][j] = 0; continue; }
                uint64_t fac = 1; for (int u = 0; u < t; u++) fac *= (j - u);
                A[row][j] = mulmod(fac % p, pw[j-t], p);
            }
            A[row][4] = d[t]; row++;
        }
    }
    /* infinity rows */
    uint64_t pc[4] = {0,0,0,0};   /* fixed coefficients p3, p2, p1, p0 as determined */
    if (f6_is_zero) {
        A[row][0] = A[row][1] = A[row][2] = 0; A[row][3] = 1; A[row][4] = 0; row++;   /* p3 = 0 */
    } else {
        uint64_t s = has_sqrt_f6 ? mpz_fdiv_ui(sqrt_f6z, p) : 0;
        if (sigma < 0) s = (p - s) % p;
        if (s == 0) return -1;
        uint64_t inv2s = invmod(mulmod(2, s, p), p);
        pc[3] = s;
        if (k >= 2) pc[2] = mulmod(fr[kp][5], inv2s, p);
        if (k >= 3) pc[1] = mulmod((fr[kp][4] + p - mulmod(pc[2], pc[2], p)) % p, inv2s, p);
        if (k >= 4) pc[0] = mulmod((fr[kp][3] + p - mulmod(2, mulmod(pc[2], pc[1], p), p)) % p, inv2s, p);
        for (int t = 0; t < k; t++) { for (int j = 0; j < 4; j++) A[row][j] = 0; A[row][3 - t] = 1; A[row][4] = pc[3 - t]; row++; }
    }
    if (row != 4) return -1;
    if (!solve4_mod(p, A)) return -1;
    uint64_t P[4]; for (int j = 0; j < 4; j++) P[j] = A[j][4];
    uint64_t det = A[0][0], det2 = mulmod(det, det, p);
    uint64_t r[7] = {0,0,0,0,0,0,0};
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { r[i+j] += mulmod(P[i], P[j], p); if (r[i+j] >= p) r[i+j] -= p; }
    for (int d = 0; d < 7; d++) r[d] = (r[d] + p - mulmod(det2, fr[kp][d], p)) % p;
    int deg = 6 - k;
    if (f6_is_zero) deg = 5;
    /* the top coefficients must vanish */
    for (int d = deg + 1; d <= 6; d++) if (r[d] != 0) return -1;
    for (int i = 0; i < nfin; i++) for (int t = 0; t < m[i]; t++) {
        for (int d = deg; d >= 1; d--) r[d-1] = (r[d-1] + mulmod(r[d], xs[i], p)) % p;
        for (int d = 0; d < deg; d++) r[d] = r[d+1];
        r[deg] = 0; deg--;
    }
    if (deg != 2) return -1;
    if (r[2] == 0) return -1;
    uint64_t disc = (mulmod(r[1], r[1], p) + p - mulmod(4 % p, mulmod(r[2], r[0], p), p)) % p;
    if (kp < NP) return (QR[kp][disc >> 3] >> (disc & 7)) & 1;
    return legendre_big(disc, p);
}

static int inf_exact(const int *idx, const int *m, int nfin, int eps, int k, int sigma, newpt_t *out) {
    mpq_t A[4][5], t, u, x[4], y[4];
    for (int r = 0; r < 4; r++) for (int c = 0; c < 5; c++) mpq_init(A[r][c]);
    mpq_init(t); mpq_init(u);
    int row = 0;
    for (int i = 0; i < nfin; i++) {
        mpq_init(x[i]); mpq_init(y[i]);
        mpq_set_num(x[i], pts[idx[i]].a); mpq_set_den(x[i], pts[idx[i]].c); mpq_canonicalize(x[i]);
        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, pts[idx[i]].c, 3);
        mpq_set_num(y[i], pts[idx[i]].b); mpq_set_den(y[i], c3); mpq_canonicalize(y[i]); mpz_clear(c3);
        if (eps & (1 << i)) mpq_neg(y[i], y[i]);
        mpq_t d[4], f1, f2, f3, pw[7]; for (int j = 0; j < 4; j++) mpq_init(d[j]); mpq_init(f1); mpq_init(f2); mpq_init(f3);
        for (int j = 0; j < 7; j++) { mpq_init(pw[j]); if (j == 0) mpq_set_ui(pw[j], 1, 1); else mpq_mul(pw[j], pw[j-1], x[i]); }
        mpq_set(d[0], y[i]);
        if (m[i] >= 2) {
            for (int j = 1; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j, 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-1]); mpq_add(f1, f1, t); }
            for (int j = 2; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j*(j-1), 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-2]); mpq_add(f2, f2, t); }
            for (int j = 3; j <= 6; j++) { mpq_set_z(t, f[j]); mpq_set_ui(u, j*(j-1)*(j-2), 1); mpq_mul(t, t, u); mpq_mul(t, t, pw[j-3]); mpq_add(f3, f3, t); }
            mpq_set_ui(u, 2, 1); mpq_mul(u, u, y[i]);
            mpq_div(d[1], f1, u);
            mpq_mul(t, d[1], d[1]); mpq_set_ui(d[2], 2, 1); mpq_mul(t, t, d[2]); mpq_sub(t, f2, t); mpq_div(d[2], t, u);
            mpq_mul(t, d[1], d[2]); mpq_set_ui(d[3], 6, 1); mpq_mul(t, t, d[3]); mpq_sub(t, f3, t); mpq_div(d[3], t, u);
        }
        for (int tt = 0; tt < m[i]; tt++) {
            for (int j = 0; j < 4; j++) {
                if (j < tt) { mpq_set_ui(A[row][j], 0, 1); continue; }
                unsigned long fac = 1; for (int v = 0; v < tt; v++) fac *= (j - v);
                mpq_set_ui(A[row][j], fac, 1); mpq_mul(A[row][j], A[row][j], pw[j-tt]);
            }
            mpq_set(A[row][4], d[tt]); row++;
        }
        for (int j = 0; j < 4; j++) mpq_clear(d[j]); mpq_clear(f1); mpq_clear(f2); mpq_clear(f3); for (int j = 0; j < 7; j++) mpq_clear(pw[j]);
    }
    /* infinity rows */
    mpq_t pc[4]; for (int j = 0; j < 4; j++) { mpq_init(pc[j]); }
    if (f6_is_zero) {
        for (int j = 0; j < 4; j++) mpq_set_ui(A[row][j], 0, 1);
        mpq_set_ui(A[row][3], 1, 1); mpq_set_ui(A[row][4], 0, 1); row++;
    } else {
        mpq_t s, inv2s; mpq_init(s); mpq_init(inv2s);
        mpq_set_z(s, sqrt_f6z); if (sigma < 0) mpq_neg(s, s);
        mpq_set_ui(inv2s, 2, 1); mpq_mul(inv2s, inv2s, s); mpq_inv(inv2s, inv2s);
        mpq_set(pc[3], s);
        if (k >= 2) { mpq_set_z(t, f[5]); mpq_mul(pc[2], t, inv2s); }
        if (k >= 3) { mpq_set_z(t, f[4]); mpq_mul(u, pc[2], pc[2]); mpq_sub(t, t, u); mpq_mul(pc[1], t, inv2s); }
        if (k >= 4) { mpq_set_z(t, f[3]); mpq_mul(u, pc[2], pc[1]); mpq_add(u, u, u); mpq_sub(t, t, u); mpq_mul(pc[0], t, inv2s); }
        for (int tt = 0; tt < k; tt++) { for (int j = 0; j < 4; j++) mpq_set_ui(A[row][j], 0, 1); mpq_set_ui(A[row][3 - tt], 1, 1); mpq_set(A[row][4], pc[3 - tt]); row++; }
        mpq_clear(s); mpq_clear(inv2s);
    }
    int ok = (row == 4);
    for (int col = 0; col < 4 && ok; col++) {
        int piv = -1; for (int r = col; r < 4; r++) if (mpq_sgn(A[r][col])) { piv = r; break; }
        if (piv < 0) { ok = 0; break; }
        if (piv != col) for (int c = 0; c < 5; c++) mpq_swap(A[col][c], A[piv][c]);
        for (int c = 4; c >= col; c--) mpq_div(A[col][c], A[col][c], A[col][col]);
        for (int r = 0; r < 4; r++) if (r != col && mpq_sgn(A[r][col])) { for (int c = 4; c >= col; c--) { mpq_mul(t, A[r][col], A[col][c]); mpq_sub(A[r][c], A[r][c], t); } }
    }
    int nnew = 0;
    if (ok) {
        mpq_t P[4], r[7]; for (int j = 0; j < 4; j++) { mpq_init(P[j]); mpq_set(P[j], A[j][4]); } for (int d = 0; d < 7; d++) mpq_init(r[d]);
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) { mpq_mul(t, P[i], P[j]); mpq_add(r[i+j], r[i+j], t); }
        for (int d = 0; d < 7; d++) { mpq_set_z(t, f[d]); mpq_sub(r[d], r[d], t); }
        int deg = f6_is_zero ? 5 : 6 - k;
        for (int d = deg + 1; d <= 6; d++) if (mpq_sgn(r[d]) != 0) ok = 0;
        for (int i = 0; i < nfin && ok; i++) for (int tt = 0; tt < m[i]; tt++) {
            for (int d = deg; d >= 1; d--) { mpq_mul(t, r[d], x[i]); mpq_add(r[d-1], r[d-1], t); }
            if (mpq_sgn(r[0]) != 0) ok = 0;
            for (int d = 0; d < deg; d++) mpq_set(r[d], r[d+1]);
            mpq_set_ui(r[deg], 0, 1); deg--;
        }
        if (ok && deg == 2 && mpq_sgn(r[2]) != 0) {
            mpq_t disc; mpq_init(disc);
            mpq_mul(disc, r[1], r[1]); mpq_mul(t, r[2], r[0]); mpq_set_ui(u, 4, 1); mpq_mul(t, t, u); mpq_sub(disc, disc, t);
            if (mpq_sgn(disc) >= 0 && mpz_perfect_square_p(mpq_numref(disc)) && mpz_perfect_square_p(mpq_denref(disc))) {
                mpq_t sq; mpq_init(sq); mpz_sqrt(mpq_numref(sq), mpq_numref(disc)); mpz_sqrt(mpq_denref(sq), mpq_denref(disc)); mpq_canonicalize(sq);
                for (int s2 = 0; s2 < 2; s2++) {
                    mpq_t xn, yn; mpq_init(xn); mpq_init(yn);
                    if (s2) mpq_sub(xn, r[1], sq); else mpq_add(xn, r[1], sq);
                    mpq_neg(xn, xn); mpq_set_ui(u, 2, 1); mpq_mul(t, r[2], u); mpq_div(xn, xn, t);
                    mpq_set(yn, P[3]); for (int d = 2; d >= 0; d--) { mpq_mul(yn, yn, xn); mpq_add(yn, yn, P[d]); }
                    mpq_t fx; mpq_init(fx); mpq_set_z(fx, f[6]); for (int d = 5; d >= 0; d--) { mpq_mul(fx, fx, xn); mpq_set_z(t, f[d]); mpq_add(fx, fx, t); }
                    mpq_mul(t, yn, yn);
                    if (mpq_equal(t, fx)) {
                        mpz_init_set(out[nnew].a, mpq_numref(xn)); mpz_init_set(out[nnew].c, mpq_denref(xn));
                        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, mpq_denref(xn), 3);
                        mpq_set_z(u, c3); mpq_mul(u, yn, u);
                        mpz_init_set(out[nnew].b, mpq_numref(u));
                        if (mpz_sgn(out[nnew].b) < 0) mpz_neg(out[nnew].b, out[nnew].b);
                        out[nnew].ys = mpq_sgn(yn);
                        for (int q = 0; q < 4; q++) out[nnew].idx[q] = q < nfin ? idx[q] : -1;
                        out[nnew].eps = eps | (k << 8) | ((sigma > 0 ? 1 : 2) << 12);
                        mpz_clear(c3); nnew++;
                    }
                    mpq_clear(fx); mpq_clear(xn); mpq_clear(yn);
                }
                mpq_clear(sq);
            }
            mpq_clear(disc);
        }
        for (int j = 0; j < 4; j++) mpq_clear(P[j]); for (int d = 0; d < 7; d++) mpq_clear(r[d]);
    }
    for (int j = 0; j < 4; j++) mpq_clear(pc[j]);
    for (int i = 0; i < nfin; i++) { mpq_clear(x[i]); mpq_clear(y[i]); }
    for (int r = 0; r < 4; r++) for (int c = 0; c < 5; c++) mpq_clear(A[r][c]);
    mpq_clear(t); mpq_clear(u);
    return nnew;
}

static int is_weier(int i) { return mpz_sgn(pts[i].b) == 0; }

/* try one infinity relation: modular filters then exact */
static void inf_try(const int *idx, const int *m, int nfin, int eps, int k, int sigma, newpt_t *loc, int *nloc, uint64_t *cnt) {
    for (int i = 0; i < nfin; i++) if (m[i] >= 2 && is_weier(idx[i])) return;
    int ok = 1;
    for (int kp = 0; kp < NP + NBIG && ok; kp++) { int r = inf_test_mod(kp, idx, m, nfin, eps, k, sigma); if (r == 0) ok = 0; }
    (*cnt)++;
    if (!ok) return;
    newpt_t out[2]; int nn = inf_exact(idx, m, nfin, eps, k, sigma, out);
    for (int q = 0; q < nn; q++) if (*nloc < 64) loc[(*nloc)++] = out[q];
}



/* ---------------- 6-point relations: phi = A(x) + B(x) y, deg A <= 4, deg B <= 1 (sections of 4D) ----------------
 * six known points (x_i, eps_i y_i) determine phi up to scalar; A^2 - B^2 f (degree 8) has the six x_i as roots;
 * the residual quadratic gives up to two new points with y = -A(x)/B(x).
 */
/* nullspace vector of a 6x7 matrix mod p (rows: conditions). returns 0 if rank < 6 (undecided) */
static int nullvec7_mod(uint64_t p, uint64_t M[6][7], uint64_t v[7]) {
    int piv[6]; int r = 0; int used[7] = {0,0,0,0,0,0,0};
    for (int col = 0; col < 7 && r < 6; col++) {
        int pr = -1; for (int i = r; i < 6; i++) if (M[i][col]) { pr = i; break; }
        if (pr < 0) continue;
        if (pr != r) for (int k = 0; k < 7; k++) { uint64_t t = M[r][k]; M[r][k] = M[pr][k]; M[pr][k] = t; }
        uint64_t inv = invmod(M[r][col], p);
        for (int k = 0; k < 7; k++) M[r][k] = mulmod(M[r][k], inv, p);
        for (int i = 0; i < 6; i++) if (i != r && M[i][col]) { uint64_t fct = M[i][col]; for (int k = 0; k < 7; k++) M[i][k] = (M[i][k] + p - mulmod(fct, M[r][k], p)) % p; }
        piv[r] = col; used[col] = 1; r++;
    }
    if (r < 6) return 0;
    int freec = -1; for (int c = 0; c < 7; c++) if (!used[c]) { freec = c; break; }
    for (int k = 0; k < 7; k++) v[k] = 0;
    v[freec] = 1;
    for (int i = 0; i < 6; i++) v[piv[i]] = (p - M[i][freec]) % p;
    return 1;
}
/* returns 1 = candidate, 0 = rejected, -1 = undecided */
static int six_test_mod(int kp, const int *idx, int eps) {
    uint64_t p = kp < NP ? PR[kp] : BIGP[kp-NP];
    uint64_t M[6][7]; uint64_t xs[6];
    for (int i = 0; i < 6; i++) {
        uint64_t x = xr[kp][idx[i]], y = yr[kp][idx[i]]; if (eps & (1 << i)) y = (p - y) % p;
        xs[i] = x;
        uint64_t pw = 1;
        for (int j = 0; j < 5; j++) { M[i][j] = pw; pw = mulmod(pw, x, p); }
        M[i][5] = y; M[i][6] = mulmod(x, y, p);
    }
    for (int i = 0; i < 6; i++) for (int j = i+1; j < 6; j++) if (xs[i] == xs[j]) return -1;
    uint64_t v[7];
    if (!nullvec7_mod(p, M, v)) return -1;
    /* A = v0..v4, B = v5 + v6 x ; r = A^2 - B^2 f, degree <= 8 */
    uint64_t r[9] = {0,0,0,0,0,0,0,0,0};
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) r[i+j] = (r[i+j] + mulmod(v[i], v[j], p)) % p;
    uint64_t B2[3]; B2[0] = mulmod(v[5], v[5], p); B2[1] = mulmod(2, mulmod(v[5], v[6], p), p); B2[2] = mulmod(v[6], v[6], p);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 7; j++) r[i+j] = (r[i+j] + p - mulmod(B2[i], fr[kp][j], p)) % p;
    int deg = 8;
    if (r[8] == 0) return -1;     /* zeros at infinity: skip */
    for (int i = 0; i < 6; i++) {
        for (int d = deg; d >= 1; d--) r[d-1] = (r[d-1] + mulmod(r[d], xs[i], p)) % p;
        if (r[0] != 0) return -1;  /* not a root mod p (should not happen) */
        for (int d = 0; d < deg; d++) r[d] = r[d+1];
        r[deg] = 0; deg--;
    }
    if (r[2] == 0) return -1;
    uint64_t disc = (mulmod(r[1], r[1], p) + p - mulmod(4 % p, mulmod(r[2], r[0], p), p)) % p;
    if (kp < NP) return (QR[kp][disc >> 3] >> (disc & 7)) & 1;
    return legendre_big(disc, p);
}
static int six_exact(const int *idx, int eps, newpt_t *out) {
    mpq_t M[6][8], t, u, x[6], y[6];
    for (int i = 0; i < 6; i++) for (int j = 0; j < 8; j++) mpq_init(M[i][j]);
    mpq_init(t); mpq_init(u);
    for (int i = 0; i < 6; i++) {
        mpq_init(x[i]); mpq_init(y[i]);
        mpq_set_num(x[i], pts[idx[i]].a); mpq_set_den(x[i], pts[idx[i]].c); mpq_canonicalize(x[i]);
        mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, pts[idx[i]].c, 3);
        mpq_set_num(y[i], pts[idx[i]].b); mpq_set_den(y[i], c3); mpq_canonicalize(y[i]); mpz_clear(c3);
        if (eps & (1 << i)) mpq_neg(y[i], y[i]);
        mpq_set_ui(M[i][0], 1, 1);
        for (int j = 1; j < 5; j++) mpq_mul(M[i][j], M[i][j-1], x[i]);
        mpq_set(M[i][5], y[i]); mpq_mul(M[i][6], x[i], y[i]);
    }
    /* row reduce 6x7, find nullspace */
    int piv[6]; int r = 0; int used[7] = {0,0,0,0,0,0,0};
    for (int col = 0; col < 7 && r < 6; col++) {
        int pr = -1; for (int i = r; i < 6; i++) if (mpq_sgn(M[i][col])) { pr = i; break; }
        if (pr < 0) continue;
        if (pr != r) for (int k = 0; k < 7; k++) mpq_swap(M[r][k], M[pr][k]);
        for (int k = 6; k >= col; k--) mpq_div(M[r][k], M[r][k], M[r][col]);
        for (int i = 0; i < 6; i++) if (i != r && mpq_sgn(M[i][col])) { for (int k = 6; k >= col; k--) { mpq_mul(t, M[i][col], M[r][k]); mpq_sub(M[i][k], M[i][k], t); } }
        piv[r] = col; used[col] = 1; r++;
    }
    int nnew = 0;
    if (r == 6) {
        int freec = -1; for (int c = 0; c < 7; c++) if (!used[c]) { freec = c; break; }
        mpq_t v[7]; for (int k = 0; k < 7; k++) { mpq_init(v[k]); mpq_set_ui(v[k], 0, 1); }
        mpq_set_ui(v[freec], 1, 1);
        for (int i = 0; i < 6; i++) { mpq_neg(v[piv[i]], M[i][freec]); }
        /* r = A^2 - B^2 f */
        mpq_t rr[9]; for (int d = 0; d < 9; d++) mpq_init(rr[d]);
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) { mpq_mul(t, v[i], v[j]); mpq_add(rr[i+j], rr[i+j], t); }
        mpq_t B2[3]; for (int i = 0; i < 3; i++) mpq_init(B2[i]);
        mpq_mul(B2[0], v[5], v[5]); mpq_mul(B2[1], v[5], v[6]); mpq_add(B2[1], B2[1], B2[1]); mpq_mul(B2[2], v[6], v[6]);
        for (int i = 0; i < 3; i++) for (int j = 0; j < 7; j++) { mpq_set_z(t, f[j]); mpq_mul(t, t, B2[i]); mpq_sub(rr[i+j], rr[i+j], t); }
        int deg = 8; int ok = mpq_sgn(rr[8]) != 0;
        for (int i = 0; i < 6 && ok; i++) {
            for (int d = deg; d >= 1; d--) { mpq_mul(t, rr[d], x[i]); mpq_add(rr[d-1], rr[d-1], t); }
            if (mpq_sgn(rr[0]) != 0) ok = 0;
            for (int d = 0; d < deg; d++) mpq_set(rr[d], rr[d+1]);
            mpq_set_ui(rr[deg], 0, 1); deg--;
        }
        if (ok && deg == 2 && mpq_sgn(rr[2]) != 0) {
            mpq_t disc; mpq_init(disc);
            mpq_mul(disc, rr[1], rr[1]); mpq_mul(t, rr[2], rr[0]); mpq_set_ui(u, 4, 1); mpq_mul(t, t, u); mpq_sub(disc, disc, t);
            if (mpq_sgn(disc) >= 0 && mpz_perfect_square_p(mpq_numref(disc)) && mpz_perfect_square_p(mpq_denref(disc))) {
                mpq_t sq; mpq_init(sq); mpz_sqrt(mpq_numref(sq), mpq_numref(disc)); mpz_sqrt(mpq_denref(sq), mpq_denref(disc)); mpq_canonicalize(sq);
                for (int s2 = 0; s2 < 2; s2++) {
                    mpq_t xn, yn, Av, Bv; mpq_init(xn); mpq_init(yn); mpq_init(Av); mpq_init(Bv);
                    if (s2) mpq_sub(xn, rr[1], sq); else mpq_add(xn, rr[1], sq);
                    mpq_neg(xn, xn); mpq_set_ui(u, 2, 1); mpq_mul(t, rr[2], u); mpq_div(xn, xn, t);
                    /* y = -A(x)/B(x) */
                    mpq_set(Av, v[4]); for (int d = 3; d >= 0; d--) { mpq_mul(Av, Av, xn); mpq_add(Av, Av, v[d]); }
                    mpq_mul(Bv, v[6], xn); mpq_add(Bv, Bv, v[5]);
                    if (mpq_sgn(Bv) != 0) {
                        mpq_div(yn, Av, Bv); mpq_neg(yn, yn);
                        mpq_t fx; mpq_init(fx); mpq_set_z(fx, f[6]); for (int d = 5; d >= 0; d--) { mpq_mul(fx, fx, xn); mpq_set_z(t, f[d]); mpq_add(fx, fx, t); }
                        mpq_mul(t, yn, yn);
                        if (mpq_equal(t, fx)) {
                            mpz_init_set(out[nnew].a, mpq_numref(xn)); mpz_init_set(out[nnew].c, mpq_denref(xn));
                            mpz_t c3; mpz_init(c3); mpz_pow_ui(c3, mpq_denref(xn), 3);
                            mpq_set_z(u, c3); mpq_mul(u, yn, u);
                            mpz_init_set(out[nnew].b, mpq_numref(u));
                            if (mpz_sgn(out[nnew].b) < 0) mpz_neg(out[nnew].b, out[nnew].b);
                            for (int q = 0; q < 4; q++) out[nnew].idx[q] = idx[q]; out[nnew].eps = eps | (6 << 12);
                            mpz_clear(c3); nnew++;
                        }
                        mpq_clear(fx);
                    }
                    mpq_clear(xn); mpq_clear(yn); mpq_clear(Av); mpq_clear(Bv);
                }
                mpq_clear(sq);
            }
            mpq_clear(disc);
        }
        for (int d = 0; d < 9; d++) mpq_clear(rr[d]); for (int i = 0; i < 3; i++) mpq_clear(B2[i]); for (int k = 0; k < 7; k++) mpq_clear(v[k]);
    }
    for (int i = 0; i < 6; i++) { mpq_clear(x[i]); mpq_clear(y[i]); }
    for (int i = 0; i < 6; i++) for (int j = 0; j < 8; j++) mpq_clear(M[i][j]);
    mpq_clear(t); mpq_clear(u);
    return nnew;
}
/* random sample of nsamp 6-subsets (distinct x, no Weierstrass points) with random signs */
static void six_pass(int n, uint64_t nsamp, newpt_t **found, int *nfound, int *cfound, uint64_t *nexact) {
    uint64_t cnt = 0, nex = 0;
    #pragma omp parallel reduction(+:cnt, nex)
    {
        newpt_t loc[64]; int nloc = 0;
        uint64_t rng = 0x9E3779B97F4A7C15ULL * (omp_get_thread_num() + 1) ^ 0xD1B54A32D192ED03ULL;
        #pragma omp for schedule(dynamic, 4096)
        for (uint64_t s = 0; s < nsamp; s++) {
            int idx[6]; int okset = 1;
            for (int i = 0; i < 6; i++) {
                rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                idx[i] = (int)(rng % (uint64_t)n);
                for (int j = 0; j < i; j++) if (idx[j] == idx[i]) okset = 0;
                if (is_weier(idx[i])) okset = 0;
            }
            if (!okset) continue;
            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
            int eps = (int)(rng & 31) << 1;     /* sign of point 0 fixed */
            int ok = 1, dec = 0;
            for (int kp = 0; kp < NP + NBIG && ok; kp++) { int r = six_test_mod(kp, idx, eps); if (r == 0) ok = 0; else if (r == 1) dec++; }
            cnt++;
            if (!ok || dec < 4) continue;
            nex++;
            newpt_t out[2]; int nn = six_exact(idx, eps, out);
            for (int q = 0; q < nn; q++) if (nloc < 64) loc[nloc++] = out[q];
            if (nloc >= 60) {
                #pragma omp critical
                { for (int q = 0; q < nloc; q++) { if (*nfound == *cfound) { *cfound = *cfound ? 2 * *cfound : 64; *found = realloc(*found, *cfound * sizeof(newpt_t)); } (*found)[(*nfound)++] = loc[q]; } nloc = 0; }
            }
        }
        if (nloc) {
            #pragma omp critical
            { for (int q = 0; q < nloc; q++) { if (*nfound == *cfound) { *cfound = *cfound ? 2 * *cfound : 64; *found = realloc(*found, *cfound * sizeof(newpt_t)); } (*found)[(*nfound)++] = loc[q]; } nloc = 0; }
        }
    }
    *nexact += nex;
    fprintf(stderr, "six-point sample: %lu subsets tested, %lu exact checks\n", cnt, nex);
}


/* ================= AVX-512 IFMA vectorised 4-point relation test =================
 * lanes = 8 subsets {i0,i1,i2,i3+l}, l = 0..7 (i3 consecutive); 16-bit primes.  Lazy reduction: residues are kept
 * "half-reduced" in [0, 2p) (or larger where annotated), products are accumulated with vpmadd52lo (exact below 2^52)
 * and brought back to [0, 2p) by a single Barrett step with M = floor(2^52/p) (valid for any x < 2^52); a full
 * reduction to [0, p) is only needed for the discriminant table lookup.  The quotient of p^2 - f by the monic
 * quartic Q = prod (x - x_i) only needs the coefficients of x^4..x^6 of p^2 - f (the remainder is zero), i.e. only
 * p_1, p_2, p_3 and Q_3, Q_2, which halves the interpolation work.
 */
#define NVP 8                      /* primes handled in the vector stage */
typedef struct { uint64_t p, M; __m512i vp, v2p, v8p, vM; } vprime_t;
static vprime_t VP[NVP];

static inline __m512i vred(__m512i x, const vprime_t *P) {   /* x < 2^52  ->  x mod p, in [0, 2p) */
    __m512i q = _mm512_madd52hi_epu64(_mm512_setzero_si512(), x, P->vM);
    return _mm512_sub_epi64(x, _mm512_madd52lo_epu64(_mm512_setzero_si512(), q, P->vp));
}
static inline __m512i vfull(__m512i x, const vprime_t *P) {  /* [0, 2p)  ->  [0, p) */
    __mmask8 ge = _mm512_cmpge_epu64_mask(x, P->vp);
    return _mm512_mask_sub_epi64(x, ge, x, P->vp);
}
static inline __m512i vmul(__m512i a, __m512i b, const vprime_t *P) {   /* a * b < 2^52  ->  [0, 2p) */
    return vred(_mm512_madd52lo_epu64(_mm512_setzero_si512(), a, b), P);
}
static inline __m512i vmla(__m512i c, __m512i a, __m512i b) { return _mm512_madd52lo_epu64(c, a, b); }   /* c + a*b, no reduction */

typedef struct { int32_t i0, i1, i2, i3; int32_t eps; } surv_t;

/* scalar tail for a survivor of the vector stage: primes NVP.. then exact */
static int tail_check(int i0, int i1, int i2, int i3, int eps, newpt_t *out) {
    int idx[4] = {i0, i1, i2, i3};
    for (int k = NVP; k < NP + NBIG; k++) {
        sub_t S; if (!prep_subset(k, idx, &S)) continue;
        uint64_t disc; if (!test_pattern(k, &S, eps, &disc)) return -1;
    }
    return exact_check(idx, eps, out);
}

/* per-lane data of one prime for a block of 8 subsets: half-reduced y_i L_i (degrees 1..3) and their negatives,
   -Q_3, -Q_2 (half-reduced), and the collision mask */
typedef struct { __m512i yL[4][4], nyL[4][4], NQ3, NQ2; __mmask8 ok; } vblock_t;

static inline void vblock_setup(vblock_t *B, int k, int i0, int i1, int i2, int i3, __mmask8 lanes,
                                const uint64_t c3[4], const uint64_t cq[3][3], const uint64_t dinv3[3], const uint64_t y012[3]) {
    const vprime_t *P = &VP[k];
    __m512i X3 = _mm512_maskz_loadu_epi64(lanes, xr[k] + i3), Y3 = _mm512_maskz_loadu_epi64(lanes, yr[k] + i3);
    __m512i inv0 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i0 * ncap + i3);
    __m512i inv1 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i1 * ncap + i3);
    __m512i inv2 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i2 * ncap + i3);
    B->ok = _mm512_test_epi64_mask(inv0, inv0) & _mm512_test_epi64_mask(inv1, inv1) & _mm512_test_epi64_mask(inv2, inv2);
    __m512i NX3 = _mm512_sub_epi64(P->vp, X3);                       /* -x3, in (0, p] */
    __m512i c0 = _mm512_set1_epi64(c3[0]), c1 = _mm512_set1_epi64(c3[1]), c2 = _mm512_set1_epi64(c3[2]);
    /* Q = C (x - x3): Q3 = c2 - x3, Q2 = c1 - c2 x3 */
    B->NQ3 = _mm512_sub_epi64(P->v2p, _mm512_add_epi64(c2, NX3));      /* -(c2 + NX3), c2 + NX3 in (0, 2p) */
    B->NQ2 = _mm512_sub_epi64(P->v2p, vred(vmla(c1, c2, NX3), P));
    __m512i invv[3] = {inv0, inv1, inv2};
    for (int i = 0; i < 3; i++) {
        __m512i q0 = _mm512_set1_epi64(cq[i][0]), q1 = _mm512_set1_epi64(cq[i][1]), q2 = _mm512_set1_epi64(cq[i][2]);
        /* L_i = (C/(x - x_i)) (x - x3) / (prod_{j<3, j != i} (x_i - x_j) (x_i - x3)), scaled by y_i; degrees 1..3 */
        __m512i L3 = q2, L2 = vred(vmla(q1, q2, NX3), P), L1 = vred(vmla(q0, q1, NX3), P);
        __m512i sc = vmul(vmul(_mm512_set1_epi64(dinv3[i]), invv[i], P), _mm512_set1_epi64(y012[i]), P);
        B->yL[i][3] = vmul(L3, sc, P); B->yL[i][2] = vmul(L2, sc, P); B->yL[i][1] = vmul(L1, sc, P);
    }
    {   /* L_3 = C / prod_j (x3 - x_j) = -inv0 inv1 inv2 C, scaled by y3 */
        __m512i d3 = _mm512_sub_epi64(P->v2p, vmul(vmul(inv0, inv1, P), inv2, P));
        __m512i sc = vmul(d3, Y3, P);
        B->yL[3][3] = sc; B->yL[3][2] = vmul(c2, sc, P); B->yL[3][1] = vmul(c1, sc, P);
    }
    for (int i = 0; i < 4; i++) for (int d = 1; d < 4; d++) B->nyL[i][d] = _mm512_sub_epi64(P->v2p, B->yL[i][d]);
    (void)c0;
}

/* core of the modular tests: p_1, p_2, p_3 (< 2^24 each) are the top coefficients of the interpolating cubic,
   NQ3, NQ2 = -Q_3, -Q_2 (half-reduced) of the monic quartic Q = prod (x - x_i) (with multiplicity); computes the
   quotient q = (p^2 - f) / Q and returns the mask of lanes whose disc(q) is a square mod p or whose q_2 is 0 mod p */
static inline __mmask8 vdisc_test(int k, __m512i P1, __m512i P2, __m512i P3, __m512i NQ3, __m512i NQ2, __mmask8 lanes) {
    const vprime_t *P = &VP[k];
    __m512i P3x2 = _mm512_add_epi64(P3, P3);
    __m512i NF6 = _mm512_set1_epi64(P->p - fr[k][6]), NF5 = _mm512_set1_epi64(P->p - fr[k][5]), NF4 = _mm512_set1_epi64(P->p - fr[k][4]);
    __m512i r6 = vmla(NF6, P3, P3);                                   /* coefficients of x^6, x^5, x^4 of p^2 - f, < 2^51 */
    __m512i r5 = vmla(NF5, P3x2, P2);
    __m512i r4 = vmla(vmla(NF4, P3x2, P1), P2, P2);
    __m512i qq2 = vred(r6, P);                                        /* quotient by the monic Q: q2 x^2 + q1 x + q0 */
    r5 = vmla(r5, qq2, NQ3); r4 = vmla(r4, qq2, NQ2);
    __m512i qq1 = vred(r5, P);
    r4 = vmla(r4, qq1, NQ3);
    __m512i qq0 = vred(r4, P);
    __m512i u = vmul(qq2, qq0, P);                                    /* disc = q1^2 - 4 q2 q0 */
    __m512i w = _mm512_sub_epi64(P->v8p, _mm512_slli_epi64(u, 2));
    __m512i disc = vfull(vred(vmla(w, qq1, qq1), P), P);
    __mmask8 q2z = _mm512_cmpeq_epi64_mask(qq2, _mm512_setzero_si512()) | _mm512_cmpeq_epi64_mask(qq2, P->vp);
    __m256i g = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lanes, disc, (const int*)QRB[k], 1);
    __mmask8 pass = _mm256_test_epi32_mask(g, _mm256_set1_epi32(0xff));
    return pass | q2z;
}

/* discriminant test of sign pattern eps (bit i = sign of point i, i = 0..2; point 3 positive) for all lanes:
   returns the mask of lanes that pass (QR), are undecided (leading coefficient 0 mod p) or collide */
static inline __mmask8 vblock_test(const vblock_t *B, int k, int eps, __mmask8 lanes) {
    __m512i Pc[4];
    for (int d = 1; d < 4; d++) {
        __m512i acc = B->yL[3][d];
        acc = _mm512_add_epi64(acc, (eps & 1) ? B->nyL[0][d] : B->yL[0][d]);
        acc = _mm512_add_epi64(acc, (eps & 2) ? B->nyL[1][d] : B->yL[1][d]);
        acc = _mm512_add_epi64(acc, (eps & 4) ? B->nyL[2][d] : B->yL[2][d]);
        Pc[d] = acc;                                                 /* < 8p */
    }
    return vdisc_test(k, Pc[1], Pc[2], Pc[3], B->NQ3, B->NQ2, lanes) | (__mmask8)~B->ok;
}

static int KRES = 0;   /* -K, four-point pass only: i1 < KRES (i0 < i1 are then two of the KRES lowest-height seed points; KRES <= seed size, see main); 0 = unrestricted */

/* vectorised pass over all 4-subsets with i3 >= max(nstart, i2+1) (subsets with at least one new point) */
static void four_pass_vec(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *nsub_out, uint64_t *nsurv_out, uint64_t *nexact_out) {
    uint64_t nsub = 0, nsurv = 0, nexact = 0;
    #pragma omp parallel reduction(+:nsub, nsurv, nexact)
    {
        newpt_t loc[64]; int nloc = 0;
        surv_t *sv = malloc(sizeof(surv_t) * 8192); int nsv = 0;
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 < n; i0++) {
            for (int i1 = i0 + 1; i1 < n && (KRES <= 0 || i1 < KRES); i1++) for (int i2 = i1 + 1; i2 < n; i2++) {
                int i3start = i2 + 1; if (i3start < nstart) i3start = nstart;   /* all four old -> skip */
                if (i3start >= n) continue;
                /* per-prime scalar data for the triple */
                uint64_t c3[NVP][4], cq[NVP][3][3], dinv3[NVP][3], y012[NVP][3];
                int tripok[NVP];   /* per prime: 0 when two of the triple's x-values coincide mod p or a point has no residues mod p -- that prime is undecided for the triple's subsets, the others still filter (a single flag used to drop the whole triple: Codex round 3, PR #15) */
                for (int k = 0; k < NVP; k++) {
                    uint64_t p = VP[k].p; double pinv = PINV[k];
                    uint64_t x[3] = {xr[k][i0], xr[k][i1], xr[k][i2]};
                    y012[k][0] = yr[k][i0]; y012[k][1] = yr[k][i1]; y012[k][2] = yr[k][i2];
                    /* cubic C = (x-x0)(x-x1)(x-x2): coefficients c3[k][0..3] */
                    uint64_t C[4] = {1, 0, 0, 0}; int dq = 0;
                    for (int i = 0; i < 3; i++) { uint64_t nx = x[i] ? p - x[i] : 0; for (int d = dq + 1; d >= 1; d--) C[d] = redp(C[d-1] + C[d] * nx, p, pinv); C[0] = redp(C[0] * nx, p, pinv); dq++; }
                    memcpy(c3[k], C, sizeof C);
                    /* quadratics C/(x - x_i) and partial inverse denominators over the triple */
                    int idx3[3] = {i0, i1, i2};
                    tripok[k] = 1;
                    for (int i = 0; i < 3; i++) {
                        uint64_t carry = C[3]; cq[k][i][2] = carry;
                        for (int d = 2; d >= 1; d--) { carry = redp(C[d] + carry * x[i], p, pinv); cq[k][i][d-1] = carry; }
                        uint64_t di = 1;
                        for (int j = 0; j < 3; j++) if (j != i) { uint64_t iv = invd[k][(size_t)idx3[i] * ncap + idx3[j]]; if (!iv) tripok[k] = 0; di = redp(di * iv, p, pinv); }
                        dinv3[k][i] = di;
                    }
                }
                for (int i3 = i3start; i3 < n; i3 += 8) {
                    int nl = n - i3 < 8 ? n - i3 : 8;
                    __mmask8 lanes = (__mmask8)((1u << nl) - 1);
                    __mmask8 alive[8]; for (int e = 0; e < 8; e++) alive[e] = lanes;
                    nsub += nl;
                    for (int k = 0; k < NVP && (alive[0]|alive[1]|alive[2]|alive[3]|alive[4]|alive[5]|alive[6]|alive[7]); k++) {
                        if (!tripok[k]) continue;   /* undecided for this prime */
                        vblock_t B;
                        vblock_setup(&B, k, i0, i1, i2, i3, lanes, c3[k], cq[k], dinv3[k], y012[k]);
                        if (!B.ok) continue;   /* all lanes collide mod this prime: undecided for this prime */
                        #pragma GCC unroll 8
                        for (int eps = 0; eps < 8; eps++) {
                            if (!alive[eps]) continue;
                            alive[eps] &= vblock_test(&B, k, eps, lanes);
                        }
                    }
                    /* survivors -> scalar tail */
                    for (int eps = 0; eps < 8; eps++) {
                        __mmask8 m = alive[eps];
                        while (m) { int l = __builtin_ctz(m); m &= m - 1; if (nsv < 8192) { sv[nsv].i0 = i0; sv[nsv].i1 = i1; sv[nsv].i2 = i2; sv[nsv].i3 = i3 + l; sv[nsv].eps = eps; nsv++; } }
                    }
                    if (nsv >= 8000) {
                        for (int q = 0; q < nsv; q++) { newpt_t out[2]; int nn = tail_check(sv[q].i0, sv[q].i1, sv[q].i2, sv[q].i3, sv[q].eps, out); nsurv++; nexact += nn >= 0; for (int u = 0; u < nn; u++) if (nloc < 64) loc[nloc++] = out[u]; if (nloc >= 62) { FLUSH_LOC(); } }
                        nsv = 0;
                    }
                }
            }
            for (int q = 0; q < nsv; q++) { newpt_t out[2]; int nn = tail_check(sv[q].i0, sv[q].i1, sv[q].i2, sv[q].i3, sv[q].eps, out); nsurv++; nexact += nn >= 0; for (int u = 0; u < nn; u++) if (nloc < 64) loc[nloc++] = out[u]; if (nloc >= 62) { FLUSH_LOC(); } }
            nsv = 0;
            if (nloc) { FLUSH_LOC(); }
        }
        free(sv);
    }
    *nsub_out += nsub; *nsurv_out += nsurv; *nexact_out += nexact;
}

/* ---------------- index sieve for the four-point pass ----------------
 * Let p be an odd prime of good reduction at which f has a root r.  Over F_p the curve is s^2 = g(t) with g of degree 5
 * (t = 1/(x - r), s = y t^3; the Weierstrass point above r goes to the point inf' at infinity), and P -> [P - inf'] maps
 * its points to the Jacobian J(F_p).  The six points of the section by a cubic y = p(x) satisfy sum [P_i - inf'] = 0, since
 * their sum is linearly equivalent to 3K ~ 6 inf'.  Let H be the subgroup generated by the classes [P - inf'], P in
 * C(F_p), and W = { [A - inf'] + [B - inf'] : A, B in C(F_p) }, a subset of H stable under negation.  A 4-subset with
 * signs can have two rational residual points only if the sum of its four classes lies in W.  When H is generated by
 * two classes G1, G2 (G1 of order M1, m2 G2 = t G1 with m2 > 0 minimal), every class is x1 G1 + x2 G2 with 0 <= x1 < M1,
 * 0 <= x2 < m2 in one way, and the classes are numbered x1 + M1 x2.  For a fixed triple with signs the condition is, for
 * all fourth points and both of their signs at once, one row of bits that depends only on the number of the class of the
 * triple; the rows of the primes are ANDed.  Nothing is undecided modulo p: the statement holds for every 4-subset,
 * whatever coincides modulo p.  The survivors go to the exact check.  A relation among six known points is verified
 * exactly once: the other 14 of its 4-subsets have the same cubic, so their residual points are known.  The tangent
 * relations and the relations through the points at infinity are sieved in the same way (a point of contact of order m
 * counts m times). */
typedef struct { int deg; uint32_t c[14]; } fpoly;   /* deg = -1: the zero polynomial; coefficients in [0, p) */
static inline void fp_trim(fpoly *a) { while (a->deg >= 0 && a->c[a->deg] == 0) a->deg--; }
static void fp_mul(fpoly *r, const fpoly *a, const fpoly *b, uint32_t p) {
    fpoly t; if (a->deg < 0 || b->deg < 0) { r->deg = -1; return; }
    t.deg = a->deg + b->deg; for (int i = 0; i <= t.deg; i++) t.c[i] = 0;
    for (int i = 0; i <= a->deg; i++) for (int j = 0; j <= b->deg; j++) t.c[i+j] = (t.c[i+j] + a->c[i] * b->c[j]) % p;
    *r = t;
}
static void fp_add(fpoly *r, const fpoly *a, const fpoly *b, uint32_t p) {
    fpoly t; t.deg = a->deg > b->deg ? a->deg : b->deg;
    for (int i = 0; i <= t.deg; i++) t.c[i] = ((i <= a->deg ? a->c[i] : 0) + (i <= b->deg ? b->c[i] : 0)) % p;
    fp_trim(&t); *r = t;
}
static void fp_sub(fpoly *r, const fpoly *a, const fpoly *b, uint32_t p) {
    fpoly t; t.deg = a->deg > b->deg ? a->deg : b->deg;
    for (int i = 0; i <= t.deg; i++) t.c[i] = ((i <= a->deg ? a->c[i] : 0) + p - (i <= b->deg ? b->c[i] : 0)) % p;
    fp_trim(&t); *r = t;
}
static void fp_scale(fpoly *a, uint32_t s, uint32_t p) { for (int i = 0; i <= a->deg; i++) a->c[i] = a->c[i] * s % p; }
/* a = q b + r (b != 0); q or r may be NULL, and may be a or b */
static void fp_divrem(fpoly *q, fpoly *r, const fpoly *a, const fpoly *b, uint32_t p, const uint32_t *inv) {
    fpoly rr = *a, qq; int db = b->deg; qq.deg = a->deg - db;
    if (qq.deg < 0) qq.deg = -1;
    uint32_t li = inv[b->c[db]];
    for (int d = rr.deg; d >= db; d--) {
        uint32_t t = rr.c[d] * li % p; qq.c[d - db] = t;
        if (t) for (int j = 0; j <= db; j++) rr.c[d - db + j] = (rr.c[d - db + j] + (p - t) * b->c[j]) % p;
    }
    if (rr.deg >= db) rr.deg = db - 1;
    fp_trim(&rr);
    if (q) *q = qq;
    if (r) *r = rr;
}
/* g = s a + t b, g monic (g = 0 only if a = b = 0) */
static void fp_xgcd(fpoly *g, fpoly *s, fpoly *t, const fpoly *a, const fpoly *b, uint32_t p, const uint32_t *inv) {
    fpoly r0 = *a, r1 = *b, s0, s1, t0, t1; s0.deg = 0; s0.c[0] = 1; s1.deg = -1; t0.deg = -1; t1.deg = 0; t1.c[0] = 1;
    while (r1.deg >= 0) {
        fpoly q, r, x; fp_divrem(&q, &r, &r0, &r1, p, inv);
        r0 = r1; r1 = r;
        fp_mul(&x, &q, &s1, p); fp_sub(&x, &s0, &x, p); s0 = s1; s1 = x;
        fp_mul(&x, &q, &t1, p); fp_sub(&x, &t0, &x, p); t0 = t1; t1 = x;
    }
    if (r0.deg >= 0) { uint32_t li = inv[r0.c[r0.deg]]; fp_scale(&r0, li, p); fp_scale(&s0, li, p); fp_scale(&t0, li, p); }
    *g = r0; *s = s0; *t = t0;
}
/* divisor classes of s^2 = g(t), deg g = 5, in Mumford form (u monic of degree <= 2, deg v < deg u, u | v^2 - g); the
 * class 0 is (1, 0).  Cantor's composition and reduction. */
typedef struct { fpoly u, v; } jdiv;
static void jac_add(jdiv *R, const jdiv *A, const jdiv *B, const fpoly *g, uint32_t p, const uint32_t *inv) {
    fpoly d1, e1, e2, d, c1, c2, vs, u, v, t1, t2, t3, s1, s2;
    fp_xgcd(&d1, &e1, &e2, &A->u, &B->u, p, inv);
    fp_add(&vs, &A->v, &B->v, p);
    fp_xgcd(&d, &c1, &c2, &d1, &vs, p, inv);
    fp_mul(&s1, &c1, &e1, p); fp_mul(&s2, &c1, &e2, p);
    fp_mul(&u, &A->u, &B->u, p);
    if (d.deg > 0) { fp_divrem(&u, NULL, &u, &d, p, inv); fp_divrem(&u, NULL, &u, &d, p, inv); }
    fp_mul(&t1, &s1, &A->u, p); fp_mul(&t1, &t1, &B->v, p);
    fp_mul(&t2, &s2, &B->u, p); fp_mul(&t2, &t2, &A->v, p);
    fp_mul(&t3, &A->v, &B->v, p); fp_add(&t3, &t3, g, p); fp_mul(&t3, &c2, &t3, p);
    fp_add(&v, &t1, &t2, p); fp_add(&v, &v, &t3, p);
    if (d.deg > 0) fp_divrem(&v, NULL, &v, &d, p, inv);
    fp_divrem(NULL, &v, &v, &u, p, inv);
    while (u.deg > 2) {
        fpoly w, nv; fp_mul(&w, &v, &v, p); fp_sub(&w, g, &w, p); fp_divrem(&u, NULL, &w, &u, p, inv);
        fp_scale(&u, inv[u.c[u.deg]], p);
        nv.deg = -1; fp_sub(&nv, &nv, &v, p); fp_divrem(NULL, &v, &nv, &u, p, inv);
    }
    R->u = u; R->v = v;
}

#define IXMAX 40                 /* primes used at most */
typedef struct {
    uint32_t p, r;               /* the prime, the root of f */
    uint32_t M, m2, t, N;        /* H = { x1 G1 + x2 G2 }: G1 of order M, m2 G2 = t G1, N = M m2 classes, numbered x1 + M x2 */
    uint32_t g[6];               /* the quintic */
    uint32_t *inv;               /* inverses mod p */
    uint16_t *sq;                /* sq[v]: the smaller square root of v mod p, 0xFFFF if v is not a square */
    uint32_t *dlpt;              /* dlpt[2 t + e]: number of the class of the point (t, s), s = sq[g(t)] for e = 0 and -s for e = 1 */
    uint8_t *W; uint32_t stride; /* W[x2 * stride + x1], x1 in [0, 2M): the class (x1 mod M, x2) lies in W */
    double dens;                 /* |W| / N */
    uint32_t dinf[2];            /* classes of the points at infinity (1 : s : 0), (1 : -s : 0), s^2 = f6 (if f6 is a square) */
    uint32_t *d;                 /* d[i]: class of the i-th rational point, with the sign of the stored y */
    uint8_t *rows;               /* row c (64 nw bytes): bit 2 i + e = [c + d_i in W] (e = 0), [c - d_i in W] (e = 1) */
} ixp_t;
static ixp_t IX[IXMAX]; static int nix = 0;
static inline uint32_t ixg_add(const ixp_t *X, uint32_t c, uint32_t d) {
    uint32_t M = X->M, s1 = c % M + d % M, s2 = c / M + d / M;
    if (s2 >= X->m2) { s2 -= X->m2; s1 += X->t; }
    return s1 % M + M * s2;
}
static inline uint32_t ixg_neg(const ixp_t *X, uint32_t c) {
    uint32_t M = X->M, x1 = c % M, x2 = c / M;
    return x2 == 0 ? (M - x1) % M : (2 * M - x1 - X->t) % M + M * (X->m2 - x2);
}
/* a class in Mumford form: u = t^2 + u1 t + u0 (deg 2), t + u0 (deg 1), 1 (deg 0); v = v1 t + v0 */
typedef struct { int deg; uint32_t u1, u0, v1, v0; } jst_t;
typedef struct { uint32_t p, a, b, g5, g4, g3, ig5; const uint32_t *inv; const fpoly *g; jdiv Pd; } jpt_t;   /* a point (a, b) to add, with the constants of the curve */
static void jpt_set(jpt_t *Q, uint32_t p, const uint32_t *inv, const fpoly *g, uint32_t a, uint32_t b) {
    Q->p = p; Q->a = a; Q->b = b; Q->inv = inv; Q->g = g; Q->g5 = g->c[5]; Q->g4 = g->c[4]; Q->g3 = g->c[3]; Q->ig5 = inv[g->c[5]];
    Q->Pd.u.deg = 1; Q->Pd.u.c[0] = (p - a) % p; Q->Pd.u.c[1] = 1; Q->Pd.v.deg = b ? 0 : -1; Q->Pd.v.c[0] = b;
}
/* S += [Q - inf'] */
static inline void jst_addpt(jst_t *S, const jpt_t *Q) {
    uint32_t p = Q->p, a = Q->a, u1 = S->u1, u0 = S->u0, v1 = S->v1, v0 = S->v0;
    uint32_t ua = S->deg == 2 ? (a * a + u1 * a + u0) % p : 0;
    if (ua) {
        /* u (t - a) has degree 3 and v' = v + kk u takes the value b at a; one reduction step: u'' = (g - v'^2) / (u (t - a)),
         * made monic, v'' = -v' mod u'' */
        uint32_t va = (v1 * a + v0) % p, kk = (Q->b + p - va) * Q->inv[ua] % p;
        uint32_t V1 = (v1 + kk * u1) % p, V0 = (v0 + kk * u0) % p;
        uint32_t c2 = (u1 + p - a) % p, c1 = (u0 + p - a * u1 % p) % p;
        uint32_t w4 = (Q->g4 + p - kk * kk % p) % p, w3 = (Q->g3 + 2 * (p - kk * V1 % p)) % p;
        uint32_t q1 = (w4 + p - Q->g5 * c2 % p) % p, q0 = (w3 + 2 * p - Q->g5 * c1 % p - q1 * c2 % p) % p;
        uint32_t U1 = q1 * Q->ig5 % p, U0 = q0 * Q->ig5 % p;
        S->v1 = (kk * U1 + p - V1) % p; S->v0 = (kk * U0 + p - V0) % p; S->u1 = U1; S->u0 = U0;
    } else {
        jdiv D, R; D.u.deg = S->deg; if (S->deg == 2) { D.u.c[0] = u0; D.u.c[1] = u1; D.u.c[2] = 1; } else if (S->deg == 1) { D.u.c[0] = u0; D.u.c[1] = 1; } else D.u.c[0] = 1;
        D.v.deg = 1; D.v.c[0] = v0; D.v.c[1] = v1; fp_trim(&D.v);
        jac_add(&R, &D, &Q->Pd, Q->g, p, Q->inv);
        S->deg = R.u.deg; S->u1 = R.u.deg == 2 ? R.u.c[1] : 0; S->u0 = R.u.deg >= 1 ? R.u.c[0] : 1;
        S->v1 = R.v.deg >= 1 ? R.v.c[1] : 0; S->v0 = R.v.deg >= 0 ? R.v.c[0] : 0;
    }
}
static inline uint64_t jst_key(const jst_t *S, uint32_t p) { return ((((uint64_t)S->deg * p + S->u1) * p + S->u0) * p + S->v1) * p + S->v0 + 1; }
static inline uint64_t ix_mix(uint64_t k) { k ^= k >> 33; k *= UINT64_C(0xff51afd7ed558ccd); k ^= k >> 33; k *= UINT64_C(0xc4ceb9fe1a85ec53); k ^= k >> 33; return k; }

/* the structures of one prime; returns 0 if the prime cannot be used */
static int ix_prime(ixp_t *X, uint32_t p) {
    uint32_t fm[7]; for (int d = 0; d < 7; d++) fm[d] = mpz_fdiv_ui(f[d], p);
    if (fm[6] == 0) return 0;
    uint32_t *inv = malloc(p * sizeof(uint32_t)); inv[0] = 0; inv[1] = 1; for (uint32_t a = 2; a < p; a++) inv[a] = (p - (p / a) * inv[p % a] % p) % p;
    fpoly F, dF, gg, s_, t_; F.deg = 6; for (int d = 0; d < 7; d++) F.c[d] = fm[d];
    dF.deg = 5; for (int d = 0; d < 6; d++) dF.c[d] = (d + 1) * fm[d + 1] % p; fp_trim(&dF);
    if (dF.deg < 0) { free(inv); return 0; }
    fp_xgcd(&gg, &s_, &t_, &F, &dF, p, inv);
    if (gg.deg != 0) { free(inv); return 0; }                         /* bad reduction */
    int r = -1; for (uint32_t x = 0; x < p && r < 0; x++) { uint32_t v = 0; for (int d = 6; d >= 0; d--) v = (v * x + fm[d]) % p; if (v == 0) r = (int)x; }
    if (r < 0) { free(inv); return 0; }                               /* no root */
    uint32_t tc[7]; memcpy(tc, fm, sizeof tc);                        /* Taylor coefficients of f at r */
    for (int i = 0; i < 6; i++) for (int j = 5; j >= i; j--) tc[j] = (tc[j] + (uint32_t)r * tc[j + 1]) % p;
    fpoly g; g.deg = 5; for (int j = 0; j <= 5; j++) g.c[j] = tc[6 - j];   /* g(t) = t^6 f(r + 1/t) */
    X->p = p; X->r = (uint32_t)r; X->inv = inv; for (int j = 0; j <= 5; j++) X->g[j] = g.c[j];
    X->sq = malloc(p * sizeof(uint16_t)); for (uint32_t v = 0; v < p; v++) X->sq[v] = 0xFFFF;
    for (uint32_t s = 0; s <= p / 2; s++) X->sq[s * s % p] = (uint16_t)s;
    X->dlpt = malloc(2 * p * sizeof(uint32_t));
    /* the points */
    uint32_t *pt = malloc(2 * p * sizeof(uint32_t)), *ps = malloc(2 * p * sizeof(uint32_t)); int npt = 0;
    for (uint32_t t = 0; t < p; t++) {
        uint32_t v = 0; for (int d = 5; d >= 0; d--) v = (v * t + g.c[d]) % p;
        if (X->sq[v] == 0xFFFF) continue;
        pt[npt] = t; ps[npt] = X->sq[v]; npt++;
        if (v) { pt[npt] = t; ps[npt] = p - X->sq[v]; npt++; }
    }
    double sp = 1.0 + __builtin_sqrt((double)p); uint32_t maxord = (uint32_t)(sp * sp * sp * sp) + 8;
    uint32_t hsz = 64; while (hsz < 2 * maxord) hsz *= 2;
    uint64_t *hk = malloc(hsz * sizeof(uint64_t)); uint32_t *hv = malloc(hsz * sizeof(uint32_t));
    int ok = 0; uint64_t seed = 0x9E3779B97F4A7C15ULL * p;
    #define IX_RECORD(S, num) do { if ((S).deg == 1) { uint32_t t0_ = (p - (S).u0) % p; if ((S).v0 == 0) { X->dlpt[2 * t0_] = X->dlpt[2 * t0_ + 1] = (num); } else X->dlpt[2 * t0_ + ((S).v0 > p / 2)] = (num); cnt++; } } while (0)
    for (int attempt = 0; attempt < 6 && !ok && npt >= 2; attempt++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL; int i1 = (int)((seed >> 33) % (uint64_t)npt);
        if (ps[i1] == 0) continue;
        /* the multiples of the class G1 of a point, until the class 0 comes back; the classes of degree 1 met on the way
         * are the points */
        jpt_t Q1; jpt_set(&Q1, p, inv, &g, pt[i1], ps[i1]);
        for (uint32_t i = 0; i < 2 * p; i++) X->dlpt[i] = UINT32_MAX;
        memset(hk, 0, hsz * sizeof(uint64_t));
        jst_t S; S.deg = 1; S.u1 = 0; S.u0 = (p - pt[i1]) % p; S.v1 = 0; S.v0 = ps[i1];
        uint32_t k = 1; int cnt = 0;
        while (S.deg > 0 && k <= maxord) {
            IX_RECORD(S, k);
            { uint64_t key = jst_key(&S, p); uint32_t h = (uint32_t)ix_mix(key) & (hsz - 1); while (hk[h]) h = (h + 1) & (hsz - 1); hk[h] = key; hv[h] = k; }
            jst_addpt(&S, &Q1); k++;
        }
        if (S.deg != 0) continue;
        uint32_t M = k, m2 = 1, t = 0;
        if (cnt < npt) {
            /* a second generator: a point G2 outside <G1>; the cosets x2 G2 + <G1> until x2 G2 lies in <G1> */
            int i2 = -1; for (int i = 0; i < npt && i2 < 0; i++) if (ps[i] && X->dlpt[2 * pt[i] + (ps[i] > p / 2)] == UINT32_MAX) i2 = i;
            if (i2 < 0) continue;
            jpt_t Q2; jpt_set(&Q2, p, inv, &g, pt[i2], ps[i2]);
            jst_t T; T.deg = 1; T.u1 = 0; T.u0 = (p - pt[i2]) % p; T.v1 = 0; T.v0 = ps[i2];
            int found = 0;
            for (uint32_t x2 = 1; (uint64_t)M * x2 <= maxord && !found; x2++) {
                if (T.deg == 0) { m2 = x2; t = 0; found = 1; break; }
                uint64_t key = jst_key(&T, p); uint32_t h = (uint32_t)ix_mix(key) & (hsz - 1); int in = 0;
                while (hk[h]) { if (hk[h] == key) { in = 1; break; } h = (h + 1) & (hsz - 1); }
                if (in) { m2 = x2; t = hv[h]; found = 1; break; }
                S = T; for (uint32_t x1 = 0; x1 < M; x1++) { IX_RECORD(S, x1 + M * x2); jst_addpt(&S, &Q1); }
                jst_addpt(&T, &Q2);
            }
            if (!found || cnt < npt) continue;
        }
        ok = 1; X->M = M; X->m2 = m2; X->t = t; X->N = M * m2;
    }
    #undef IX_RECORD
    free(hk); free(hv);
    if (!ok) { free(pt); free(ps); free(X->sq); free(X->dlpt); free(inv); return 0; }
    uint32_t M = X->M, m2 = X->m2, N = X->N;
    X->stride = 2 * M + 8; X->W = calloc((size_t)m2 * X->stride, 1);
    uint32_t *dl = malloc((npt + 1) * sizeof(uint32_t)); dl[0] = 0;
    for (int i = 0; i < npt; i++) dl[i + 1] = X->dlpt[2 * pt[i] + (ps[i] > p / 2)];
    for (int i = 0; i <= npt; i++) for (int j = i; j <= npt; j++) { uint32_t w = ixg_add(X, dl[i], dl[j]); X->W[(w / M) * X->stride + w % M] = 1; }
    uint32_t nw = 0; for (uint32_t x2 = 0; x2 < m2; x2++) for (uint32_t x1 = 0; x1 < M; x1++) { uint8_t b = X->W[x2 * X->stride + x1]; nw += b; X->W[x2 * X->stride + M + x1] = b; }
    X->dens = (double)nw / N; X->d = NULL; X->rows = NULL; X->dinf[0] = X->dinf[1] = 0;
    if (has_sqrt_f6) { uint32_t sv = mpz_fdiv_ui(sqrt_f6z, p), e = sv > p / 2; X->dinf[0] = X->dlpt[e]; X->dinf[1] = X->dlpt[1 - e]; }   /* t = 0, s = +-sqrt(f6) */
    free(pt); free(ps); free(dl);
    return 1;
}
/* class of the rational point i modulo the prime of X */
static uint32_t ix_dl(const ixp_t *X, int i) {
    uint32_t p = X->p, a = mpz_fdiv_ui(pts[i].a, p), b = mpz_fdiv_ui(pts[i].b, p), c = mpz_fdiv_ui(pts[i].c, p);
    uint32_t den = (a + p - X->r * c % p) % p;
    if (den == 0) return 0;                                            /* the point above the root: inf' */
    uint32_t di = X->inv[den], t = c * di % p, s = b * di % p * di % p * di % p;
    uint32_t v = 0; for (int d = 5; d >= 0; d--) v = (v * t + X->g[d]) % p;
    if (s * s % p != v) { fprintf(stderr, "index sieve: point %d is not on the curve modulo %u\n", i, p); abort(); }
    uint32_t dl = X->dlpt[2 * t + (s > p / 2)];
    if (dl == UINT32_MAX) { fprintf(stderr, "index sieve: point %d has no class modulo %u\n", i, p); abort(); }
    return dl;
}
/* the primes (once): in increasing order, until the product of the densities is below target */
static void ix_init(int pmax, double target, int verbose) {
    double prod = 1;
    for (uint32_t p = 11; p <= (uint32_t)pmax && nix < IXMAX && prod > target; p += 2) {
        int pr = 1; for (uint32_t q = 3; q * q <= p; q += 2) if (p % q == 0) pr = 0;
        if (!pr) continue;
        if (ix_prime(&IX[nix], p)) { prod *= IX[nix].dens; nix++; }
    }
    if (verbose) { fprintf(stderr, "index sieve: %d primes:", nix); for (int k = 0; k < nix; k++) { if (IX[k].m2 > 1) fprintf(stderr, " %u(%ux%u,%.2f)", IX[k].p, IX[k].M, IX[k].m2, IX[k].dens); else fprintf(stderr, " %u(%u,%.2f)", IX[k].p, IX[k].M, IX[k].dens); } fprintf(stderr, "; product of the densities %.3g\n", prod); }
}
/* the classes of the points and the rows, for n points (every pass) */
static int ix_nw = 0, ix_sh = 0;   /* rows are 64 ix_nw bytes = 1 << ix_sh */
static int IXnb = 0;               /* 16 primes per vector */
static __m512i *IXD1p = NULL, *IXD2p = NULL, *IXD1m = NULL, *IXD2m = NULL;   /* (x1, x2) of the classes of the points and of their negatives */
static __m512i IXMv[(IXMAX + 15) / 16], IXm2v[(IXMAX + 15) / 16], IXtv[(IXMAX + 15) / 16], IXI1[2][(IXMAX + 15) / 16], IXI2[2][(IXMAX + 15) / 16];
static const uint8_t *IXrows[IXMAX]; static uint8_t *IXarena = NULL; static size_t IXarena_sz = 0;
static void ix_rows(int n) {
    ix_nw = 1; ix_sh = 6; while (512 * ix_nw < 2 * n) { ix_nw *= 2; ix_sh++; }
    size_t rb = (size_t)64 * ix_nw; int ncol = 512 * ix_nw;
    /* one arena for the rows of all primes, on huge pages when the system gives them */
    size_t need = 0, offs[IXMAX]; for (int k = 0; k < nix; k++) { offs[k] = need; need += ((size_t)IX[k].N + 4) * rb; }
    need = (need + ((size_t)1 << 21) - 1) & ~(((size_t)1 << 21) - 1);
    if (need > IXarena_sz) { free(IXarena); IXarena = aligned_alloc((size_t)1 << 21, need); IXarena_sz = need; madvise(IXarena, need, MADV_HUGEPAGE); }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < nix; k++) {
        ixp_t *X = &IX[k]; uint32_t M = X->M, m2 = X->m2;
        X->d = realloc(X->d, n * sizeof(uint32_t)); for (int i = 0; i < n; i++) X->d[i] = ix_dl(X, i);
        X->rows = IXarena + offs[k]; IXrows[k] = X->rows;
        /* column j = 2 i + e of the row of the class (x1, x2) is W at the class (x1, x2) + y_j, y_j the class of the point i or
         * its negative: for a fixed x2 that is W[x2' * stride + r_j + x1], so one gather reads it for 16 columns and 4 rows */
        uint32_t *y = malloc(ncol * sizeof(uint32_t)), *delta = aligned_alloc(64, ncol * sizeof(uint32_t));
        for (int j = 0; j < ncol; j++) { int i = j >> 1; y[j] = i < n ? ((j & 1) ? ixg_neg(X, X->d[i]) : X->d[i]) : 0; }
        __mmask16 *valid = malloc((ncol / 16) * sizeof(__mmask16));
        for (int q = 0; q < ncol / 16; q++) { unsigned m = 0; for (int j = 0; j < 16; j++) if ((16 * q + j) >> 1 < n) m |= 1u << j; valid[q] = (__mmask16)m; }
        __m512i b0 = _mm512_set1_epi32(1), b1 = _mm512_set1_epi32(1 << 8), b2 = _mm512_set1_epi32(1 << 16), b3 = _mm512_set1_epi32(1 << 24);
        for (uint32_t x2 = 0; x2 < m2; x2++) {
            for (int j = 0; j < ncol; j++) { uint32_t s2 = x2 + y[j] / M, r1 = y[j] % M; if (s2 >= m2) { s2 -= m2; r1 = (r1 + X->t) % M; } delta[j] = s2 * X->stride + r1; }
            for (uint32_t x1 = 0; x1 < M; x1 += 4) {
                uint16_t *r0 = (uint16_t*)(X->rows + ((size_t)x2 * M + x1) * rb), *r1 = (uint16_t*)((uint8_t*)r0 + rb), *r2 = (uint16_t*)((uint8_t*)r1 + rb), *r3 = (uint16_t*)((uint8_t*)r2 + rb);
                for (int q = 0; q < ncol / 16; q++) {
                    __m512i x = _mm512_i32gather_epi32(_mm512_load_si512(delta + 16 * q), X->W + x1, 1);
                    r0[q] = _mm512_test_epi32_mask(x, b0) & valid[q]; r1[q] = _mm512_test_epi32_mask(x, b1) & valid[q];
                    r2[q] = _mm512_test_epi32_mask(x, b2) & valid[q]; r3[q] = _mm512_test_epi32_mask(x, b3) & valid[q];
                }
            }
        }
        free(y); free(delta); free(valid);
    }
    int nb = IXnb = (nix + 15) / 16;
    free(IXD1p); free(IXD2p); free(IXD1m); free(IXD2m);
    IXD1p = aligned_alloc(64, (size_t)n * nb * 64); IXD2p = aligned_alloc(64, (size_t)n * nb * 64); IXD1m = aligned_alloc(64, (size_t)n * nb * 64); IXD2m = aligned_alloc(64, (size_t)n * nb * 64);
    for (int b = 0; b < nb; b++) {
        uint32_t m[16], mm[16], tt[16], i1[2][16], i2[2][16];
        for (int j = 0; j < 16; j++) { int k = 16 * b + j; const ixp_t *X = &IX[k < nix ? k : 0]; int u = k < nix;
            m[j] = u ? X->M : 1; mm[j] = u ? X->m2 : 1; tt[j] = u ? X->t : 0;
            for (int si = 0; si < 2; si++) { i1[si][j] = u ? X->dinf[si] % X->M : 0; i2[si][j] = u ? X->dinf[si] / X->M : 0; } }
        IXMv[b] = _mm512_loadu_si512(m); IXm2v[b] = _mm512_loadu_si512(mm); IXtv[b] = _mm512_loadu_si512(tt);
        for (int si = 0; si < 2; si++) { IXI1[si][b] = _mm512_loadu_si512(i1[si]); IXI2[si][b] = _mm512_loadu_si512(i2[si]); }
    }
    for (int i = 0; i < n; i++) for (int b = 0; b < nb; b++) {
        uint32_t p1[16], p2[16], m1[16], m2[16];
        for (int j = 0; j < 16; j++) { int k = 16 * b + j; if (k < nix) { uint32_t d = IX[k].d[i], e = ixg_neg(&IX[k], d), M = IX[k].M; p1[j] = d % M; p2[j] = d / M; m1[j] = e % M; m2[j] = e / M; } else p1[j] = p2[j] = m1[j] = m2[j] = 0; }
        IXD1p[(size_t)i * nb + b] = _mm512_loadu_si512(p1); IXD2p[(size_t)i * nb + b] = _mm512_loadu_si512(p2); IXD1m[(size_t)i * nb + b] = _mm512_loadu_si512(m1); IXD2m[(size_t)i * nb + b] = _mm512_loadu_si512(m2);
    }
}
/* (r1, r2) = (a1, a2) + (b1, b2) in the groups of the primes of block b */
static inline void ix_vadd(int b, __m512i a1, __m512i a2, __m512i b1, __m512i b2, __m512i *r1, __m512i *r2) {
    __m512i s2 = _mm512_add_epi32(a2, b2); __mmask16 cy = _mm512_cmpge_epu32_mask(s2, IXm2v[b]); s2 = _mm512_mask_sub_epi32(s2, cy, s2, IXm2v[b]);
    __m512i s1 = _mm512_add_epi32(a1, b1); s1 = _mm512_min_epu32(s1, _mm512_sub_epi32(s1, IXMv[b]));
    s1 = _mm512_mask_add_epi32(s1, cy, s1, IXtv[b]); s1 = _mm512_min_epu32(s1, _mm512_sub_epi32(s1, IXMv[b]));
    *r1 = s1; *r2 = s2;
}
/* the class of a point with sign e (0: as stored, 1: negative) */
#define IXPT1(i, e, b) ((e) ? IXD1m[(size_t)(i) * IXnb + (b)] : IXD1p[(size_t)(i) * IXnb + (b)])
#define IXPT2(i, e, b) ((e) ? IXD2m[(size_t)(i) * IXnb + (b)] : IXD2p[(size_t)(i) * IXnb + (b)])
/* byte offsets of the rows of the classes (x1, x2) */
static inline __m512i ix_off(int b, __m512i x1, __m512i x2, __m128i vsh) { return _mm512_sll_epi32(_mm512_add_epi32(x1, _mm512_mullo_epi32(x2, IXMv[b])), vsh); }
/* the AND over the primes of word w (512 bits) of the rows at the offsets o, under a mask; returns 0 if it is empty */
static inline int ix_chain(const uint32_t *o, int w, const uint64_t *msk, uint64_t *wd) {
    int nq = nix;
    __m512i acc = _mm512_and_si512(_mm512_load_si512(IXrows[0] + o[0] + 64 * w), _mm512_loadu_si512(msk));
    int k = 1;
    for (; k + 3 < nq; k += 4) {
        acc = _mm512_and_si512(_mm512_and_si512(acc, _mm512_and_si512(_mm512_load_si512(IXrows[k] + o[k] + 64 * w), _mm512_load_si512(IXrows[k+1] + o[k+1] + 64 * w))),
                               _mm512_and_si512(_mm512_load_si512(IXrows[k+2] + o[k+2] + 64 * w), _mm512_load_si512(IXrows[k+3] + o[k+3] + 64 * w)));
        if (_mm512_test_epi64_mask(acc, acc) == 0) return 0;
    }
    for (; k < nq; k++) acc = _mm512_and_si512(acc, _mm512_load_si512(IXrows[k] + o[k] + 64 * w));
    if (_mm512_test_epi64_mask(acc, acc) == 0) return 0;
    _mm512_store_si512(wd, acc); return 1;
}
/* mask of word w: the bits of the points i >= first, without the point excl (if excl >= 0) */
static inline void ix_mask(uint64_t *msk, int w, int first, int excl) {
    int lo = 2 * first - 512 * w;
    for (int q = 0; q < 8; q++) msk[q] = lo >= 64 * (q + 1) ? 0 : lo <= 64 * q ? ~UINT64_C(0) : ~UINT64_C(0) << (lo - 64 * q);
    if (excl >= 0) { int b = 2 * excl - 512 * w; if (b >= 0 && b < 512) msk[b >> 6] &= ~(UINT64_C(3) << (b & 63)); }
}
/* scalar test of one relation: kinf times the point at infinity si, plus the points pt[j] with multiplicities mult[j] and the
 * signs given by the bits of eps */
static inline int ix_test(int nterm, const int *pt, const int *mult, int eps, int kinf, int si) {
    for (int k = 0; k < nix; k++) {
        const ixp_t *X = &IX[k]; uint32_t w = 0;
        for (int q = 0; q < kinf; q++) w = ixg_add(X, w, X->dinf[si]);
        for (int j = 0; j < nterm; j++) { uint32_t d = X->d[pt[j]]; if ((eps >> j) & 1) d = ixg_neg(X, d); for (int q = 0; q < mult[j]; q++) w = ixg_add(X, w, d); }
        if (!X->W[(w / X->M) * X->stride + w % X->M]) return 0;
    }
    return 1;
}

/* set of 4-subsets with signs whose residual points are known (insert-only, shared by the threads) */
static uint64_t *KQ = NULL; static uint64_t KQmask = 0; static long KQn = 0, KQmax = 0;
static inline uint64_t kq_key(int t0, int t1, int t2, int t3, int eps) { return (UINT64_C(1) << 63) | ((uint64_t)t0 << 45) | ((uint64_t)t1 << 31) | ((uint64_t)t2 << 17) | ((uint64_t)t3 << 3) | (uint64_t)eps; }
static inline int kq_has(uint64_t k) { uint64_t h = ix_mix(k) & KQmask; for (;;) { uint64_t v = __atomic_load_n(&KQ[h], __ATOMIC_RELAXED); if (v == k) return 1; if (v == 0) return 0; h = (h + 1) & KQmask; } }
static inline void kq_add(uint64_t k) {
    if (__atomic_load_n(&KQn, __ATOMIC_RELAXED) >= KQmax) return;
    uint64_t h = ix_mix(k) & KQmask;
    for (;;) {
        uint64_t v = __atomic_load_n(&KQ[h], __ATOMIC_RELAXED);
        if (v == k) return;
        if (v == 0) { uint64_t e = 0; if (__atomic_compare_exchange_n(&KQ[h], &e, k, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) { __atomic_fetch_add(&KQn, 1, __ATOMIC_RELAXED); return; } if (e == k) return; }
        h = (h + 1) & KQmask;
    }
}
/* index of a known x-value (hash on the low limbs, exact comparison), -1 if it is not known */
static int *KX = NULL; static uint32_t KXmask = 0;
static inline uint32_t kx_hash(mpz_t a, mpz_t c) { uint64_t h = mpz_get_ui(a) * UINT64_C(0x9E3779B97F4A7C15) ^ (mpz_get_ui(c) + (mpz_sgn(a) < 0 ? UINT64_C(0x632BE59BD9B4E019) : 0)) * UINT64_C(0xC2B2AE3D27D4EB4F); return (uint32_t)(h >> 32); }
static void kx_build(int n) {
    uint32_t sz = 1024; while (sz < 4u * (uint32_t)n) sz *= 2;
    free(KX); KX = malloc(sz * sizeof(int)); for (uint32_t i = 0; i < sz; i++) KX[i] = -1; KXmask = sz - 1;
    for (int i = 0; i < n; i++) { uint32_t h = kx_hash(pts[i].a, pts[i].c) & KXmask; while (KX[h] >= 0) h = (h + 1) & KXmask; KX[h] = i; }
}
static inline int kx_find(mpz_t a, mpz_t c) {
    uint32_t h = kx_hash(a, c) & KXmask;
    for (;;) { int i = KX[h]; if (i < 0) return -1; if (mpz_cmp(pts[i].a, a) == 0 && mpz_cmp(pts[i].c, c) == 0) return i; h = (h + 1) & KXmask; }
}

static int ix_pretest = 0;   /* number of table primes tried on a survivor before the exact check */
/* a survivor of the sieve, given by four points and their signs (1 = negative): exact check unless its residual points are
 * already known; new points go to loc */
static inline void ix_survivor(const int *pi, const int *pe, newpt_t *loc, int *nloc, uint64_t *nexact, uint64_t *nskip) {
    int idx[4] = {pi[0], pi[1], pi[2], pi[3]}, sg[4] = {pe[0], pe[1], pe[2], pe[3]};
    for (int a = 1; a < 4; a++) { int ti = idx[a], ts = sg[a], b = a - 1; while (b >= 0 && idx[b] > ti) { idx[b+1] = idx[b]; sg[b+1] = sg[b]; b--; } idx[b+1] = ti; sg[b+1] = ts; }
    int eps = (sg[0] != sg[3] ? 1 : 0) | (sg[1] != sg[3] ? 2 : 0) | (sg[2] != sg[3] ? 4 : 0);   /* the exact check has y_3 positive */
    if (kq_has(kq_key(idx[0], idx[1], idx[2], idx[3], eps))) { (*nskip)++; return; }
    for (int k = 0; k < ix_pretest; k++) { sub_t S; uint64_t disc; if (prep_subset(k, idx, &S) && !test_pattern(k, &S, eps, &disc)) return; }   /* a weak sieve leaves many survivors: a few of the table primes first */
    newpt_t out[2]; int nn = exact_check(idx, eps, out); (*nexact)++;
    int who[2] = {-1, -1};
    for (int u = 0; u < nn; u++) who[u] = kx_find(out[u].a, out[u].c);
    if (nn == 2 && who[0] >= 0 && who[1] >= 0 && who[0] != who[1] && out[0].ys && out[1].ys) {
        /* six known points with distinct x and y != 0 on one cubic: every 4-subset of them has known residual points */
        int id[6] = {idx[0], idx[1], idx[2], idx[3], who[0], who[1]}, s6[6] = {eps & 1 ? -1 : 1, eps & 2 ? -1 : 1, eps & 4 ? -1 : 1, 1, out[0].ys, out[1].ys};
        int distinct = 1; for (int a = 0; a < 4; a++) if (id[a] == who[0] || id[a] == who[1] || is_weier(id[a])) distinct = 0;
        if (distinct) {
            for (int a = 1; a < 6; a++) { int ti = id[a], ts = s6[a], b = a - 1; while (b >= 0 && id[b] > ti) { id[b+1] = id[b]; s6[b+1] = s6[b]; b--; } id[b+1] = ti; s6[b+1] = ts; }
            for (int a = 0; a < 6; a++) for (int b = a + 1; b < 6; b++) {   /* leave out a and b */
                int t[4], s[4], m = 0; for (int q = 0; q < 6; q++) if (q != a && q != b) { t[m] = id[q]; s[m] = s6[q]; m++; }
                kq_add(kq_key(t[0], t[1], t[2], t[3], (s[0] != s[3] ? 1 : 0) | (s[1] != s[3] ? 2 : 0) | (s[2] != s[3] ? 4 : 0)));
            }
        }
    }
    for (int u = 0; u < nn; u++) {
        if (who[u] < 0 && *nloc < 64) loc[(*nloc)++] = out[u];
        else { mpz_clear(out[u].a); mpz_clear(out[u].b); mpz_clear(out[u].c); }
    }
}

/* the four-point pass with the index sieve: all 4-subsets with at least one point of index >= nstart, all signs.  The
 * points are taken in the order "new points first": a subset whose first point in that order is new is met once. */
static void four_pass_ix(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *nsub_out, uint64_t *nsurv_out, uint64_t *nexact_out) {
    double tr1 = omp_get_wtime();
    { uint64_t sz = UINT64_C(1) << 22; while (sz < (uint64_t)n * n * 24 && sz < (UINT64_C(1) << 27)) sz *= 2;
      free(KQ); KQ = calloc(sz, sizeof(uint64_t)); KQmask = sz - 1; KQn = 0; KQmax = (long)(sz * 6 / 10); }
    int nb = IXnb, nw = ix_nw, npf = nix < 16 ? nix : 16, nnew = n - nstart; __m128i vsh = _mm_cvtsi32_si128(ix_sh);
    int *ord = malloc(n * sizeof(int)); for (int a = 0; a < n; a++) ord[a] = a < nnew ? nstart + a : a - nnew;
    /* after[c]: the bits of the points that come after position c in the order */
    uint64_t *after = aligned_alloc(64, (size_t)n * 64 * nw); memset(after, 0, (size_t)n * 64 * nw);
    for (int c = n - 2; c >= 0; c--) { uint64_t *m = after + (size_t)c * 8 * nw; memcpy(m, after + (size_t)(c + 1) * 8 * nw, 64 * nw); int i = ord[c + 1]; m[(2 * i) >> 6] |= UINT64_C(3) << ((2 * i) & 63); }
    uint64_t nsub = 0, nsurv = 0, nexact = 0, nskip = 0;
    #pragma omp parallel reduction(+:nsub, nsurv, nexact, nskip)
    {
        newpt_t loc[64]; int nloc = 0;
        uint32_t offb[2][4][64] __attribute__((aligned(64)));
        #pragma omp for schedule(dynamic, 1)
        for (int a = 0; a < nnew; a++) {
            int i0 = ord[a];
            for (int b = a + 1; b < n - 2; b++) {
                int i1 = ord[b];
                __m512i s1[2][(IXMAX + 15) / 16], s2[2][(IXMAX + 15) / 16];
                for (int q = 0; q < nb; q++) for (int e1 = 0; e1 < 2; e1++) ix_vadd(q, IXPT1(i0, 0, q), IXPT2(i0, 0, q), IXPT1(i1, e1, q), IXPT2(i1, e1, q), &s1[e1][q], &s2[e1][q]);
                /* the offsets of a third point are computed, and their rows prefetched, one step before they are used */
                int prev = -1, cur = 0;
                for (int c = b + 1; c <= n - 1; c++) {
                    int have = c < n - 1;
                    if (have) {
                        int i2 = ord[c]; uint32_t (*off)[64] = offb[cur];
                        for (int q = 0; q < nb; q++) for (int e1 = 0; e1 < 2; e1++) for (int e2 = 0; e2 < 2; e2++) {
                            __m512i r1, r2; ix_vadd(q, s1[e1][q], s2[e1][q], IXPT1(i2, e2, q), IXPT2(i2, e2, q), &r1, &r2);
                            _mm512_store_si512(off[2 * e1 + e2] + 16 * q, ix_off(q, r1, r2, vsh));
                        }
                        for (int pat = 0; pat < 4; pat++) for (int k = 0; k < npf; k++) _mm_prefetch((const char*)(IXrows[k] + off[pat][k]), _MM_HINT_T0);
                    }
                    if (prev >= 0) {
                        int pc = prev, i2 = ord[pc]; uint32_t (*off)[64] = offb[cur ^ 1];
                        const uint64_t *am = after + (size_t)pc * 8 * nw;
                        nsub += (uint64_t)(n - 1 - pc);
                        for (int pat = 0; pat < 4; pat++) {                           /* pat = 2 e1 + e2: signs of the second and third points (the first positive) */
                            for (int w = 0; w < nw; w++) {
                                uint64_t wd[8] __attribute__((aligned(64)));
                                if (!ix_chain(off[pat], w, am + 8 * w, wd)) continue;
                                for (int q = 0; q < 8; q++) { uint64_t x = wd[q];
                                    while (x) {
                                        int bit = 512 * w + 64 * q + __builtin_ctzll(x); x &= x - 1;
                                        int pi[4] = {i0, i1, i2, bit >> 1}, pe[4] = {0, pat >> 1, pat & 1, bit & 1};
                                        nsurv++;
                                        ix_survivor(pi, pe, loc, &nloc, &nexact, &nskip);
                                        if (nloc >= 62) { FLUSH_LOC(); }
                                    }
                                }
                            }
                        }
                    }
                    if (have) { prev = c; cur ^= 1; } else prev = -1;
                }
            }
            if (nloc) { FLUSH_LOC(); }
        }
    }
    free(ord); free(after);
    if (getenv("CUBIC_IXSTAT")) fprintf(stderr, "index sieve: four-point pass %.3fs; %lu survivors, %lu exact checks, %lu skipped as known relations (%ld 4-subsets recorded)\n", omp_get_wtime() - tr1, nsurv, nexact, nskip, KQn);
    *nsub_out += nsub; *nsurv_out += nsurv; *nexact_out += nexact;
}

/* ---------------- tangent (Hermite) pass ----------------
 * pattern (2,1,1) (i0 tangent, i1 simple, lanes = i2 consecutive, 4 sign patterns for y1, y2) is vectorised:
 * Newton form p = y0 + y0' (x - x0) + c2 (x - x0)^2 + c3 (x - x0)^2 (x - x1) with the divided differences
 * c2 = [x0,x0,x1] = (([x0,x1] - y0') / (x1 - x0)), c3 = [x0,x0,x1,x2] = ([x0,x1,x2] - c2) / (x2 - x0) taken from the
 * inverse tables (no linear solve); p3 = c3, p2 = c2 - c3 (2x0 + x1), p1 = y0' - 2 x0 c2 + c3 (x0^2 + 2 x0 x1),
 * Q = (x-x0)^2 (x-x1) (x-x2): Q3 = -(2x0 + x1 + x2), Q2 = x0^2 + 2 x0 x1 + (2x0 + x1) x2.  Survivors of the vector
 * primes go through hermite_test_mod for the remaining primes, then hermite_exact.  The patterns (2,2), (3,1), (4)
 * (O(n^2) tests) stay scalar.  A prime is skipped (undecided) for a lane if some x-values coincide mod p or y0 = 0 mod p. */
static void hermite_pass(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *ntests) {
    static const int pats[4][4] = { {2,1,1,0}, {2,2,0,0}, {3,1,0,0}, {4,0,0,0} };
    static const int npat_k[4] = {3, 2, 2, 1};
    uint64_t cnt = 0;
    #pragma omp parallel reduction(+:cnt)
    {
        newpt_t loc[64]; int nloc = 0;
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 < n; i0++) {
            if (!is_weier(i0)) for (int i1 = 0; i1 < n; i1++) {   /* pattern (2,1,1) */
                if (i1 == i0) continue;
                int i2start = i1 + 1; if (i0 < nstart && i1 < nstart && i2start < nstart) i2start = nstart;
                if (i2start >= n) continue;
                struct { uint64_t A, NA, Bc, c2[2], d01[2], E[2], Y1[2]; int ok; } H[NVP];
                for (int k = 0; k < NVP; k++) {
                    uint64_t p = PR[k]; double pinv = PINV[k];
                    uint64_t x0 = xr[k][i0], y0 = yr[k][i0], x1 = xr[k][i1], y1 = yr[k][i1], inv10 = invd[k][(size_t)i1 * ncap + i0], yd0 = yd1[k][i0];
                    H[k].ok = inv10 != 0 && inv2y[k][i0] != 0;
                    H[k].A = (2 * x0 + x1) % p; H[k].NA = (p - H[k].A) % p;
                    H[k].Bc = redp(x0 * ((x0 + 2 * x1) % p), p, pinv);
                    for (int s1 = 0; s1 < 2; s1++) {
                        uint64_t Y1 = s1 ? (p - y1) % p : y1;
                        uint64_t d01 = redp(((Y1 + p - y0) % p) * inv10, p, pinv);
                        uint64_t c2 = redp(((d01 + p - yd0) % p) * inv10, p, pinv);
                        H[k].Y1[s1] = Y1; H[k].d01[s1] = d01; H[k].c2[s1] = c2;
                        H[k].E[s1] = (yd0 + p - redp(((2 * x0) % p) * c2, p, pinv)) % p;   /* y0' - 2 x0 c2 */
                    }
                }
                for (int i2 = i2start; i2 < n; i2 += 8) {
                    int nl = n - i2 < 8 ? n - i2 : 8;
                    __mmask8 lanes = (__mmask8)((1u << nl) - 1);
                    if (i0 >= i2 && i0 < i2 + nl) lanes &= (__mmask8)~(1u << (i0 - i2));
                    if (!lanes) continue;
                    __mmask8 alive[4]; for (int e = 0; e < 4; e++) alive[e] = lanes;
                    cnt += 4 * (uint64_t)__builtin_popcount(lanes);
                    for (int k = 0; k < NVP && (alive[0]|alive[1]|alive[2]|alive[3]); k++) {
                        if (!H[k].ok) continue;
                        const vprime_t *P = &VP[k];
                        __m512i X2 = _mm512_maskz_loadu_epi64(lanes, xr[k] + i2), Y2 = _mm512_maskz_loadu_epi64(lanes, yr[k] + i2);
                        __m512i ninv21 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i1 * ncap + i2);   /* 1/(x1 - x2) */
                        __m512i ninv20 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i0 * ncap + i2);   /* 1/(x0 - x2) */
                        __mmask8 ok = _mm512_test_epi64_mask(ninv21, ninv21) & _mm512_test_epi64_mask(ninv20, ninv20);
                        if (!ok) continue;
                        __m512i inv21 = _mm512_sub_epi64(P->vp, ninv21), inv20 = _mm512_sub_epi64(P->vp, ninv20);   /* in (0, p] */
                        __m512i vA = _mm512_set1_epi64(H[k].A), vNA = _mm512_set1_epi64(H[k].NA), vBc = _mm512_set1_epi64(H[k].Bc);
                        __m512i NQ3 = _mm512_add_epi64(vA, X2);                                     /* -Q3 = 2x0 + x1 + x2, < 2p */
                        __m512i NQ2 = _mm512_sub_epi64(P->v2p, vred(vmla(vBc, vA, X2), P));         /* -Q2 */
                        __m512i NY2 = _mm512_sub_epi64(P->v2p, Y2);
                        #pragma GCC unroll 4
                        for (int eps = 0; eps < 4; eps++) {   /* bit 0: sign of y1, bit 1: sign of y2 */
                            if (!alive[eps]) continue;
                            int s1 = eps & 1;
                            __m512i Y2s = (eps & 2) ? NY2 : Y2;
                            __m512i d12 = vmul(_mm512_add_epi64(Y2s, _mm512_set1_epi64(2 * P->p - H[k].Y1[s1])), inv21, P);      /* [x1,x2] */
                            __m512i d012 = vmul(_mm512_add_epi64(d12, _mm512_set1_epi64(2 * P->p - H[k].d01[s1])), inv20, P);   /* [x0,x1,x2] */
                            __m512i c3 = vmul(_mm512_add_epi64(d012, _mm512_set1_epi64(2 * P->p - H[k].c2[s1])), inv20, P);     /* [x0,x0,x1,x2] */
                            __m512i P3 = c3;
                            __m512i P2 = _mm512_add_epi64(_mm512_set1_epi64(H[k].c2[s1]), vmul(c3, vNA, P));
                            __m512i P1 = _mm512_add_epi64(_mm512_set1_epi64(H[k].E[s1]), vmul(c3, vBc, P));
                            alive[eps] &= vdisc_test(k, P1, P2, P3, NQ3, NQ2, lanes) | (__mmask8)~ok;
                        }
                    }
                    for (int eps = 0; eps < 4; eps++) {
                        __mmask8 m = alive[eps];
                        while (m) {
                            int l = __builtin_ctz(m); m &= m - 1;
                            int idx[4] = {i0, i1, i2 + l, 0}; int ok = 1;
                            for (int kp = NVP; kp < NP + NBIG && ok; kp++) { uint64_t disc; if (hermite_test_mod(kp, idx, pats[0], 3, eps << 1, &disc) == 0) ok = 0; }
                            if (ok) { newpt_t out[2]; int nn = hermite_exact(idx, pats[0], 3, eps << 1, out); for (int q = 0; q < nn; q++) if (nloc < 64) loc[nloc++] = out[q]; }
                            if (nloc >= 62) { FLUSH_LOC(); }
                        }
                    }
                }
            }
            for (int pt = 1; pt < 4; pt++) {
                int k = npat_k[pt]; const int *m = pats[pt];
                int idx[4]; idx[0] = i0;
                if (k == 1) {
                    if (i0 < nstart || is_weier(i0)) continue;
                    uint64_t disc; int ok = 1;
                    for (int kp = 0; kp < NP + NBIG && ok; kp++) if (hermite_test_mod(kp, idx, m, k, 0, &disc) == 0) ok = 0;
                    cnt++;
                    if (ok) { newpt_t out[2]; int nn = hermite_exact(idx, m, k, 0, out); for (int q = 0; q < nn; q++) if (nloc < 64) loc[nloc++] = out[q]; }
                    if (nloc >= 62) { FLUSH_LOC(); }
                } else {   /* k == 2: (2,2) with i1 > i0, (3,1) with i1 != i0 */
                    for (int i1 = 0; i1 < n; i1++) {
                        if (i1 == i0) continue;
                        if (pt == 1 && i1 < i0) continue;
                        if (i0 < nstart && i1 < nstart) continue;
                        if (is_weier(i0) || (pt == 1 && is_weier(i1))) continue;
                        idx[1] = i1;
                        for (int eps = 0; eps < 2; eps++) {
                            uint64_t disc; int ok = 1;
                            for (int kp = 0; kp < NP + NBIG && ok; kp++) if (hermite_test_mod(kp, idx, m, k, eps << 1, &disc) == 0) ok = 0;
                            cnt++;
                            if (ok) { newpt_t out[2]; int nn = hermite_exact(idx, m, k, eps << 1, out); for (int q = 0; q < nn; q++) if (nloc < 64) loc[nloc++] = out[q]; }
                            if (nloc >= 62) { FLUSH_LOC(); }
                        }
                    }
                }
            }
            if (nloc) { FLUSH_LOC(); }
        }
    }
    *ntests += cnt;
}

/* quotient of a degree-5 polynomial r5 x^5 + r4 x^4 + r3 x^3 + ... (accumulators < 2^52) by the monic cubic
   x^3 + Qc2 x^2 + Qc1 x + Qc0 (NQc2 = -Qc2, NQc1 = -Qc1 half-reduced): q2 x^2 + q1 x + q0; returns the mask of lanes
   whose disc(q) is a square mod p or whose q2 is 0 mod p */
static inline __mmask8 vdisc3(int k, __m512i r5, __m512i r4, __m512i r3, __m512i NQc2, __m512i NQc1, __mmask8 lanes) {
    const vprime_t *P = &VP[k];
    __m512i qq2 = vred(r5, P);
    r4 = vmla(r4, qq2, NQc2);
    __m512i qq1 = vred(r4, P);
    r3 = vmla(vmla(r3, qq2, NQc1), qq1, NQc2);
    __m512i qq0 = vred(r3, P);
    __m512i u = vmul(qq2, qq0, P);
    __m512i w = _mm512_sub_epi64(P->v8p, _mm512_slli_epi64(u, 2));
    __m512i disc = vfull(vred(vmla(w, qq1, qq1), P), P);
    __mmask8 q2z = _mm512_cmpeq_epi64_mask(qq2, _mm512_setzero_si512()) | _mm512_cmpeq_epi64_mask(qq2, P->vp);
    __m256i g = _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), lanes, disc, (const int*)QRB[k], 1);
    __mmask8 pass = _mm256_test_epi32_mask(g, _mm256_set1_epi32(0xff));
    return pass | q2z;
}

/* pass over the relations through the points at infinity.  The dominant pattern k = 1, (1,1,1) (contact order 1 at
 * inf_sigma, three finite simple points i0 < i1 < i2, all 8 sign patterns, lanes = i2) is vectorised: p = sigma s x^3 + c(x)
 * with c the quadratic through (x_i, eps_i y_i - sigma s x_i^3) (Newton form from the inverse tables), p^2 - f has degree
 * <= 5 and is divided by the monic cubic prod (x - x_i).  For quintics (f6 = 0) the same code runs with s = 0 (conics).
 * The O(n^2) patterns stay scalar.  Candidates are flushed after every test (a full local buffer would drop points). */
static void inf_pass(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *ntests) {
    if (!f6_is_zero && !has_sqrt_f6) return;
    int nsig = f6_is_zero ? 1 : 2;
    uint64_t cnt = 0;
    #pragma omp parallel reduction(+:cnt)
    {
        newpt_t loc[64]; int nloc = 0;
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 < n; i0++) {
            for (int si = 0; si < nsig; si++) {
                int sigma = si ? -1 : 1;
                int idx[4], m[4];
                idx[0] = i0;
                /* k = 1, (1,1,1): i0 < i1 < i2, i2 >= nstart; vectorised over i2 */
                for (int i1 = i0 + 1; i1 < n; i1++) {
                    int i2start = i1 + 1; if (i2start < nstart) i2start = nstart;
                    if (i2start >= n) continue;
                    struct { uint64_t p3, p3x2, NS01, P01, c1[4], g1[4], E1[4]; int ok; } H[NVP];
                    for (int k = 0; k < NVP; k++) {
                        uint64_t p = PR[k]; double pinv = PINV[k];
                        uint64_t x0 = xr[k][i0], y0 = yr[k][i0], x1 = xr[k][i1], y1 = yr[k][i1], inv10 = invd[k][(size_t)i1 * ncap + i0];
                        uint64_t sk = has_sqrt_f6 ? mpz_fdiv_ui(sqrt_f6z, p) : 0;
                        if (sigma < 0) sk = (p - sk) % p;
                        H[k].ok = inv10 != 0 && (f6_is_zero || sk != 0);
                        H[k].p3 = sk; H[k].p3x2 = (2 * sk) % p;
                        H[k].NS01 = (2 * p - x0 - x1) % p; H[k].P01 = redp(x0 * x1, p, pinv);
                        uint64_t sx0 = sigma < 0 ? (p - sx3[k][i0]) % p : sx3[k][i0], sx1 = sigma < 0 ? (p - sx3[k][i1]) % p : sx3[k][i1];
                        for (int e = 0; e < 4; e++) {   /* bit 0: sign of y0, bit 1: sign of y1 */
                            uint64_t g0 = (((e & 1) ? (p - y0) % p : y0) + p - sx0) % p;
                            uint64_t g1 = (((e & 2) ? (p - y1) % p : y1) + p - sx1) % p;
                            uint64_t c1 = redp(((g1 + p - g0) % p) * inv10, p, pinv);
                            H[k].g1[e] = g1; H[k].c1[e] = c1;
                            H[k].E1[e] = (g0 + p - redp(c1 * x0, p, pinv)) % p;   /* c0 - c1 x0 */
                        }
                    }
                    for (int i2 = i2start; i2 < n; i2 += 8) {
                        int nl = n - i2 < 8 ? n - i2 : 8;
                        __mmask8 lanes = (__mmask8)((1u << nl) - 1);
                        __mmask8 alive[8]; for (int e = 0; e < 8; e++) alive[e] = lanes;
                        cnt += 8 * (uint64_t)nl;
                        for (int k = 0; k < NVP && (alive[0]|alive[1]|alive[2]|alive[3]|alive[4]|alive[5]|alive[6]|alive[7]); k++) {
                            if (!H[k].ok) continue;
                            const vprime_t *P = &VP[k];
                            __m512i X2 = _mm512_maskz_loadu_epi64(lanes, xr[k] + i2), Y2 = _mm512_maskz_loadu_epi64(lanes, yr[k] + i2);
                            __m512i SX2 = _mm512_maskz_loadu_epi64(lanes, sx3[k] + i2);
                            __m512i ninv21 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i1 * ncap + i2);   /* 1/(x1 - x2) */
                            __m512i ninv20 = _mm512_maskz_loadu_epi64(lanes, invd[k] + (size_t)i0 * ncap + i2);   /* 1/(x0 - x2) */
                            __mmask8 ok = _mm512_test_epi64_mask(ninv21, ninv21) & _mm512_test_epi64_mask(ninv20, ninv20);
                            if (!ok) continue;
                            __m512i inv21 = _mm512_sub_epi64(P->vp, ninv21), inv20 = _mm512_sub_epi64(P->vp, ninv20);
                            if (sigma < 0) SX2 = _mm512_sub_epi64(P->vp, SX2);          /* sigma s x2^3, in [0, p] */
                            __m512i NSX2 = _mm512_sub_epi64(P->vp, SX2);                 /* -sigma s x2^3, in [0, p] */
                            __m512i vNS01 = _mm512_set1_epi64(H[k].NS01), vP01 = _mm512_set1_epi64(H[k].P01), v2p3 = _mm512_set1_epi64(H[k].p3x2);
                            /* monic cubic (x-x0)(x-x1)(x-x2): -Qc2 = x0 + x1 + x2, -Qc1 = -(x0 x1 + (x0 + x1) x2) */
                            __m512i NQc2 = _mm512_add_epi64(_mm512_sub_epi64(P->v2p, vNS01), X2);   /* (x0 + x1) + x2, < 4p */
                            __m512i NQc1 = _mm512_sub_epi64(P->v2p, vred(vmla(vP01, _mm512_sub_epi64(P->v2p, vNS01), X2), P));
                            __m512i NF5 = _mm512_set1_epi64(P->p - fr[k][5]), NF4 = _mm512_set1_epi64(P->p - fr[k][4]), NF3 = _mm512_set1_epi64(P->p - fr[k][3]);
                            __m512i g2[2] = { _mm512_add_epi64(Y2, NSX2), _mm512_add_epi64(_mm512_sub_epi64(P->v2p, Y2), NSX2) };   /* +-y2 - sigma s x2^3, < 3p */
                            #pragma GCC unroll 8
                            for (int eps = 0; eps < 8; eps++) {   /* bit i = sign of y_i */
                                if (!alive[eps]) continue;
                                int e01 = eps & 3, s2 = eps >> 2;
                                __m512i d12 = vmul(_mm512_add_epi64(g2[s2], _mm512_set1_epi64(2 * P->p - H[k].g1[e01])), inv21, P);     /* [x1,x2] */
                                __m512i c2 = vmul(_mm512_add_epi64(d12, _mm512_set1_epi64(2 * P->p - H[k].c1[e01])), inv20, P);       /* [x0,x1,x2] */
                                __m512i p2 = c2;
                                __m512i p1 = _mm512_add_epi64(_mm512_set1_epi64(H[k].c1[e01]), vmul(c2, vNS01, P));                   /* c1 - c2 (x0 + x1) */
                                __m512i p0 = _mm512_add_epi64(_mm512_set1_epi64(H[k].E1[e01]), vmul(c2, vP01, P));                    /* c0 - c1 x0 + c2 x0 x1 */
                                __m512i r5 = vmla(NF5, v2p3, p2);
                                __m512i r4 = vmla(vmla(NF4, v2p3, p1), p2, p2);
                                __m512i r3 = vmla(vmla(NF3, v2p3, p0), _mm512_add_epi64(p2, p2), p1);
                                alive[eps] &= vdisc3(k, r5, r4, r3, NQc2, NQc1, lanes) | (__mmask8)~ok;
                            }
                        }
                        for (int eps = 0; eps < 8; eps++) {
                            __mmask8 mm = alive[eps];
                            while (mm) {
                                int l = __builtin_ctz(mm); mm &= mm - 1;
                                int id[4] = {i0, i1, i2 + l, 0}; static const int m111[4] = {1, 1, 1, 0}; int ok = 1;
                                for (int kp = NVP; kp < NP + NBIG && ok; kp++) if (inf_test_mod(kp, id, m111, 3, eps, 1, sigma) == 0) ok = 0;
                                if (ok) { newpt_t out[2]; int nn = inf_exact(id, m111, 3, eps, 1, sigma, out); for (int q = 0; q < nn; q++) if (nloc < 64) loc[nloc++] = out[q]; }
                                if (nloc >= 62) { FLUSH_LOC(); }
                            }
                        }
                    }
                }
                /* k = 1: (2,1) with i0 tangent, i1 any ; (3) i0 */
                for (int i1 = 0; i1 < n; i1++) {
                    if (i1 == i0 || (i0 < nstart && i1 < nstart)) continue;
                    idx[1] = i1; m[0] = 2; m[1] = 1;
                    for (int eps = 0; eps < 4; eps++) { inf_try(idx, m, 2, eps, 1, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
                }
                if (i0 >= nstart) { m[0] = 3; for (int eps = 0; eps < 2; eps++) { inf_try(idx, m, 1, eps, 1, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                if (f6_is_zero) continue;
                /* k = 2: (1,1) i0 < i1 ; (2) i0 */
                for (int i1 = i0 + 1; i1 < n; i1++) {
                    if (i1 < nstart) continue;
                    idx[1] = i1; m[0] = m[1] = 1;
                    for (int eps = 0; eps < 4; eps++) { inf_try(idx, m, 2, eps, 2, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
                }
                if (i0 >= nstart) { m[0] = 2; for (int eps = 0; eps < 2; eps++) { inf_try(idx, m, 1, eps, 2, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                /* k = 3: (1) */
                if (i0 >= nstart) { m[0] = 1; for (int eps = 0; eps < 2; eps++) { inf_try(idx, m, 1, eps, 3, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                /* k = 4: no finite points (once per pass, by thread of i0 == 0 in pass 1) */
                if (i0 == 0 && nstart == 0) { inf_try(idx, m, 0, 0, 4, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
            }
            if (nloc) { FLUSH_LOC(); }
        }
    }
    *ntests += cnt;
}

/* candidate points that are already known are dropped at once */
#define IX_KEEP_NEW(out, nn) do { for (int u_ = 0; u_ < (nn); u_++) { if (kx_find((out)[u_].a, (out)[u_].c) < 0 && nloc < 64) loc[nloc++] = (out)[u_]; else { mpz_clear((out)[u_].a); mpz_clear((out)[u_].b); mpz_clear((out)[u_].c); } } } while (0)
/* the tangent pass with the index sieve: the same relations as hermite_pass, in the same order.  A cubic that meets the
 * curve with multiplicity m_j at the point j gives sum m_j [P_j - inf'] + (two residual points) = 0. */
static void hermite_pass_ix(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *ntests) {
    static const int pats[4][4] = { {2,1,1,0}, {2,2,0,0}, {3,1,0,0}, {4,0,0,0} };
    static const int npat_k[4] = {3, 2, 2, 1};
    int nb = IXnb, nw = ix_nw; __m128i vsh = _mm_cvtsi32_si128(ix_sh);
    uint64_t cnt = 0;
    #pragma omp parallel reduction(+:cnt)
    {
        newpt_t loc[64]; int nloc = 0;
        uint32_t off[64] __attribute__((aligned(64)));
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 < n; i0++) {
            if (!is_weier(i0)) {   /* pattern (2,1,1): i0 tangent, i1 any, i2 > i1; the row gives all i2 with both signs */
                __m512i t1[(IXMAX + 15) / 16], t2[(IXMAX + 15) / 16]; for (int b = 0; b < nb; b++) ix_vadd(b, IXPT1(i0, 0, b), IXPT2(i0, 0, b), IXPT1(i0, 0, b), IXPT2(i0, 0, b), &t1[b], &t2[b]);
                for (int i1 = 0; i1 < n; i1++) {
                    if (i1 == i0) continue;
                    int i2start = i1 + 1; if (i0 < nstart && i1 < nstart && i2start < nstart) i2start = nstart;
                    if (i2start >= n) continue;
                    cnt += 4 * (uint64_t)(n - i2start - (i0 >= i2start ? 1 : 0));
                    for (int s1 = 0; s1 < 2; s1++) {
                        for (int b = 0; b < nb; b++) { __m512i r1, r2; ix_vadd(b, t1[b], t2[b], IXPT1(i1, s1, b), IXPT2(i1, s1, b), &r1, &r2); _mm512_store_si512(off + 16 * b, ix_off(b, r1, r2, vsh)); }
                        for (int w = (2 * i2start) >> 9; w < nw; w++) {
                            uint64_t msk[8], wd[8] __attribute__((aligned(64))); ix_mask(msk, w, i2start, i0);
                            if (!ix_chain(off, w, msk, wd)) continue;
                            for (int q = 0; q < 8; q++) { uint64_t x = wd[q];
                                while (x) {
                                    int bit = 512 * w + 64 * q + __builtin_ctzll(x); x &= x - 1;
                                    int idx[4] = {i0, i1, bit >> 1, 0}, eps = s1 | ((bit & 1) << 1);   /* bit 0: sign of y1, bit 1: sign of y2 */
                                    newpt_t out[2]; int nn = hermite_exact(idx, pats[0], 3, eps << 1, out);
                                    IX_KEEP_NEW(out, nn);
                                    if (nloc >= 62) { FLUSH_LOC(); }
                                }
                            }
                        }
                    }
                }
            }
            for (int pt = 1; pt < 4; pt++) {
                int k = npat_k[pt]; const int *m = pats[pt];
                int idx[4]; idx[0] = i0;
                if (k == 1) {
                    if (i0 < nstart || is_weier(i0)) continue;
                    cnt++;
                    if (ix_test(1, idx, m, 0, 0, 0)) { newpt_t out[2]; int nn = hermite_exact(idx, m, k, 0, out); IX_KEEP_NEW(out, nn); }
                    if (nloc >= 62) { FLUSH_LOC(); }
                } else {   /* k == 2: (2,2) with i1 > i0, (3,1) with i1 != i0 */
                    for (int i1 = 0; i1 < n; i1++) {
                        if (i1 == i0) continue;
                        if (pt == 1 && i1 < i0) continue;
                        if (i0 < nstart && i1 < nstart) continue;
                        if (is_weier(i0) || (pt == 1 && is_weier(i1))) continue;
                        idx[1] = i1;
                        for (int eps = 0; eps < 2; eps++) {
                            cnt++;
                            if (ix_test(2, idx, m, eps << 1, 0, 0)) { newpt_t out[2]; int nn = hermite_exact(idx, m, k, eps << 1, out); IX_KEEP_NEW(out, nn); }
                            if (nloc >= 62) { FLUSH_LOC(); }
                        }
                    }
                }
            }
            if (nloc) { FLUSH_LOC(); }
        }
    }
    *ntests += cnt;
}

/* the relations through the points at infinity (f6 a nonzero square) with the index sieve.  Contact of order k with the
 * point at infinity (1 : s : 0) contributes k times its class.  Only that point is used: the hyperelliptic involution
 * takes a relation through (1 : -s : 0) to one through (1 : s : 0) with all signs changed and the same residual x-values,
 * and all signs are tried.  A relation with three simple points among known points is verified once: the cubic with
 * leading coefficient s through any three of its five finite points is the same. */
static inline void inf_try_ix(const int *idx, const int *m, int nfin, int eps, int k, int sigma, newpt_t *loc, int *nloc_, uint64_t *cnt) {
    for (int i = 0; i < nfin; i++) if (m[i] >= 2 && is_weier(idx[i])) return;
    (*cnt)++;
    if (!ix_test(nfin, idx, m, eps, k, sigma < 0)) return;
    newpt_t out[2]; int nn = inf_exact(idx, m, nfin, eps, k, sigma, out); int nloc = *nloc_;
    IX_KEEP_NEW(out, nn); *nloc_ = nloc;
}
static void inf_pass_ix(int n, int nstart, newpt_t **found, int *nfound, int *cfound, uint64_t *ntests) {
    if (!has_sqrt_f6) return;
    int nb = IXnb, nw = ix_nw; __m128i vsh = _mm_cvtsi32_si128(ix_sh);
    uint64_t cnt = 0;
    #pragma omp parallel reduction(+:cnt)
    {
        newpt_t loc[64]; int nloc = 0;
        uint32_t off[64] __attribute__((aligned(64)));
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 < n; i0++) {
            for (int si = 0; si < 1; si++) {
                int sigma = 1;
                int idx[4], m[4];
                idx[0] = i0;
                /* k = 1, (1,1,1): i0 < i1 < i2, i2 >= nstart; the row gives all i2 with both signs */
                for (int e0 = 0; e0 < 2; e0++) {
                    __m512i a1[(IXMAX + 15) / 16], a2[(IXMAX + 15) / 16]; for (int b = 0; b < nb; b++) ix_vadd(b, IXI1[si][b], IXI2[si][b], IXPT1(i0, e0, b), IXPT2(i0, e0, b), &a1[b], &a2[b]);
                    for (int i1 = i0 + 1; i1 < n; i1++) {
                        int i2start = i1 + 1; if (i2start < nstart) i2start = nstart;
                        if (i2start >= n) continue;
                        if (e0 == 0) cnt += 8 * (uint64_t)(n - i2start);
                        for (int e1 = 0; e1 < 2; e1++) {
                            for (int b = 0; b < nb; b++) { __m512i r1, r2; ix_vadd(b, a1[b], a2[b], IXPT1(i1, e1, b), IXPT2(i1, e1, b), &r1, &r2); _mm512_store_si512(off + 16 * b, ix_off(b, r1, r2, vsh)); }
                            for (int w = (2 * i2start) >> 9; w < nw; w++) {
                                uint64_t msk[8], wd[8] __attribute__((aligned(64))); ix_mask(msk, w, i2start, -1);
                                if (!ix_chain(off, w, msk, wd)) continue;
                                for (int q = 0; q < 8; q++) { uint64_t x = wd[q];
                                    while (x) {
                                        int bit = 512 * w + 64 * q + __builtin_ctzll(x); x &= x - 1;
                                        int id[4] = {i0, i1, bit >> 1, 0}, eps = e0 | (e1 << 1) | ((bit & 1) << 2); static const int m111[4] = {1, 1, 1, 0};   /* bit i = sign of y_i */
                                        if (kq_has(kq_key(0, id[0], id[1], id[2], eps) | (UINT64_C(1) << 62))) continue;
                                        newpt_t out[2]; int nn = inf_exact(id, m111, 3, eps, 1, sigma, out);
                                        int who[2] = {-1, -1}; for (int u = 0; u < nn; u++) who[u] = kx_find(out[u].a, out[u].c);
                                        if (nn == 2 && who[0] >= 0 && who[1] >= 0 && who[0] != who[1] && out[0].ys && out[1].ys) {
                                            int i5[5] = {id[0], id[1], id[2], who[0], who[1]}, s5[5] = {eps & 1, (eps >> 1) & 1, (eps >> 2) & 1, out[0].ys < 0, out[1].ys < 0};
                                            int distinct = 1; for (int a = 0; a < 3; a++) if (i5[a] == who[0] || i5[a] == who[1] || is_weier(i5[a])) distinct = 0;
                                            if (distinct) {
                                                for (int a = 1; a < 5; a++) { int ti = i5[a], ts = s5[a], b = a - 1; while (b >= 0 && i5[b] > ti) { i5[b+1] = i5[b]; s5[b+1] = s5[b]; b--; } i5[b+1] = ti; s5[b+1] = ts; }
                                                for (int a = 0; a < 5; a++) for (int b = a + 1; b < 5; b++) {   /* leave out a and b */
                                                    int t[3], s[3], mm = 0; for (int q2 = 0; q2 < 5; q2++) if (q2 != a && q2 != b) { t[mm] = i5[q2]; s[mm] = s5[q2]; mm++; }
                                                    kq_add(kq_key(0, t[0], t[1], t[2], s[0] | (s[1] << 1) | (s[2] << 2)) | (UINT64_C(1) << 62));
                                                }
                                            }
                                        }
                                        for (int u = 0; u < nn; u++) { if (who[u] < 0 && nloc < 64) loc[nloc++] = out[u]; else { mpz_clear(out[u].a); mpz_clear(out[u].b); mpz_clear(out[u].c); } }
                                        if (nloc >= 62) { FLUSH_LOC(); }
                                    }
                                }
                            }
                        }
                    }
                }
                /* k = 1: (2,1) with i0 tangent, i1 any ; (3) i0 */
                for (int i1 = 0; i1 < n; i1++) {
                    if (i1 == i0 || (i0 < nstart && i1 < nstart)) continue;
                    idx[1] = i1; m[0] = 2; m[1] = 1;
                    for (int eps = 0; eps < 4; eps++) { inf_try_ix(idx, m, 2, eps, 1, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
                }
                if (i0 >= nstart) { m[0] = 3; for (int eps = 0; eps < 2; eps++) { inf_try_ix(idx, m, 1, eps, 1, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                /* k = 2: (1,1) i0 < i1 ; (2) i0 */
                for (int i1 = i0 + 1; i1 < n; i1++) {
                    if (i1 < nstart) continue;
                    idx[1] = i1; m[0] = m[1] = 1;
                    for (int eps = 0; eps < 4; eps++) { inf_try_ix(idx, m, 2, eps, 2, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
                }
                if (i0 >= nstart) { m[0] = 2; for (int eps = 0; eps < 2; eps++) { inf_try_ix(idx, m, 1, eps, 2, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                /* k = 3: (1) */
                if (i0 >= nstart) { m[0] = 1; for (int eps = 0; eps < 2; eps++) { inf_try_ix(idx, m, 1, eps, 3, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } } }
                /* k = 4: no finite points (once per pass, by thread of i0 == 0 in pass 1) */
                if (i0 == 0 && nstart == 0) { inf_try_ix(idx, m, 0, 0, 4, sigma, loc, &nloc, &cnt); if (nloc >= 62) { FLUSH_LOC(); } }
            }
            if (nloc) { FLUSH_LOC(); }
        }
    }
    *ntests += cnt;
}

/* seed points as read (for -K): sorted by height max(|a|, |c|), ties by a then c, so the seed order is deterministic */
typedef struct { mpz_t a, b, c; } inpt_t;
static int inpt_cmp(const void *u, const void *v) {
    const inpt_t *P = u, *Q = v;
    mpz_srcptr hp = mpz_cmpabs(P->a, P->c) >= 0 ? P->a : P->c, hq = mpz_cmpabs(Q->a, Q->c) >= 0 ? Q->a : Q->c;
    int r = mpz_cmpabs(hp, hq); if (r) return r;
    r = mpz_cmp(P->a, Q->a); if (r) return r;
    return mpz_cmp(P->c, Q->c);
}

/* bring (a : b : c) (x = a/c, y = b/c^3) to gcd(a, c) = 1, c > 0, and check that the point is on the curve */
static int point_normalize(mpz_t a, mpz_t b, mpz_t c) {
    mpz_t g, F, cp; mpz_init(g); mpz_init(F); mpz_init(cp);
    mpz_gcd(g, a, c);
    if (mpz_cmp_ui(g, 1) != 0) {   /* b/c^3 = (b/g^3)/(c/g)^3: b must be divisible by g^3 */
        mpz_pow_ui(F, g, 3);
        if (!mpz_divisible_p(b, F)) { mpz_clear(g); mpz_clear(F); mpz_clear(cp); return 0; }
        mpz_divexact(a, a, g); mpz_divexact(c, c, g); mpz_divexact(b, b, F);
    }
    if (mpz_sgn(c) < 0) { mpz_neg(a, a); mpz_neg(b, b); mpz_neg(c, c); }
    mpz_set(F, f[6]); mpz_set_ui(cp, 1);
    for (int d = 5; d >= 0; d--) { mpz_mul(cp, cp, c); mpz_mul(F, F, a); mpz_addmul(F, f[d], cp); }
    mpz_mul(g, b, b);
    int ok = mpz_cmp(g, F) == 0;
    mpz_clear(g); mpz_clear(F); mpz_clear(cp);
    return ok;
}

/* Codex round 48: integer options/arguments are validated as long long against their range BEFORE any narrowing cast */
static long long chk_ll(const char *what, const char *s, long long lo, long long hi) { char *end = NULL; errno = 0; long long v = strtoll(s, &end, 10); if (errno != 0 || end == s || *end != 0) { fprintf(stderr, "%s: invalid integer '%s'\n", what, s); exit(1); } if (v < lo || v > hi) { fprintf(stderr, "%s: %lld is outside [%lld, %lld]\n", what, v, lo, hi); exit(1); } return v; }
static int arg_int(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return (int)chk_ll(opt, argv[++*i], lo, hi); }
static __attribute__((unused)) long long arg_llr(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return chk_ll(opt, argv[++*i], lo, hi); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s 'f0 .. f6' points.txt [-t threads] [-maxpass k] [-notangent] [-noinf] [-six n] [-K k]\n", argv[0]); return 1; }
    for (int d = 0; d < 7; d++) mpz_init(f[d]);
    { char *s = strdup(argv[1]); int i = 0; char *tok = strtok(s, " ,"); while (tok) { if (i == 7 || mpz_set_str(f[i], tok, 10) != 0) { fprintf(stderr, "model: give exactly seven integer coefficients f0 .. f6, not '%s'\n", argv[1]); return 1; } i++; tok = strtok(NULL, " ,"); } if (i != 7) { fprintf(stderr, "model: give exactly seven integer coefficients f0 .. f6, not '%s'\n", argv[1]); return 1; } free(s); }
    int nthreads = omp_get_max_threads(), maxpass = 10, use_tangent = 1, use_inf = 1; uint64_t six = 0; int nofilter = getenv("CUBIC_NOFILTER") != NULL;
    check_exact = getenv("CUBIC_CHECKEXACT") != NULL;
    for (int i = 3; i < argc; i++) { if (!strcmp(argv[i], "-t")) nthreads = arg_int(argc, argv, &i, 1, 4096); else if (!strcmp(argv[i], "-maxpass")) maxpass = arg_int(argc, argv, &i, 0, 1000); else if (!strcmp(argv[i], "-notangent")) use_tangent = 0; else if (!strcmp(argv[i], "-noinf")) use_inf = 0; else if (!strcmp(argv[i], "-six")) six = (uint64_t)arg_llr(argc, argv, &i, 0, 1LL << 62); else if (!strcmp(argv[i], "-K")) KRES = arg_int(argc, argv, &i, 0, 1 << 20); else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; } }
    omp_set_num_threads(nthreads);
    for (int k = 0; k < NP+NBIG; k++) { uint64_t p = k < NP ? PR[k] : BIGP[k-NP]; PINV[k] = 1.0 / (double)p; for (int d = 0; d < 7; d++) fr[k][d] = mpz_fdiv_ui(f[d], p); }
    for (int k = 0; k < NVP; k++) { VP[k].p = PR[k]; VP[k].M = (UINT64_C(1) << 52) / PR[k]; VP[k].vp = _mm512_set1_epi64(PR[k]); VP[k].vM = _mm512_set1_epi64(VP[k].M); VP[k].v2p = _mm512_set1_epi64(2 * PR[k]); VP[k].v8p = _mm512_set1_epi64(8 * PR[k]); }
    f6_is_zero = (mpz_sgn(f[6]) == 0);
    mpz_init(sqrt_f6z);
    if (!f6_is_zero && mpz_sgn(f[6]) > 0 && mpz_perfect_square_p(f[6])) { mpz_sqrt(sqrt_f6z, f[6]); has_sqrt_f6 = 1; }
    fprintf(stderr, "infinity relations: %s\n", f6_is_zero ? "quintic (conics)" : has_sqrt_f6 ? "two rational points at infinity" : "none");
    for (int k = 0; k < NP; k++) { QR[k] = calloc(PR[k]/8 + 1, 1); QRB[k] = calloc(PR[k] + 8, 1); for (uint64_t x = 0; x < PR[k]; x++) { uint64_t s = mulmod(x, x, PR[k]); QR[k][s >> 3] |= 1 << (s & 7); QRB[k][s] = 1; } }
    /* read points (with -K: collected first and sorted by height max(|a|, |c|), so that indices < K are the K lowest) */
    FILE *fp = fopen(argv[2], "r"); if (!fp) { perror("points"); return 1; }
    char line[4096]; mpz_t a, b, c; mpz_init(a); mpz_init(b); mpz_init(c); int nbad = 0;
    inpt_t *inp = NULL; int ninp = 0, cinp = 0;
    while (fgets(line, sizeof line, fp)) {
        char *s = line; while (*s == '(' || *s == ' ') s++;
        char *sa = strtok(s, " :)\n"), *sb = strtok(NULL, " :)\n"), *sc = strtok(NULL, " :)\n");
        if (!sa || !sb || !sc) continue;
        mpz_set_str(a, sa, 10); mpz_set_str(b, sb, 10); mpz_set_str(c, sc, 10);
        if (mpz_sgn(c) == 0) continue;                   /* skip points at infinity */
        if (mpz_sgn(b) < 0) continue;                    /* one representative per x */
        if (!point_normalize(a, b, c)) { nbad++; continue; }
        if (KRES <= 0) { if (!known_x(a, c)) add_pt(a, b, c); continue; }
        if (ninp == cinp) { cinp = cinp ? 2 * cinp : 1024; inp = realloc(inp, cinp * sizeof(inpt_t)); }
        mpz_init_set(inp[ninp].a, a); mpz_init_set(inp[ninp].b, b); mpz_init_set(inp[ninp].c, c); ninp++;
    }
    fclose(fp);
    if (KRES > 0) {
        qsort(inp, ninp, sizeof(inpt_t), inpt_cmp);
        for (int i = 0; i < ninp; i++) { if (!known_x(inp[i].a, inp[i].c)) add_pt(inp[i].a, inp[i].b, inp[i].c); mpz_clear(inp[i].a); mpz_clear(inp[i].b); mpz_clear(inp[i].c); }
        free(inp);
    }
    if (nbad) fprintf(stderr, "warning: %d input points are not on the curve (skipped)\n", nbad);
    fprintf(stderr, "read %d x-values\n", npts);
    if (KRES > 0) {
        /* indices >= npts are points found by the closure: never let them into the two restricted slots */
        if (KRES > npts) { fprintf(stderr, "seed has only %d x-values: K reduced from %d to %d\n", npts, KRES, npts); KRES = npts; }
        fprintf(stderr, "restricted closure: four-point pass limited to subsets with two points among the %d lowest-height seed points (tangent and infinity passes unrestricted)\n", KRES);
    }
    /* the index sieve is used unless the closure is restricted (-K), it is switched off (CUBIC_OLD), or too few primes are usable */
    int use_ix = 0;
    if (KRES <= 0 && !f6_is_zero && !getenv("CUBIC_OLD") && !getenv("CUBIC_SCALAR")) {
        double tix = omp_get_wtime();
        ix_init(getenv("CUBIC_IXPMAX") ? atoi(getenv("CUBIC_IXPMAX")) : 257, getenv("CUBIC_IXTARGET") ? atof(getenv("CUBIC_IXTARGET")) : 1e-5, 1);
        double prod = 1; for (int k = 0; k < nix; k++) prod *= IX[k].dens;
        use_ix = nix >= 4 && prod < 1e-3; ix_pretest = prod > 3e-5 ? 4 : 0;
        fprintf(stderr, "index sieve %s (set up in %.3fs)\n", use_ix ? "used" : "not used", omp_get_wtime() - tix);
    }
    int nstart = 0;    /* points with index >= nstart are "new" in the current pass */
    for (int pass = 1; pass <= maxpass; pass++) {
        int n = npts;
        double t0 = omp_get_wtime();
        uint64_t nsub = 0, nsurv1 = 0, nexact = 0; int nfound = 0;
        newpt_t *found = NULL; int cfound = 0;
        if (use_ix) { ix_rows(n); kx_build(n); four_pass_ix(n, nstart, &found, &nfound, &cfound, &nsub, &nsurv1, &nexact); }
        else if (!getenv("CUBIC_SCALAR")) four_pass_vec(n, nstart, &found, &nfound, &cfound, &nsub, &nsurv1, &nexact); else {
        #pragma omp parallel reduction(+:nsub, nsurv1, nexact)
        {
            newpt_t loc[64]; int nloc = 0;
            #pragma omp for schedule(dynamic, 1)
            for (int i0 = 0; i0 < n; i0++) {
                for (int i1 = i0 + 1; i1 < n && (KRES <= 0 || i1 < KRES); i1++) for (int i2 = i1 + 1; i2 < n; i2++) for (int i3 = i2 + 1; i3 < n; i3++) {
                    if (i3 < nstart) continue;           /* all four old: done in an earlier pass */
                    int idx[4] = {i0, i1, i2, i3};
                    sub_t S0; int p0ok = prep_subset(0, idx, &S0);
                    nsub++;
                    for (int eps = 0; eps < 8; eps++) {   /* sign of y_0 fixed */
                        uint64_t disc;
                        if (!nofilter && p0ok && !test_pattern(0, &S0, eps, &disc)) continue;
                        nsurv1++;
                        int ok = 1;
                        if (!nofilter) for (int k = 1; k < NP + NBIG && ok; k++) { sub_t S; if (!prep_subset(k, idx, &S)) continue; ok = test_pattern(k, &S, eps, &disc); }
                        if (!ok) continue;
                        nexact++;
                        newpt_t out[2]; int nn = exact_check(idx, eps, out);
                        for (int q = 0; q < nn; q++) { if (nloc < 64) loc[nloc++] = out[q]; }
                        if (nloc >= 62) {
                            #pragma omp critical
                            { for (int q = 0; q < nloc; q++) { if (nfound == cfound) { cfound = cfound ? 2*cfound : 64; found = realloc(found, cfound * sizeof(newpt_t)); } found[nfound++] = loc[q]; } nloc = 0; }
                        }
                    }
                }
                if (nloc) {
                    #pragma omp critical
                    { for (int q = 0; q < nloc; q++) { if (nfound == cfound) { cfound = cfound ? 2*cfound : 64; found = realloc(found, cfound * sizeof(newpt_t)); } found[nfound++] = loc[q]; } nloc = 0; }
                }
            }
        }
        }
        uint64_t nherm = 0; double t1 = omp_get_wtime();
        if (use_tangent) { if (use_ix && !getenv("CUBIC_OLDTAN")) hermite_pass_ix(n, nstart, &found, &nfound, &cfound, &nherm); else hermite_pass(n, nstart, &found, &nfound, &cfound, &nherm); }
        double t2 = omp_get_wtime();
        if (use_inf) { if (use_ix && !getenv("CUBIC_OLDTAN")) inf_pass_ix(n, nstart, &found, &nfound, &cfound, &nherm); else inf_pass(n, nstart, &found, &nfound, &cfound, &nherm); }
        double t3 = omp_get_wtime();
        /* add genuinely new points */
        int added = 0;
        for (int q = 0; q < nfound; q++) {
            if (!known_x(found[q].a, found[q].c)) {
                add_pt(found[q].a, found[q].b, found[q].c); added++;
                gmp_printf("(%Zd : %Zd : %Zd)  # pass %d from x-indices %d %d %d %d eps %d\n", found[q].a, found[q].b, found[q].c, pass, found[q].idx[0], found[q].idx[1], found[q].idx[2], found[q].idx[3], found[q].eps);
                fflush(stdout);
            }
        }
        fprintf(stderr, "pass %d: %lu subsets, %lu survived %s, %lu exact checks, %lu tangent tests, %d candidate points, %d new x-values; total %d; %.1fs\n", pass, nsub, nsurv1, use_ix ? "the index sieve" : getenv("CUBIC_SCALAR") ? "prime 1" : "the vector primes", nexact, nherm, nfound, added, npts, omp_get_wtime() - t0);
        if (getenv("CUBIC_IXSTAT")) fprintf(stderr, "  seconds: four-point %.3f, tangent %.3f, infinity %.3f, new points %.3f\n", t1 - t0, t2 - t1, t3 - t2, omp_get_wtime() - t3);
        free(found);
        if (added == 0) break;
        nstart = n;
    }
    if (six) {
        newpt_t *found = NULL; int nfound = 0, cfound = 0; uint64_t nex = 0;
        double t0 = omp_get_wtime();
        six_pass(npts, six, &found, &nfound, &cfound, &nex);
        int added = 0;
        for (int q = 0; q < nfound; q++) if (!known_x(found[q].a, found[q].c)) { add_pt(found[q].a, found[q].b, found[q].c); added++; gmp_printf("(%Zd : %Zd : %Zd)  # six-point relation\n", found[q].a, found[q].b, found[q].c); }
        fprintf(stderr, "six-point pass: %d candidate points, %d new x-values; total %d; %.1fs\n", nfound, added, npts, omp_get_wtime() - t0);
        free(found);
    }
    { int nw = 0; for (int i = 0; i < npts; i++) nw += is_weier(i); int ninfty = f6_is_zero ? 1 : (has_sqrt_f6 ? 2 : 0);
      fprintf(stderr, "final: %d x-values (%d rational points: %d finite, %d at infinity)\n", npts, 2*npts - nw + ninfty, 2*npts - nw, ninfty); }
    return 0;
}

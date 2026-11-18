/* conicmodp.c -- all smooth plane conics over F_p on which a given plane sextic S is a constant times a square.
 *
 *   conicmodp sextic.txt p [-t threads] > conics_mod_p.txt
 *
 * Sextic file: 28 integer coefficients in the order i = 6..0, j = 6-i..0, k = 6-i-j for x0^i x1^j x2^k.
 * Every conic q = c00 x0^2 + c01 x0 x1 + c02 x0 x2 + c11 x1^2 + c12 x1 x2 + c22 x2^2 over F_p, normalised so that its
 * first nonzero coefficient (in that order) is 1, with det(Gram) != 0, is parametrised through one of its points and
 * S(Q(s,t)) is evaluated at the p+1 points of P^1(F_p): the conic survives if all nonzero values have the same
 * quadratic character (early exit at the first mismatch, ~2.5 evaluations on average); survivors are then tested
 * exactly (S(Q(s,t)) = c G(s,t)^2 as a binary form over F_p, all root multiplicities even).
 * Output lines: "c00 c01 c02 c11 c12 c22 chi" with chi = +1 if the constant c is a square mod p, -1 otherwise.
 * (p must be an odd prime, p <= 1000.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

static int64_t p;
static int64_t Sc[28]; static int SI[28], SJ[28], SK[28];
static int8_t *chi;        /* quadratic character table */
static int64_t *sqrt_tab;  /* sqrt_tab[a] = a square root of a, or -1 */
static int64_t *inv_tab;

static int64_t md(int64_t a) { a %= p; return a < 0 ? a + p : a; }
static int64_t pw(int64_t a, int64_t e) { int64_t r = 1; a = md(a); while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; } return r; }

/* S at a point (mod p) */
static inline int64_t evalS(int64_t x0, int64_t x1, int64_t x2) {
    int64_t p0[7], p1[7], p2[7];
    p0[0] = p1[0] = p2[0] = 1;
    for (int e = 1; e <= 6; e++) { p0[e] = p0[e-1] * x0 % p; p1[e] = p1[e-1] * x1 % p; p2[e] = p2[e-1] * x2 % p; }
    int64_t s = 0;
    for (int n = 0; n < 28; n++) if (Sc[n]) s = (s + Sc[n] * (p0[SI[n]] * p1[SJ[n]] % p) % p * p2[SK[n]]) % p;
    return s;
}

/* ---- small polynomial arithmetic mod p (dense, coefficient arrays, degree <= 12) ---- */
#define MAXD 30
typedef struct { int d; int64_t c[MAXD + 1]; } poly;   /* d = -1 for zero */

static void pnorm(poly *a) { while (a->d >= 0 && a->c[a->d] == 0) a->d--; }
static poly pmul(poly *a, poly *b) {
    poly r; r.d = -1; if (a->d < 0 || b->d < 0) return r;
    r.d = a->d + b->d; for (int i = 0; i <= r.d; i++) r.c[i] = 0;
    for (int i = 0; i <= a->d; i++) for (int j = 0; j <= b->d; j++) r.c[i+j] = (r.c[i+j] + a->c[i] * b->c[j]) % p;
    pnorm(&r); return r;
}
static poly pderiv(poly *a) { poly r; r.d = a->d - 1; for (int i = 1; i <= a->d; i++) r.c[i-1] = a->c[i] * i % p; pnorm(&r); if (a->d <= 0) r.d = -1; return r; }
/* a mod b, b != 0 */
static poly pmod(poly a, poly *b) {
    int64_t ib = inv_tab[b->c[b->d]];
    while (a.d >= b->d) {
        int64_t f = a.c[a.d] * ib % p; int sh = a.d - b->d;
        for (int i = 0; i <= b->d; i++) a.c[i+sh] = md(a.c[i+sh] - f * b->c[i]);
        pnorm(&a);
    }
    return a;
}
static poly pdiv(poly a, poly *b) {   /* exact quotient assumed (remainder discarded) */
    poly q; q.d = a.d - b->d; if (q.d < 0) { q.d = -1; return q; }
    for (int i = 0; i <= q.d; i++) q.c[i] = 0;
    int64_t ib = inv_tab[b->c[b->d]];
    while (a.d >= b->d) {
        int64_t f = a.c[a.d] * ib % p; int sh = a.d - b->d; q.c[sh] = f;
        for (int i = 0; i <= b->d; i++) a.c[i+sh] = md(a.c[i+sh] - f * b->c[i]);
        pnorm(&a);
    }
    pnorm(&q); return q;
}
static poly pgcd(poly a, poly b) {
    while (b.d >= 0) { poly r = pmod(a, &b); a = b; b = r; }
    return a;
}
/* is the univariate polynomial a (= F(1,t)) a constant times a square, with (12 - deg) even (the multiplicity at t = infinity)? */
static int is_csquare(poly F) {
    if (F.d < 0) return 0;              /* identically zero: exclude */
    if ((12 - F.d) & 1) return 0;
    while (F.d > 0) {
        poly D = pderiv(&F);
        if (D.d < 0) return 0;          /* cannot happen for p > 12 */
        poly G = pgcd(F, D);            /* prod P_i^(m_i-1) */
        poly R = pdiv(F, &G);           /* radical */
        if (2 * R.d > F.d) return 0;
        poly R2 = pmul(&R, &R);
        poly rem = pmod(F, &R2);
        if (rem.d >= 0) return 0;
        F = pdiv(F, &R2);
    }
    return 1;
}

/* parametrisation data of a conic through point P: Q(s,t) = q(M) P - 2 B(P,M) M, M = s M1 + t M2 */
typedef struct { int64_t A0[3], A1[3], A2[3]; } param;   /* Q_i(s,t) = A0[i] s^2 + A1[i] s t + A2[i] t^2 */

static inline int64_t qf(const int64_t *c, const int64_t *x) {
    return md(c[0]*x[0]%p*x[0] + c[1]*x[0]%p*x[1] + c[2]*x[0]%p*x[2] + c[3]*x[1]%p*x[1] + c[4]*x[1]%p*x[2] + c[5]*x[2]%p*x[2]);
}
/* polar form: 2B(x,y) = q(x+y) - q(x) - q(y) */
static inline int64_t twoB(const int64_t *c, const int64_t *x, const int64_t *y) {
    int64_t z[3] = {md(x[0]+y[0]), md(x[1]+y[1]), md(x[2]+y[2])};
    return md(qf(c, z) - qf(c, x) - qf(c, y));
}

/* find a point on the conic; returns 0 if none found (should not happen for nondegenerate conics) */
static int find_point(const int64_t *c, int64_t *P) {
    if (c[0] == 0) { P[0] = 1; P[1] = 0; P[2] = 0; return 1; }
    /* c00 x0^2 + (c01 x1 + c02 x2) x0 + (c11 x1^2 + c12 x1 x2 + c22 x2^2) = 0 */
    int64_t i2a = inv_tab[md(2 * c[0])];
    for (int64_t k = 0; k <= p; k++) {
        int64_t x1, x2;
        if (k < p) { x1 = 1; x2 = k; } else { x1 = 0; x2 = 1; }
        int64_t b = md(c[1]*x1 + c[2]*x2), cc = md(c[3]*x1%p*x1 + c[4]*x1%p*x2 + c[5]*x2%p*x2);
        int64_t disc = md(b*b - 4*c[0]%p*cc);
        int64_t r = sqrt_tab[disc];
        if (r >= 0) { P[0] = md((-b + r) * i2a); P[1] = x1; P[2] = x2; return 1; }
    }
    return 0;
}

static void make_param(const int64_t *c, const int64_t *P, param *pr) {
    int64_t M1[3] = {0,0,0}, M2[3] = {0,0,0};
    if (P[0]) { M1[1] = 1; M2[2] = 1; } else if (P[1]) { M1[0] = 1; M2[2] = 1; } else { M1[0] = 1; M2[1] = 1; }
    /* q(s M1 + t M2) = q(M1) s^2 + 2B(M1,M2) s t + q(M2) t^2 ; 2B(P, s M1 + t M2) = 2B(P,M1) s + 2B(P,M2) t */
    int64_t qa = qf(c, M1), qb = twoB(c, M1, M2), qc = qf(c, M2);
    int64_t b1 = twoB(c, P, M1), b2 = twoB(c, P, M2);
    for (int i = 0; i < 3; i++) {
        /* Q_i = (qa s^2 + qb s t + qc t^2) P_i - (b1 s + b2 t)(M1_i s + M2_i t) */
        pr->A0[i] = md(qa * P[i] - b1 * M1[i]);
        pr->A1[i] = md(qb * P[i] - b1 * M2[i] - b2 * M1[i]);
        pr->A2[i] = md(qc * P[i] - b2 * M2[i]);
    }
}

static inline void param_pt(const param *pr, int64_t s, int64_t t, int64_t *x) {
    int64_t s2 = s*s%p, st = s*t%p, t2 = t*t%p;
    for (int i = 0; i < 3; i++) x[i] = md(pr->A0[i]*s2 + pr->A1[i]*st + pr->A2[i]*t2);
}

/* exact polynomial F(1,t) = S(Q(1,t)) via polynomial arithmetic */
static poly Fpoly(const param *pr) {
    poly Q[3]; poly P0[7], P1[7], P2[7];
    for (int i = 0; i < 3; i++) { Q[i].d = 2; Q[i].c[0] = pr->A0[i]; Q[i].c[1] = pr->A1[i]; Q[i].c[2] = pr->A2[i]; pnorm(&Q[i]); }
    poly one; one.d = 0; one.c[0] = 1;
    P0[0] = P1[0] = P2[0] = one;
    for (int e = 1; e <= 6; e++) { P0[e] = pmul(&P0[e-1], &Q[0]); P1[e] = pmul(&P1[e-1], &Q[1]); P2[e] = pmul(&P2[e-1], &Q[2]); }
    poly F; F.d = 12; for (int i = 0; i <= 12; i++) F.c[i] = 0;
    for (int n = 0; n < 28; n++) {
        if (!Sc[n]) continue;
        poly m = pmul(&P0[SI[n]], &P1[SJ[n]]); m = pmul(&m, &P2[SK[n]]);
        int64_t co = md(Sc[n]);
        for (int i = 0; i <= m.d; i++) F.c[i] = (F.c[i] + co * m.c[i]) % p;
    }
    pnorm(&F); return F;
}

static int test_conic(const int64_t *c, int *chiout) {
    int64_t P[3]; param pr;
    if (!find_point(c, P)) return 0;
    make_param(c, P, &pr);
    int sign = 0;
    for (int64_t k = 0; k <= p; k++) {
        int64_t s, t; if (k < p) { s = 1; t = k; } else { s = 0; t = 1; }
        int64_t x[3]; param_pt(&pr, s, t, x);
        int64_t v = evalS(x[0], x[1], x[2]);
        if (v == 0) continue;
        int ch = chi[v];
        if (sign == 0) sign = ch; else if (ch != sign) return 0;
    }
    if (sign == 0) return 0;   /* S vanishes identically on the conic: conic component of the branch curve; ignore */
    /* exact test */
    poly F = Fpoly(&pr);
    if (!is_csquare(F)) return 0;
    *chiout = sign;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: conicmodp sextic.txt p [-t threads]\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("sextic"); return 1; }
    { int n = 0; for (int i = 6; i >= 0; i--) for (int j = 6 - i; j >= 0; j--) { SI[n] = i; SJ[n] = j; SK[n] = 6 - i - j; n++; } }
    for (int n = 0; n < 28; n++) if (fscanf(f, "%ld", &Sc[n]) != 1) { fprintf(stderr, "bad sextic file\n"); return 1; }
    fclose(f);
    p = atol(argv[2]);
    int nthreads = 8;
    for (int i = 3; i < argc; i++) if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
    omp_set_num_threads(nthreads);
    if (p < 13 || p > 1000 || (p & 1) == 0) { fprintf(stderr, "p must be an odd prime in [13,1000]\n"); return 1; }
    for (int n = 0; n < 28; n++) Sc[n] = md(Sc[n]);
    chi = malloc(p); sqrt_tab = malloc(p * sizeof(int64_t)); inv_tab = malloc(p * sizeof(int64_t));
    for (int64_t a = 0; a < p; a++) { chi[a] = 0; sqrt_tab[a] = -1; inv_tab[a] = 0; }
    for (int64_t a = 1; a < p; a++) { int64_t s = a * a % p; chi[s] = 1; sqrt_tab[s] = a; inv_tab[a] = pw(a, p - 2); }
    for (int64_t a = 1; a < p; a++) if (chi[a] == 0) chi[a] = -1;
    sqrt_tab[0] = 0; chi[0] = 0;

    int64_t total = 0, found = 0;
    /* normalised conics: first nonzero coefficient 1 at position i0; enumerate remaining 5 - i0 coefficients */
    for (int i0 = 0; i0 < 6; i0++) {
        int nfree = 5 - i0;
        int64_t count = 1; for (int i = 0; i < nfree; i++) count *= p;
        int64_t nf = 0;
        #pragma omp parallel for schedule(dynamic, 1024) reduction(+:nf)
        for (int64_t idx = 0; idx < count; idx++) {
            int64_t c[6]; int64_t r = idx;
            for (int i = 0; i < i0; i++) c[i] = 0;
            c[i0] = 1;
            for (int i = i0 + 1; i < 6; i++) { c[i] = r % p; r /= p; }
            /* determinant of the Gram matrix [[2c00,c01,c02],[c01,2c11,c12],[c02,c12,2c22]] */
            int64_t a = md(2*c[0]), b = c[1], cc = c[2], d = md(2*c[3]), e = c[4], g = md(2*c[5]);
            int64_t det = md(a*(d*g%p - e*e%p) - b*(b*g%p - e*cc%p) + cc*(b*e%p - d*cc%p));
            if (det == 0) continue;
            int chiv;
            if (test_conic(c, &chiv)) {
                nf++;
                #pragma omp critical
                printf("%ld %ld %ld %ld %ld %ld %d\n", c[0], c[1], c[2], c[3], c[4], c[5], chiv);
            }
        }
        total += count; found += nf;
    }
    fprintf(stderr, "p=%ld: %ld normalised conics scanned, %ld smooth conics with S = c*square\n", p, total, found);
    return 0;
}

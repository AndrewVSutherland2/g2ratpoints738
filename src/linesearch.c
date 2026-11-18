/*
 * linesearch - F_p prefilter over lines A x0 + B x1 + C x2 = 0 in P^2 for a K3 double plane y^2 = S(x0,x1,x2).
 * The restriction of S to a line (parametrised by a primitive integer basis of the line) is a binary sextic,
 * i.e. a genus 2 curve; its score is sum_{p<=B, good} log(#C(F_p)/p) computed entirely mod p.
 *
 * usage: linesearch sextic.txt L [-B 200] [-K 1000] [-t threads] [-through a b c] [-exact A B C]
 *   sextic.txt: 28 integers = coefficients of x0^i x1^j x2^k in the order (i,j,k) lexicographic with i+j+k=6:
 *               for i=6..0, for j=6-i..0, k=6-i-j   (i.e. x0^6, x0^5 x1, x0^5 x2, x0^4 x1^2, x0^4 x1 x2, x0^4 x2^2, ...)
 *   L: enumerate primitive (A,B,C) with max(|A|,|B|,|C|) <= L, A >= 0 (sign normalised)
 *   -through a b c: only lines through the point (a:b:c)
 *   -exact A B C: print the integer binary sextic for that line (GMP) and exit
 * output: top K lines: score  nagao  A B C
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>

typedef unsigned __int128 u128;
static int64_t S[28]; static int SI[28], SJ[28], SK[28];   /* exponents */
static int primes[64], nprimes = 0; static int8_t *leg[64]; static uint64_t Smod[64][28]; static uint64_t pmagic[64];
static inline uint64_t fmod_p(uint64_t x, int ip) { uint64_t low = pmagic[ip] * x; return (uint64_t)(((u128)low * (uint64_t)primes[ip]) >> 64); }

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) { return (a * b) % p; }
static int64_t egcd(int64_t a, int64_t b, int64_t *x, int64_t *y) { if (b == 0) { *x = 1; *y = 0; return a; } int64_t x1, y1; int64_t g = egcd(b, a % b, &x1, &y1); *x = y1; *y = x1 - (a / b) * y1; return g; }
static int64_t gcd3(int64_t a, int64_t b, int64_t c) { int64_t x, y; return llabs(egcd(llabs(egcd(llabs(a), llabs(b), &x, &y)), llabs(c), &x, &y)); }

/* primitive integer basis P, Q of the line A x0 + B x1 + C x2 = 0 (gcd(A,B,C) = 1) */
static int use_sym = 1;   /* -nosym disables the x0 -> -x0 symmetry filter in -enum */

static void line_basis(int64_t A, int64_t B, int64_t C, int64_t P[3], int64_t Q[3]) {
    int64_t x, y;
    if (B == 0 && C == 0) { P[0] = 0; P[1] = 1; P[2] = 0; Q[0] = 0; Q[1] = 0; Q[2] = 1; return; }
    int64_t g = llabs(egcd(llabs(B), llabs(C), &x, &y));       /* g = gcd(|B|, |C|) */
    Q[0] = 0; Q[1] = C / g; Q[2] = -B / g;
    /* u, v with (B/g) u + (C/g) v = 1: solve on absolute values, then fix signs */
    int64_t u, v; int64_t g2 = egcd(llabs(B / g), llabs(C / g), &u, &v);
    if (g2 < 0) { u = -u; v = -v; }
    if (B < 0) u = -u;
    if (C < 0) v = -v;
    P[0] = g; P[1] = -A * u; P[2] = -A * v;
    /* sanity: both points on the line */
    if (A * P[0] + B * P[1] + C * P[2] != 0 || A * Q[0] + B * Q[1] + C * Q[2] != 0) { fprintf(stderr, "line_basis error for (%ld,%ld,%ld)\n", A, B, C); exit(1); }
}

/* binary sextic (7 coefficients of t^6..s^6 in ascending powers of s? we use coefficient of t^(6-m) s^m, m=0..6) mod p */
static void restrict_mod(int ip, const int64_t P[3], const int64_t Q[3], uint64_t out[7]) {
    uint64_t p = primes[ip]; (void)p;
    uint64_t Pm[3], Qm[3];
    for (int i = 0; i < 3; i++) { Pm[i] = ((P[i] % (int64_t)p) + p) % p; Qm[i] = ((Q[i] % (int64_t)p) + p) % p; }
    /* powers of the linear forms (t Pm[i] + s Qm[i])^e as polynomials in s (coefficient of t^(e-m) s^m) */
    uint64_t pw[3][7][7];   /* pw[i][e][m] */
    for (int i = 0; i < 3; i++) {
        pw[i][0][0] = 1;
        for (int e = 1; e <= 6; e++) {
            for (int m = 0; m <= e; m++) {
                uint64_t v = 0;
                if (m < e) v += pw[i][e-1][m] * Pm[i];
                if (m > 0) v += pw[i][e-1][m-1] * Qm[i];
                pw[i][e][m] = fmod_p(v, ip);
            }
        }
    }
    for (int m = 0; m < 7; m++) out[m] = 0;
    for (int n = 0; n < 28; n++) {
        uint64_t c = Smod[ip][n]; if (!c) continue;
        int i = SI[n], j = SJ[n], k = SK[n];
        /* product of pw[0][i], pw[1][j], pw[2][k] */
        uint64_t t1[7] = {0,0,0,0,0,0,0};
        for (int a = 0; a <= i; a++) for (int b = 0; b <= j; b++) t1[a+b] += pw[0][i][a] * pw[1][j][b];
        for (int a = 0; a <= i + j; a++) t1[a] = fmod_p(t1[a], ip);
        for (int a = 0; a <= i + j; a++) if (t1[a]) { uint64_t ct = fmod_p(c * t1[a], ip); for (int b = 0; b <= k; b++) out[a+b] = fmod_p(out[a+b] + ct * pw[2][k][b], ip); }
    }
}

/* a^(p-2) mod p (p prime, a != 0 mod p) */
static uint64_t inv_mod_p(uint64_t a, int ip) {
    uint64_t p = primes[ip], r = 1, e = p - 2; a = fmod_p(a, ip);
    while (e) { if (e & 1) r = fmod_p(r * a, ip); a = fmod_p(a * a, ip); e >>= 1; }
    return r;
}
/* degree of gcd(g, g') over F_p for g of degree <= 6 (coefficients g[0..6] of x^0..x^6, reduced mod p); -1 if g = 0 mod p.
 * Euclid on dense polynomials; the degree is all the bad-reduction test needs (Codex round 41: the earlier test only looked
 * for common roots in F_p and missed repeated irreducible factors of degree >= 2, e.g. (x^2+1)^2 (x^2+2) mod 3). */
static int gcd_deg_mod_p(const uint64_t g[7], int ip) {
    uint64_t a[8], b[8]; int da = -1, db = -1;
    for (int j = 0; j <= 6; j++) { a[j] = fmod_p(g[j], ip); if (a[j]) da = j; }
    if (da < 0) return -1;
    for (int j = 0; j <= 5; j++) { b[j] = fmod_p((uint64_t)(j + 1) * g[j + 1], ip); if (b[j]) db = j; }
    b[6] = 0; a[7] = b[7] = 0;
    if (db < 0) return da;   /* g' = 0 mod p: g is a p-th power (p <= 3 only) or constant -> every root repeated */
    while (db >= 0) {
        /* a = a mod b */
        uint64_t ib = inv_mod_p(b[db], ip);
        while (da >= db) {
            uint64_t q = fmod_p(a[da] * ib, ip);
            for (int j = 0; j <= db; j++) { uint64_t t = fmod_p(q * b[j], ip); a[da - db + j] = fmod_p(a[da - db + j] + primes[ip] - t, ip); }
            while (da >= 0 && a[da] == 0) da--;
        }
        /* swap */
        uint64_t tmp[8]; memcpy(tmp, a, sizeof tmp); memcpy(a, b, sizeof tmp); memcpy(b, tmp, sizeof tmp); int td = da; da = db; db = td;
    }
    return da;
}
/* score of the binary sextic f(t,s) = sum out[m] t^(6-m) s^m: points over F_p of y^2 = f: x = t/s */
static int count_pts(int ip, const uint64_t f[7], int *bad) {
    uint64_t p = primes[ip]; int8_t *L = leg[ip];
    /* as polynomial in x = t (s = 1): coefficient of x^(6-m) is f[m] -> g[d] with d = 6-m */
    uint64_t g[7]; for (int m = 0; m < 7; m++) g[6-m] = f[m];
    /* bad reduction: degree drop by >= 2 (repeated root at infinity), or deg gcd(g, g') >= 1 over F_p (a repeated factor,
     * of any degree) */
    if (g[6] == 0 && g[5] == 0) { *bad = 1; return 0; }
    if (gcd_deg_mod_p(g, ip) != 0) { *bad = 1; return 0; }
    int n = 0;
    for (uint64_t x = 0; x < p; x++) {
        uint64_t r = g[6];
        for (int j = 5; j >= 0; j--) r = fmod_p(r * x + g[j], ip);
        n += 1 + L[r];
    }
    n += g[6] ? 1 + L[g[6]] : 1;
    *bad = 0;
    return n;
}

/* exact binary sextic f(t,s) = S(tP + sQ), out[m] = coefficient of t^(6-m) s^m.  Horner scheme in x0 and x1:
 *   S = sum_i x0^i G_i(x1, x2),  G_i = sum_j x1^j S_{i,j,6-i-j} x2^(6-i-j),
 * with x2 = P2 t + Q2 s expanded once (binomial powers pw2[e][m]); each Horner step multiplies a binary form by the
 * linear form (P t + Q s) (2 mpz_mul per coefficient).  ~300 mpz operations per line instead of ~1400, and all
 * temporaries are allocated once per thread. */
static void restrict_exact(const int64_t P[3], const int64_t Q[3], mpz_t out[7]) {
    static __thread int init = 0;
    static __thread mpz_t pw2[7][7], G[7], H[7], t1;
    if (!init) { for (int e = 0; e < 7; e++) for (int m = 0; m < 7; m++) mpz_init(pw2[e][m]); for (int m = 0; m < 7; m++) { mpz_init(G[m]); mpz_init(H[m]); } mpz_init(t1); init = 1; }
    /* (P2 t + Q2 s)^e */
    mpz_set_ui(pw2[0][0], 1);
    for (int e = 1; e <= 6; e++) for (int m = 0; m <= e; m++) {
        mpz_set_ui(pw2[e][m], 0);
        if (m < e) { mpz_mul_si(t1, pw2[e-1][m], P[2]); mpz_add(pw2[e][m], pw2[e][m], t1); }
        if (m > 0) { mpz_mul_si(t1, pw2[e-1][m-1], Q[2]); mpz_add(pw2[e][m], pw2[e][m], t1); }
    }
    for (int m = 0; m < 7; m++) mpz_set_ui(out[m], 0);
    int n = 0;   /* S[n] runs over (i, j) in the linesearch order: i = 6..0, j = 6-i..0 */
    for (int i = 6; i >= 0; i--) {
        /* G_i = sum_j x1^j S_{ij} x2^(6-i-j) by Horner in x1 = P1 t + Q1 s, starting from j = 6-i (degree 0 in x2) */
        int dG = -1;   /* current degree of G (as a binary form in t, s) */
        for (int j = 6 - i; j >= 0; j--) {
            int k = 6 - i - j; int64_t c = S[n++];
            /* G <- (P1 t + Q1 s) * G + c * pw2[k] */
            if (dG >= 0) {
                for (int m = 0; m <= dG + 1; m++) {
                    mpz_set_ui(H[m], 0);
                    if (m <= dG) { mpz_mul_si(t1, G[m], P[1]); mpz_add(H[m], H[m], t1); }
                    if (m > 0) { mpz_mul_si(t1, G[m-1], Q[1]); mpz_add(H[m], H[m], t1); }
                }
                dG++;
                for (int m = 0; m <= dG; m++) mpz_swap(G[m], H[m]);
            } else { for (int m = 0; m <= k; m++) mpz_set_ui(G[m], 0); dG = k; }
            if (c) for (int m = 0; m <= k; m++) { mpz_mul_si(t1, pw2[k][m], c); mpz_add(G[m], G[m], t1); }
        }
        /* out <- (P0 t + Q0 s) * out + G_i   (out has degree 6 - i - 1 before, 6 - i after) */
        if (i < 6) {
            int d = 5 - i;
            for (int m = 0; m <= d + 1; m++) {
                mpz_set_ui(H[m], 0);
                if (m <= d) { mpz_mul_si(t1, out[m], P[0]); mpz_add(H[m], H[m], t1); }
                if (m > 0) { mpz_mul_si(t1, out[m-1], Q[0]); mpz_add(H[m], H[m], t1); }
            }
            for (int m = 0; m <= d + 1; m++) mpz_swap(out[m], H[m]);
        }
        for (int m = 0; m <= dG; m++) mpz_add(out[m], out[m], G[m]);
    }
}

static void print_exact(int64_t A, int64_t B, int64_t C) {
    static __thread int init = 0; static __thread mpz_t out[7];
    if (!init) { for (int m = 0; m < 7; m++) mpz_init(out[m]); init = 1; }
    int64_t P[3], Q[3]; line_basis(A, B, C, P, Q);
    restrict_exact(P, Q, out);
    printf("%ld %ld %ld", A, B, C);
    for (int d = 0; d <= 6; d++) gmp_printf(" %Zd", out[6-d]);
    printf("\n");
}

typedef struct { double score, nagao; int64_t A, B, C; } res_t;
static int cmp_res(const void *a, const void *b) { const res_t *x = a, *y = b; return (x->score < y->score) - (x->score > y->score); }

/* Codex round 48: integer options/arguments are validated as long long against their range BEFORE any narrowing cast */
static long long chk_ll(const char *what, const char *s, long long lo, long long hi) { char *end = NULL; errno = 0; long long v = strtoll(s, &end, 10); if (errno != 0 || end == s || *end != 0) { fprintf(stderr, "%s: invalid integer '%s'\n", what, s); exit(1); } if (v < lo || v > hi) { fprintf(stderr, "%s: %lld is outside [%lld, %lld]\n", what, v, lo, hi); exit(1); } return v; }
static int arg_int(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return (int)chk_ll(opt, argv[++*i], lo, hi); }
static __attribute__((unused)) long long arg_llr(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return chk_ll(opt, argv[++*i], lo, hi); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s sextic.txt L [-B 200] [-K 1000] [-t threads] [-through a b c] [-exact A B C] [-B1 60 -frac 0.2]\n", argv[0]); return 1; }
    FILE *fp = fopen(argv[1], "r"); if (!fp) { perror("sextic"); return 1; }
    { int n = 0; char buf[65536]; size_t len = fread(buf, 1, sizeof buf - 1, fp); buf[len] = 0; fclose(fp);
      char *tok = strtok(buf, " ,\t\r\n[]");   /* Codex round 49: exactly 28 complete int64 coefficients, nothing else */
      while (tok) { if (n == 28) { fprintf(stderr, "%s: more than 28 coefficients\n", argv[1]); return 1; } char *end = NULL; errno = 0; long long v = strtoll(tok, &end, 10); if (errno != 0 || end == tok || *end != 0) { fprintf(stderr, "%s: coefficient %d is not an int64 integer: '%s'\n", argv[1], n + 1, tok); return 1; } S[n++] = v; tok = strtok(NULL, " ,\t\r\n[]"); }
      if (n != 28) { fprintf(stderr, "need 28 coefficients, got %d\n", n); return 1; } }
    { int n = 0; for (int i = 6; i >= 0; i--) for (int j = 6 - i; j >= 0; j--) { SI[n] = i; SJ[n] = j; SK[n] = 6 - i - j; n++; } }
    int64_t L = chk_ll("L", argv[2], 1, 1000000000); int B = 200, K = 1000, nthreads = omp_get_max_threads(); int64_t th[3] = {0,0,0}; int through = 0; int exact = 0; int64_t ex[3]; int pencil = 0; int64_t pc[3] = {0,0,0}; int enum_only = 0, exactbatch = 0;
    int B1 = 60; double frac = 0.25;
    for (int i = 3; i < argc; i++) {
        /* Codex round 48: every option argument is checked (presence, syntax, range); unknown options are errors — a
         * misspelt -nosym used to be ignored silently, i.e. a sweep would run WITH the symmetry filter */
        if (!strcmp(argv[i], "-B")) B = arg_int(argc, argv, &i, 3, 313); else if (!strcmp(argv[i], "-K")) K = arg_int(argc, argv, &i, 1, 100000000); else if (!strcmp(argv[i], "-t")) nthreads = arg_int(argc, argv, &i, 1, 4096);
        else if (!strcmp(argv[i], "-through")) { through = 1; for (int k = 0; k < 3; k++) th[k] = arg_llr(argc, argv, &i, -(1LL << 40), 1LL << 40); }
        else if (!strcmp(argv[i], "-nosym")) use_sym = 0; else if (!strcmp(argv[i], "-enum")) enum_only = 1; else if (!strcmp(argv[i], "-exactbatch")) exactbatch = 1;
        else if (!strcmp(argv[i], "-pencil")) { pencil = 1; for (int k = 0; k < 3; k++) pc[k] = arg_llr(argc, argv, &i, -(1LL << 40), 1LL << 40); }
        else if (!strcmp(argv[i], "-exact")) { exact = 1; for (int k = 0; k < 3; k++) ex[k] = arg_llr(argc, argv, &i, -(1LL << 40), 1LL << 40); }
        else if (!strcmp(argv[i], "-B1")) B1 = arg_int(argc, argv, &i, 3, 313);
        else if (!strcmp(argv[i], "-frac")) { if (i + 1 >= argc) { fprintf(stderr, "option -frac needs an argument\n"); return 1; } char *end = NULL; errno = 0; frac = strtod(argv[++i], &end); if (errno != 0 || end == argv[i] || *end != 0 || !(frac > 0 && frac <= 1)) { fprintf(stderr, "-frac: give a number in (0, 1], not '%s'\n", argv[i]); return 1; } }
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (exact) {
        int64_t P[3], Q[3]; line_basis(ex[0], ex[1], ex[2], P, Q);
        /* print model only (ascending in x = t/s) */
        mpz_t out[7]; for (int m = 0; m < 7; m++) mpz_init(out[m]);
        restrict_exact(P, Q, out);
        for (int d = 0; d <= 6; d++) gmp_printf("%Zd%s", out[6-d], d < 6 ? " " : "\n");
        fprintf(stderr, "basis P=(%ld,%ld,%ld) Q=(%ld,%ld,%ld); x = t/s\n", P[0], P[1], P[2], Q[0], Q[1], Q[2]);
        return 0;
    }
    if (exactbatch) {
        int64_t A, B, C;
        while (scanf("%ld %ld %ld", &A, &B, &C) == 3) print_exact(A, B, C);
        return 0;
    }
    if (enum_only) {
        if (pencil) {
            int64_t PB[3], QB[3]; line_basis(pc[0], pc[1], pc[2], PB, QB);
            for (int64_t u = 0; u <= L; u++) for (int64_t v = (u == 0 ? 1 : -L); v <= L; v++) {
                int64_t A = u * PB[0] + v * QB[0], Bq = u * PB[1] + v * QB[1], C = u * PB[2] + v * QB[2];
                if (A < 0 || (A == 0 && Bq < 0) || (A == 0 && Bq == 0 && C < 0)) { A = -A; Bq = -Bq; C = -C; }
                if (gcd3(A, Bq, C) != 1) continue;
                if (use_sym && A > 0 && Bq < 0) continue;
                printf("%ld %ld %ld\n", A, Bq, C);
            }
        } else {
            for (int64_t A = 0; A <= L; A++) for (int64_t Bq = (A == 0 ? 0 : (use_sym ? 0 : -L)); Bq <= L; Bq++) for (int64_t C = ((A == 0 && Bq == 0) ? 1 : -L); C <= L; C++) {
                /* A > 0: B >= 0 by symmetry; A == 0: (0,B,C) with B >= 0 (sign) */
                if (gcd3(A, Bq, C) != 1) continue;
                if (A == 0 && Bq == 0) { if (C != 1) continue; }
                printf("%ld %ld %ld\n", A, Bq, C);
            }
        }
        return 0;
    }
    /* primes */
    if (B < 3 || B > 313) { fprintf(stderr, "-B must be between 3 and 313 (at most 64 odd primes fit the fixed tables; Codex round 16)\n"); return 1; }
    { char *sv = calloc(B + 1, 1); for (int p = 3; p <= B; p += 2) if (!sv[p]) {   /* odd primes only (the step-1 loop let 4 through) */ if (nprimes == 64) { fprintf(stderr, "too many primes\n"); return 1; } primes[nprimes++] = p; for (int q = 2*p; q <= B; q += p) sv[q] = 1; } free(sv); }
    for (int ip = 0; ip < nprimes; ip++) { int p = primes[ip]; leg[ip] = malloc(p); memset(leg[ip], -1, p); leg[ip][0] = 0; for (int64_t x = 1; x < p; x++) leg[ip][(x*x) % p] = 1;
        for (int n = 0; n < 28; n++) Smod[ip][n] = ((S[n] % p) + p) % p; pmagic[ip] = UINT64_C(0xFFFFFFFFFFFFFFFF) / p + 1; }
    int np1 = 0; while (np1 < nprimes && primes[np1] <= B1) np1++;
    omp_set_num_threads(nthreads);
    /* stage 1 threshold: estimated by sampling? simpler: keep everything with score1 >= thr1 where thr1 chosen so that ~frac survive; do a calibration on random lines */
    double thr1 = -1e9;
    {
        int nsamp = 20000; double *sc = malloc(nsamp * sizeof(double)); uint64_t rng = 12345;
        for (int s = 0; s < nsamp; s++) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; int64_t A = rng % (2*L+1) - L, Bq = (rng >> 20) % (2*L+1) - L, C = (rng >> 40) % (2*L+1) - L;
            if (gcd3(A, Bq, C) != 1 || (A == 0 && Bq == 0 && C == 0)) { sc[s] = -1e9; continue; }
            int64_t P[3], Q[3]; line_basis(A, Bq, C, P, Q); double v = 0; int nb = 0;
            for (int ip = 0; ip < np1; ip++) { uint64_t f[7]; restrict_mod(ip, P, Q, f); int bad; int n = count_pts(ip, f, &bad); if (bad) { nb++; continue; } v += log((double)n / primes[ip]); }
            sc[s] = nb > np1 / 2 ? -1e9 : v; }
        /* find the (1-frac) quantile */
        double *tmp = malloc(nsamp * sizeof(double)); memcpy(tmp, sc, nsamp * sizeof(double));
        int cmpd(const void *a, const void *b) { double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y); }
        qsort(tmp, nsamp, sizeof(double), cmpd);
        thr1 = tmp[(int)((1 - frac) * nsamp)];
        fprintf(stderr, "stage-1 threshold (primes <= %d): %.3f (keeping ~%.0f%%)\n", B1, thr1, frac * 100);
        free(sc); free(tmp);
    }
    /* enumerate */
    res_t *top = calloc((size_t)K * nthreads * 2 + 16, sizeof(res_t)); int *ntop = calloc(nthreads, sizeof(int));
    uint64_t nlines = 0, nstage2 = 0;
    #pragma omp parallel reduction(+:nlines, nstage2)
    {
        int tid = omp_get_thread_num(); res_t *mine = top + (size_t)tid * 2 * K; int n = 0; double worst = -1e9;
        if (pencil) {
            int64_t PB[3], QB[3]; line_basis(pc[0], pc[1], pc[2], PB, QB);
            #pragma omp for schedule(dynamic, 4)
            for (int64_t u = 0; u <= L; u++) {
                for (int64_t v = (u == 0 ? 1 : -L); v <= L; v++) {
                    int64_t A = u * PB[0] + v * QB[0], Bq = u * PB[1] + v * QB[1], C = u * PB[2] + v * QB[2];
                    if (A < 0 || (A == 0 && Bq < 0) || (A == 0 && Bq == 0 && C < 0)) { A = -A; Bq = -Bq; C = -C; }
                if (gcd3(A, Bq, C) != 1) continue;
                if (use_sym && A > 0 && Bq < 0) continue;   /* symmetry x0 -> -x0 (only if the surface has it) */
                    if (through && A * th[0] + Bq * th[1] + C * th[2] != 0) continue;
                    nlines++;
                    int64_t P[3], Q[3]; line_basis(A, Bq, C, P, Q);
                    double v = 0, ng = 0; int nb = 0;
                    for (int ip = 0; ip < np1; ip++) { uint64_t f[7]; restrict_mod(ip, P, Q, f); int bad; int cnt = count_pts(ip, f, &bad); if (bad) { nb++; continue; } v += log((double)cnt / primes[ip]); ng += (cnt - primes[ip] - 1.0) * log((double)primes[ip]) / primes[ip]; }
                    if (nb > np1 / 2 || v < thr1) continue;
                    nstage2++;
                    for (int ip = np1; ip < nprimes; ip++) { uint64_t f[7]; restrict_mod(ip, P, Q, f); int bad; int cnt = count_pts(ip, f, &bad); if (bad) continue; v += log((double)cnt / primes[ip]); ng += (cnt - primes[ip] - 1.0) * log((double)primes[ip]) / primes[ip]; }
                    if (n < 2 * K) { mine[n].score = v; mine[n].nagao = ng; mine[n].A = A; mine[n].B = Bq; mine[n].C = C; n++; if (n == 2 * K) { qsort(mine, n, sizeof(res_t), cmp_res); n = K; worst = mine[K-1].score; } }
                    else if (v > worst) { mine[n].score = v; mine[n].nagao = ng; mine[n].A = A; mine[n].B = Bq; mine[n].C = C; n++; if (n == 2 * K) { qsort(mine, n, sizeof(res_t), cmp_res); n = K; worst = mine[K-1].score; } }

                }
            }
        } else {
        #pragma omp for schedule(dynamic, 4)
        for (int64_t A = 0; A <= L; A++) {
            for (int64_t Bq = (A == 0 ? 0 : -L); Bq <= L; Bq++) for (int64_t C = ((A == 0 && Bq == 0) ? 1 : -L); C <= L; C++) {
                if (gcd3(A, Bq, C) != 1) continue;
                if (use_sym && A > 0 && Bq < 0) continue;   /* symmetry x0 -> -x0 (only if the surface has it) */
                if (through && A * th[0] + Bq * th[1] + C * th[2] != 0) continue;
                nlines++;
                int64_t P[3], Q[3]; line_basis(A, Bq, C, P, Q);
                double v = 0, ng = 0; int nb = 0;
                for (int ip = 0; ip < np1; ip++) { uint64_t f[7]; restrict_mod(ip, P, Q, f); int bad; int cnt = count_pts(ip, f, &bad); if (bad) { nb++; continue; } v += log((double)cnt / primes[ip]); ng += (cnt - primes[ip] - 1.0) * log((double)primes[ip]) / primes[ip]; }
                if (nb > np1 / 2 || v < thr1) continue;
                nstage2++;
                for (int ip = np1; ip < nprimes; ip++) { uint64_t f[7]; restrict_mod(ip, P, Q, f); int bad; int cnt = count_pts(ip, f, &bad); if (bad) continue; v += log((double)cnt / primes[ip]); ng += (cnt - primes[ip] - 1.0) * log((double)primes[ip]) / primes[ip]; }
                if (n < 2 * K) { mine[n].score = v; mine[n].nagao = ng; mine[n].A = A; mine[n].B = Bq; mine[n].C = C; n++; if (n == 2 * K) { qsort(mine, n, sizeof(res_t), cmp_res); n = K; worst = mine[K-1].score; } }
                else if (v > worst) { mine[n].score = v; mine[n].nagao = ng; mine[n].A = A; mine[n].B = Bq; mine[n].C = C; n++; if (n == 2 * K) { qsort(mine, n, sizeof(res_t), cmp_res); n = K; worst = mine[K-1].score; } }
            }
        }
        }
        qsort(mine, n, sizeof(res_t), cmp_res); if (n > K) n = K; ntop[tid] = n;
    }
    /* merge */
    res_t *all = malloc((size_t)K * nthreads * sizeof(res_t)); int na = 0;
    for (int t = 0; t < nthreads; t++) { memcpy(all + na, top + (size_t)t * 2 * K, ntop[t] * sizeof(res_t)); na += ntop[t]; }
    qsort(all, na, sizeof(res_t), cmp_res); if (na > K) na = K;
    fprintf(stderr, "%lu lines, %lu passed stage 1, printing top %d\n", nlines, nstage2, na);
    for (int i = 0; i < na; i++) printf("%.4f %.4f %ld %ld %ld\n", all[i].score, all[i].nagao, all[i].A, all[i].B, all[i].C);
    return 0;
}

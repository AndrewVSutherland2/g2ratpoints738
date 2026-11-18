/*
 * g2rp - fast search for rational points on hyperelliptic curves y^2 = f(x)
 *        (deg f <= 6, coefficients f0..f6 ascending, integers of any size: GMP), height bound N:
 *        finds all (a : b : c) with x = a/c, gcd(a,c)=1, |a|,|c| <= N.
 *
 * Method: 2-D sieve over (a, c).  For each denominator c the numerators a are
 * processed in chunks of 512 (one AVX-512 register).  Sieve moduli are
 * composites M (products of prime powers, chosen by a knapsack on the
 * per-prime survival densities); for each M and each residue c mod M a
 * compact periodic bit row (M + 576 bits) holds "F(a,c) is a square mod M";
 * the 512-bit word for a chunk is cut out of the row with two unaligned
 * loads and a VBMI2 funnel shift.  Survivors go through scalar bit tables
 * for the remaining primes, a gcd test and an exact GMP square test.
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <gmp.h>
#include <omp.h>

typedef unsigned __int128 u128;

#define MAXDEG 6
#define MAXVEC 16          /* max vector moduli (two zmm of 64-bit lanes) */
#define BLOCK_CHUNKS 128
#define G 4                /* denominators processed simultaneously */
#define SMALL_N 30000      /* heights up to this get the light automatic setup (see run_main) */

static int verbose = 0;

static inline uint64_t fastmod_magic(uint32_t m) { return UINT64_C(0xFFFFFFFFFFFFFFFF) / m + 1; }
static inline uint32_t fastmod_u32(uint32_t a, uint64_t M, uint32_t m) {
    uint64_t low = M * a;
    return (uint32_t)(((u128)low * m) >> 64);
}
static int is_prime(uint32_t n) { if (n < 2) return 0; for (uint32_t d = 2; d*d <= n; d++) if (n % d == 0) return 0; return 1; }

/* ------------------------------------------------------------------ */
/* prime-power tables: ok[cbar*q + abar] */
typedef struct {
    uint32_t q, p;       /* modulus (prime power) and its prime */
    double density;
    uint8_t *ok;         /* q*q bytes */
    uint64_t *bits;      /* scalar bit table, row stride rowwords */
    uint32_t rowwords;
    uint64_t magic;
    uint32_t abase_mod;
} pptab_t;

/* Structure of the table ok[cb][ab] = "F(ab,cb) = sum f_i ab^i cb^(6-i) is a square mod q" (q = p^k), with the
 * convention ok = 0 when p | ab and p | cb (impossible for coprime (a,c)):
 *   cb a unit:  F = cb^6 f(ab/cb) and cb^6 is a unit square, so ok[cb][ab] = rowf[ab * cb^-1 mod q]
 *               where rowf[x] = "f(x) is a square mod q";
 *   p | cb, ab a unit:  F = ab^6 g(cb/ab) with g(t) = t^6 f(1/t) = sum f_i t^(6-i), so ok[cb][ab] = rowg[cb * ab^-1]
 *               where rowg[t] = "g(t) is a square mod q".
 * Hence every row is a permutation of rowf or a stepped read of rowg: the whole table costs O(q^2) byte moves (no
 * multiplications), and the exact density is O(q): (q - q/p) * (sum rowf + sum_{p | t} rowg[t]) / q^2 -- the map
 * (cb, ab) -> (cb/ab, ab) is a bijection from {p | cb, ab unit} onto {p | t, ab unit}. */
typedef struct { uint8_t *rowf, *rowg; uint32_t q; } pprows_t;

static void pp_rows(mpz_t *f, uint32_t q, pprows_t *R) {
    R->q = q; R->rowf = malloc(q); R->rowg = malloc(q);
    uint8_t *sq = calloc(q, 1);
    for (uint64_t x = 0; x < q; x++) sq[(x*x) % q] = 1;
    uint64_t magic = fastmod_magic(q);
    uint64_t fm[MAXDEG+1];
    for (int i = 0; i <= MAXDEG; i++) fm[i] = mpz_fdiv_ui(f[i], q);   /* non-negative residue of a coefficient of any size */
    for (uint32_t x = 0; x < q; x++) {
        uint64_t r = fm[MAXDEG], s = fm[0];
        for (int i = MAXDEG-1; i >= 0; i--) { r = r * x + fm[i]; uint64_t low = magic * r; r = (uint64_t)(((u128)low * q) >> 64); }
        for (int i = 1; i <= MAXDEG; i++) { s = s * x + fm[i]; uint64_t low = magic * s; s = (uint64_t)(((u128)low * q) >> 64); }
        R->rowf[x] = sq[r]; R->rowg[x] = sq[s];
    }
    free(sq);
}
static double pp_density(const pprows_t *R, uint32_t p) {
    uint32_t q = R->q; uint64_t sf = 0, sg = 0;
    for (uint32_t x = 0; x < q; x++) sf += R->rowf[x];
    for (uint32_t t = 0; t < q; t += p) sg += R->rowg[t];
    return (double)((uint64_t)(q - q / p) * (sf + sg)) / ((double)q * q);
}
static void pp_rows_free(pprows_t *R) { free(R->rowf); free(R->rowg); }

/* inverses of the units mod q (0 for non-units), by the extended Euclidean algorithm */
static uint32_t *inverse_table(uint32_t q, uint32_t p) {
    uint32_t *inv = calloc(q, sizeof(uint32_t));
    for (uint32_t u = 1; u < q; u++) {
        if (u % p == 0) continue;
        int64_t r0 = q, r1 = u, t0 = 0, t1 = 1;
        while (r1) { int64_t qq = r0 / r1, t = r0 - qq * r1; r0 = r1; r1 = t; t = t0 - qq * t1; t0 = t1; t1 = t; }
        inv[u] = (uint32_t)(t0 < 0 ? t0 + (int64_t)q : t0);
    }
    return inv;
}

/* build the q*q byte table ok[cb][ab] (see above); parallel over rows when the table is large */
static void build_pptab(pptab_t *T, uint32_t q, uint32_t p, mpz_t *f) {
    T->q = q; T->p = p; T->magic = fastmod_magic(q);
    pprows_t R; pp_rows(f, q, &R);
    T->density = pp_density(&R, p);
    T->ok = malloc((size_t)q * q);
    T->rowwords = (q + 63) / 64;
    T->bits = NULL;
    int big = (uint64_t)q * q > (1u << 22);
    /* unit rows: ok[cb][x * cb] = rowf[x] */
    #pragma omp parallel for schedule(dynamic, 16) if(big)
    for (uint32_t cb = 0; cb < q; cb++) {
        if (cb % p == 0) continue;
        uint8_t *ok = T->ok + (size_t)cb * q;
        uint32_t y = 0;
        for (uint32_t x = 0; x < q; x++) { ok[y] = R.rowf[x]; y += cb; if (y >= q) y -= q; }
    }
    /* rows with p | cb, in blocks of 64 columns: column ab (a unit) reads rowg at cb * ab^-1, which steps by
     * p * ab^-1 from one row to the next */
    uint32_t *inv = inverse_table(q, p);
    #pragma omp parallel for schedule(dynamic, 1) if(big)
    for (uint32_t ab0 = 0; ab0 < q; ab0 += 64) {
        uint32_t nb = q - ab0 < 64 ? q - ab0 : 64;
        uint32_t x[64], st[64]; uint8_t unit[64];
        for (uint32_t l = 0; l < nb; l++) { uint32_t ab = ab0 + l; unit[l] = (ab % p != 0); st[l] = unit[l] ? (uint32_t)(((uint64_t)p * inv[ab]) % q) : 0; x[l] = 0; }
        for (uint32_t cb = 0; cb < q; cb += p) {
            uint8_t *ok = T->ok + (size_t)cb * q + ab0;
            for (uint32_t l = 0; l < nb; l++) { ok[l] = unit[l] ? R.rowg[x[l]] : 0; x[l] += st[l]; if (x[l] >= q) x[l] -= q; }
        }
    }
    free(inv);
    pp_rows_free(&R);
}
/* pack a byte row (one byte per bit, values 0/1) into bits */
static inline void pack_bits(const uint8_t *ok, uint32_t n, uint64_t *w) {
    uint32_t j = 0;
    for (; j + 64 <= n; j += 64) w[j >> 6] = _mm512_cmpneq_epi8_mask(_mm512_loadu_si512(ok + j), _mm512_setzero_si512());
    if (j < n) w[j >> 6] = _mm512_mask_cmpneq_epi8_mask((__mmask64)((UINT64_C(1) << (n - j)) - 1), _mm512_maskz_loadu_epi8((__mmask64)((UINT64_C(1) << (n - j)) - 1), ok + j), _mm512_setzero_si512());
}
/* scalar bit table (used for the scalar moduli only) */
static void build_bits(pptab_t *T) {
    if (T->bits) return;
    T->bits = calloc((size_t)T->q * T->rowwords, 8);
    for (uint32_t cb = 0; cb < T->q; cb++) pack_bits(T->ok + (size_t)cb * T->q, T->q, T->bits + (size_t)cb * T->rowwords);
}

/* ------------------------------------------------------------------ */
/* vector modulus: composite M = product of prime powers; compact periodic rows */
typedef struct {
    uint32_t M;
    int nf; pptab_t *fac[8];
    double density;
    size_t rowbytes;     /* bytes per row (common stride ROWBYTES) */
    uint8_t *rows;       /* M + G rows (rows M..M+G-1 duplicate rows 0..G-1) */
    uint8_t *rowempty;
    uint64_t magic;
} vmod_t;

/* periodic extension of a packed row: bits j >= q (up to words*64) copied from bit j - q, word-wise */
static void extend_periodic(uint64_t *w, uint32_t q, size_t words) {
    size_t nb = words * 64; uint32_t Q = q;
    if (Q < 64) {   /* bit by bit up to a period Q' = q*ceil(64/q) >= 64, then word-wise with period Q' */
        Q = q * ((64 + q - 1) / q);
        for (uint32_t j = q; j < Q && j < nb; j++) w[j >> 6] |= ((w[(j - q) >> 6] >> ((j - q) & 63)) & 1) << (j & 63);
    }
    if (Q >= nb) return;
    size_t wi0 = Q >> 6; uint32_t r = Q & 63;
    if (r) w[wi0] |= w[0] << r; else w[wi0] = w[0];
    for (size_t wi = wi0 + 1; wi < words; wi++) {
        size_t s = 64 * wi - Q; uint32_t sh = (uint32_t)(s & 63); const uint64_t *src = w + (s >> 6);
        w[wi] = sh ? (src[0] >> sh) | (src[1] << (64 - sh)) : src[0];
    }
}

static void build_vmod(vmod_t *V, uint32_t M, int nf, pptab_t **fac, size_t rowbytes) {
    V->M = M; V->nf = nf; V->magic = fastmod_magic(M);
    for (int i = 0; i < nf; i++) V->fac[i] = fac[i];
    V->density = 1; for (int i = 0; i < nf; i++) V->density *= fac[i]->density;
    size_t words = rowbytes / 8;
    V->rowbytes = rowbytes;
    V->rows = aligned_alloc(64, (size_t)(M + G) * V->rowbytes + 64);
    V->rowempty = calloc(M, 1);
    if (!V->rows) { fprintf(stderr, "alloc failed for modulus %u\n", M); exit(1); }
    int big = (size_t)M * words > (1u << 18);
    /* tile each factor's rows to `words` words: tiled[i][r] has bit j = ok_q[r][j mod q] */
    uint64_t *tiled[8];
    for (int i = 0; i < nf; i++) {
        uint32_t q = fac[i]->q;
        tiled[i] = aligned_alloc(64, (size_t)q * words * 8);
        #pragma omp parallel for schedule(static) if(big)
        for (uint32_t r = 0; r < q; r++) {
            uint64_t *w = tiled[i] + (size_t)r * words;
            memset(w, 0, words * 8);
            pack_bits(fac[i]->ok + (size_t)r * q, q, w);
            extend_periodic(w, q, words);
        }
    }
    #pragma omp parallel for schedule(static) if(big)
    for (uint32_t cb = 0; cb < M; cb++) {
        uint64_t *w = (uint64_t*)(V->rows + (size_t)cb * V->rowbytes);
        const uint64_t *src = tiled[0] + (size_t)(cb % fac[0]->q) * words;
        memcpy(w, src, words * 8);
        for (int i = 1; i < nf; i++) {
            const uint64_t *t = tiled[i] + (size_t)(cb % fac[i]->q) * words;
            for (size_t k = 0; k < words; k++) w[k] &= t[k];
        }
        uint64_t any = 0; for (size_t k = 0; k < words; k++) any |= w[k];
        V->rowempty[cb] = (any == 0);
    }
    for (uint32_t g = 0; g < G; g++) memcpy(V->rows + (size_t)(M + g) * V->rowbytes, V->rows + (size_t)g * V->rowbytes, V->rowbytes);
    for (int i = 0; i < nf; i++) free(tiled[i]);
}

/* ------------------------------------------------------------------ */
/* real roots of f in long double arithmetic on the scaled coefficients of scaled_coeffs (the top 64 bits of every
 * coefficient, whatever its size).  real_roots returns -1 when the arithmetic overflowed (an absurd spread of coefficient
 * sizes: the root bound B is ~2^E for a leading coefficient E bits below the largest one, and f(B) ~ B^5 exceeds the
 * long double range once E > ~3300); the caller then searches without the real-root restriction, which is only a
 * pruning of the a-ranges */
static long double horner(long double *cf, int deg, long double x) { long double r = cf[deg]; for (int i = deg-1; i >= 0; i--) r = r*x + cf[i]; return r; }
static long double bisect(long double *cf, int deg, long double lo, long double hi) {
    long double flo = horner(cf, deg, lo);
    /* the bracket may span thousands of binary orders of magnitude (huge coefficient spreads: e.g. [-2^994, 0] with the
     * root near -2^198), so iterate to a tolerance -- 1e-18 relative or 1e-30 absolute (the intervals are padded by
     * 1e-9 (|x|+1) anyway) -- and not for a fixed count; every step halves the width, so ~16500 steps is the maximum */
    for (int it = 0; it < 40000; it++) {
        long double mid = 0.5L*(lo+hi);
        if (mid <= lo || mid >= hi) break;   /* no representable midpoint left */
        long double fm = horner(cf, deg, mid);
        if ((fm < 0) == (flo < 0)) { lo = mid; flo = fm; } else hi = mid;
        if (hi - lo <= 1e-18L * (fabsl(lo) + fabsl(hi)) || hi - lo <= 1e-30L) break;
    }
    return 0.5L*(lo+hi);
}
static int real_roots(long double *cf, int deg, long double *roots) {
    while (deg > 0 && cf[deg] == 0) deg--;
    if (deg <= 0) return 0;
    if (deg == 1) { roots[0] = -cf[0]/cf[1]; return isfinite(roots[0]) ? 1 : -1; }
    long double dcf[MAXDEG+1];
    for (int i = 1; i <= deg; i++) dcf[i-1] = i * cf[i];
    long double droots[MAXDEG+1];
    int nd = real_roots(dcf, deg-1, droots);
    if (nd < 0) return -1;
    long double B = 0;
    for (int i = 0; i < deg; i++) { long double t = fabsl(cf[i]/cf[deg]); if (t > B) B = t; }
    B += 1;
    long double pts[MAXDEG+3]; int np = 0;
    pts[np++] = -B;
    for (int i = 0; i < nd; i++) if (droots[i] > -B && droots[i] < B) pts[np++] = droots[i];
    pts[np++] = B;
    int nr = 0;
    for (int i = 0; i + 1 < np; i++) {
        long double flo = horner(cf, deg, pts[i]), fhi = horner(cf, deg, pts[i+1]);
        if (!isfinite(flo) || !isfinite(fhi)) return -1;
        if (flo == 0) { roots[nr++] = pts[i]; continue; }
        if ((flo < 0) != (fhi < 0)) roots[nr++] = bisect(cf, deg, pts[i], pts[i+1]);
    }
    if (horner(cf, deg, pts[np-1]) == 0) roots[nr++] = pts[np-1];
    return nr;
}
/* the coefficients as long doubles times a common 2^-E, E = bits of the largest |f_i| (the roots are unchanged): every
 * |cf_i| < 1, the top 64 bits of each coefficient are exact (all of them when |f_i| < 2^64), and nothing over- or
 * underflows unless the coefficient sizes span more than ~16000 bits */
static void scaled_coeffs(mpz_t *f, long double *cf) {
    size_t E = 0;
    for (int i = 0; i <= MAXDEG; i++) if (mpz_sgn(f[i]) && mpz_sizeinbase(f[i], 2) > E) E = mpz_sizeinbase(f[i], 2);
    mpz_t t; mpz_init(t);
    for (int i = 0; i <= MAXDEG; i++) {
        if (mpz_sgn(f[i]) == 0) { cf[i] = 0; continue; }
        size_t b = mpz_sizeinbase(f[i], 2); long sh = b > 64 ? (long)(b - 64) : 0;
        mpz_tdiv_q_2exp(t, f[i], sh);                                  /* |t| < 2^64, the sign of f[i] */
        long double d = (long double)mpz_get_ui(t); if (mpz_sgn(t) < 0) d = -d;
        cf[i] = ldexpl(d, (int)(sh - (long)E));
    }
    mpz_clear(t);
}
/* a coefficient token: an optional sign and decimal digits, nothing else; any size (GMP) */
static int parse_coeff(mpz_t z, const char *tok) {
    const char *s = tok + (*tok == '+' || *tok == '-');
    if (!*s) return 0;
    for (const char *p = s; *p; p++) if (*p < '0' || *p > '9') return 0;
    return mpz_set_str(z, *tok == '+' ? tok + 1 : tok, 10) == 0;
}
static void clear_coeffs(mpz_t *f) { for (int i = 0; i <= MAXDEG; i++) mpz_clear(f[i]); }

/* ------------------------------------------------------------------ */
typedef struct { int64_t a, c; mpz_t b; } point_t;
typedef struct {
    point_t *pts; size_t n, cap;
    uint64_t cand_vec, cand_scalar, cand_gmp;
    uint64_t chunks, chunks_s2, chunks_s3, c_skipped;
    uint64_t t_s1, t_s2, t_s3, t_cand, t_setup;
} thread_result_t;
static void add_point(thread_result_t *R, int64_t a, int64_t c, mpz_t b) {
    if (R->n == R->cap) { R->cap = R->cap ? 2*R->cap : 64; R->pts = realloc(R->pts, R->cap * sizeof(point_t)); }
    R->pts[R->n].a = a; R->pts[R->n].c = c; mpz_init_set(R->pts[R->n].b, b); R->n++;
}

typedef struct {
    mpz_t f[MAXDEG+1];   /* coefficients of any size */
    int64_t N;
    int64_t a_base;
    int nint; double ilo[8], ihi[8];
    /* power of 2 */
    uint32_t m2; uint8_t *W2; uint8_t *rowempty2; double d2;
    /* vector moduli */
    int nvec; vmod_t vec[MAXVEC]; size_t rowbytes;
    int K1, K2;
    int odd; int ne; int64_t elist[4096];
    /* scalar prime powers */
    int nsca; pptab_t **sca;
} config_t;

static inline uint64_t gcd64(uint64_t a, uint64_t b) {
    if (a == 0) return b;
    if (b == 0) return a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do { b >>= __builtin_ctzll(b); if (a > b) { uint64_t t = a; a = b; b = t; } b -= a; } while (b);
    return a << sh;
}

typedef struct { uint32_t *cmod; mpz_t F, t1, t2; uint8_t *scratch; } cand_ctx_t;   /* per-thread; scratch = nvec*G rows for quintic mode */

static void check_candidate(config_t *C, thread_result_t *R, cand_ctx_t *X, int64_t a, int64_t c) {
    uint32_t aa = (uint32_t)(a - C->a_base);
    for (int j = 0; j < C->nsca; j++) {
        pptab_t *T = C->sca[j];
        uint32_t am = fastmod_u32(aa, T->magic, T->q) + T->abase_mod; if (am >= T->q) am -= T->q;
        uint32_t cm = X->cmod[j];
        if (!((T->bits[(size_t)cm * T->rowwords + am/64] >> (am % 64)) & 1)) return;
    }
    R->cand_scalar++;
    if (gcd64((uint64_t)(a < 0 ? -a : a), (uint64_t)c) != 1) return;
    R->cand_gmp++;
    mpz_set(X->F, C->f[MAXDEG]);
    mpz_set_si(X->t1, 1);
    for (int i = MAXDEG-1; i >= 0; i--) {
        mpz_mul_si(X->F, X->F, a);
        mpz_mul_si(X->t1, X->t1, c);
        mpz_mul(X->t2, C->f[i], X->t1);
        mpz_add(X->F, X->F, X->t2);
    }
    if (mpz_sgn(X->F) < 0) return;
    if (mpz_perfect_square_p(X->F)) { mpz_sqrt(X->t1, X->F); add_point(R, a, c, X->t1); }
}

/* ------------------------------------------------------------------ */
/* stage buffers: surviving chunk words with a tag (block-relative chunk index) * G + g; the phases of a block's
 * chunks are stored once per chunk (they do not depend on g), not once per surviving word */
typedef struct { __m512i phi0, phi1; } phi_t;
typedef struct { __m512i *acc1, *acc2; uint32_t *tag1, *tag2; phi_t *phitab; } bufs_t;

typedef struct {
    __m512i base0, base1, inc0, inc1, mod0, mod1;
    size_t rowbytes;
} vstate_t;

/* apply moduli j0..j1-1 for denominator slot g */
static inline void apply_moduli(__m512i *acc, __m512i phi0, __m512i phi1, const vstate_t *S, int g, int j0, int j1) {
    uint64_t addr[16] __attribute__((aligned(64)));
    uint64_t ph[16] __attribute__((aligned(64)));
    __m512i mask = _mm512_set1_epi64(~7ULL);
    __m512i goff = _mm512_set1_epi64((uint64_t)g * S->rowbytes);
    _mm512_store_si512(addr, _mm512_add_epi64(_mm512_add_epi64(S->base0, goff), _mm512_and_si512(_mm512_srli_epi64(phi0, 3), mask)));
    _mm512_store_si512(addr+8, _mm512_add_epi64(_mm512_add_epi64(S->base1, goff), _mm512_and_si512(_mm512_srli_epi64(phi1, 3), mask)));
    _mm512_store_si512(ph, phi0); _mm512_store_si512(ph+8, phi1);
    __m512i a = *acc;
    for (int j = j0; j < j1; j++) {
        const uint8_t *r = (const uint8_t*)addr[j];
        __m512i w = _mm512_shrdv_epi64(_mm512_loadu_si512(r), _mm512_loadu_si512(r + 8), _mm512_set1_epi64(ph[j]));
        a = _mm512_and_si512(a, w);
    }
    *acc = a;
}

/* stage 1: chunks ib..ie for G consecutive denominators; moduli 0..K-1 */
#define DEFINE_STAGE1(K) \
static size_t stage1_##K(const uint8_t **prow, uint32_t pmask, const vstate_t *S, __m512i phi0, __m512i phi1, int64_t ib, int64_t ie, bufs_t *B) { \
    uint64_t addr[16] __attribute__((aligned(64))); uint64_t ph[16] __attribute__((aligned(64))); \
    __m512i mask = _mm512_set1_epi64(~7ULL); \
    const size_t RB = S->rowbytes; \
    size_t n = 0; \
    for (int64_t i = ib; i <= ie; i++) { \
        _mm512_store_si512(addr, _mm512_add_epi64(S->base0, _mm512_and_si512(_mm512_srli_epi64(phi0, 3), mask))); \
        if (K > 8) _mm512_store_si512(addr+8, _mm512_add_epi64(S->base1, _mm512_and_si512(_mm512_srli_epi64(phi1, 3), mask))); \
        _mm512_store_si512(ph, phi0); if (K > 8) _mm512_store_si512(ph+8, phi1); \
        __m512i acc[G]; \
        for (int g = 0; g < G; g++) acc[g] = _mm512_load_si512(prow[g] + (size_t)(i & pmask) * 64); \
        for (int j = 0; j < K; j++) { \
            const uint8_t *r = (const uint8_t*)addr[j]; \
            __m512i cnt = _mm512_set1_epi64(ph[j]); \
            for (int g = 0; g < G; g++) \
                acc[g] = _mm512_and_si512(acc[g], _mm512_shrdv_epi64(_mm512_loadu_si512(r + g*RB), _mm512_loadu_si512(r + g*RB + 8), cnt)); \
        } \
        _mm512_store_si512(&B->phitab[i - ib].phi0, phi0); _mm512_store_si512(&B->phitab[i - ib].phi1, phi1); \
        for (int g = 0; g < G; g++) { \
            _mm512_store_si512(&B->acc1[n], acc[g]); B->tag1[n] = (uint32_t)(((i - ib) << 2) | g); \
            n += (_mm512_test_epi64_mask(acc[g], acc[g]) != 0); \
        } \
        phi0 = _mm512_add_epi64(phi0, S->inc0); phi0 = _mm512_mask_sub_epi64(phi0, _mm512_cmpge_epu64_mask(phi0, S->mod0), phi0, S->mod0); \
        phi1 = _mm512_add_epi64(phi1, S->inc1); phi1 = _mm512_mask_sub_epi64(phi1, _mm512_cmpge_epu64_mask(phi1, S->mod1), phi1, S->mod1); \
    } \
    return n; \
}
DEFINE_STAGE1(0) DEFINE_STAGE1(1) DEFINE_STAGE1(2) DEFINE_STAGE1(3) DEFINE_STAGE1(4) DEFINE_STAGE1(5) DEFINE_STAGE1(6) DEFINE_STAGE1(7) DEFINE_STAGE1(8)
DEFINE_STAGE1(9) DEFINE_STAGE1(10) DEFINE_STAGE1(11) DEFINE_STAGE1(12) DEFINE_STAGE1(13) DEFINE_STAGE1(14) DEFINE_STAGE1(15)
typedef size_t (*stage1_fn)(const uint8_t **, uint32_t, const vstate_t *, __m512i, __m512i, int64_t, int64_t, bufs_t *);
static stage1_fn stage1_table[16] = { stage1_0, stage1_1, stage1_2, stage1_3, stage1_4, stage1_5, stage1_6, stage1_7, stage1_8, stage1_9, stage1_10, stage1_11, stage1_12, stage1_13, stage1_14, stage1_15 };

/* search denominators c0 .. c0+G-1 */
static void search_group(config_t *C, thread_result_t *R, cand_ctx_t *X, int64_t idx0, const int64_t *cs, bufs_t *B) {
    uint64_t ts0 = __rdtsc();
    vstate_t S;
    uint64_t base[16], inc[16], mod[16];
    for (int j = 0; j < 16; j++) { base[j] = 0; inc[j] = 0; mod[j] = ~0ULL >> 1; }
    for (int j = 0; j < C->nvec; j++) {
        vmod_t *V = &C->vec[j];
        if (!C->odd) {
            uint32_t cm = fastmod_u32((uint32_t)idx0, V->magic, V->M);
            base[j] = (uint64_t)(V->rows + (size_t)cm * V->rowbytes);
        } else {
            uint8_t *dst = X->scratch + (size_t)j * G * V->rowbytes;
            for (int g = 0; g < G; g++) {
                uint32_t cm = (uint32_t)(((u128)cs[g]) % V->M);
                memcpy(dst + (size_t)g * V->rowbytes, V->rows + (size_t)cm * V->rowbytes, V->rowbytes);
            }
            base[j] = (uint64_t)dst;
        }
        inc[j] = 512 % V->M; mod[j] = V->M;
    }
    const uint8_t *prow[G]; uint32_t pmask = C->m2/512 - 1;
    for (int g = 0; g < G; g++) prow[g] = C->W2 + (size_t)((uint32_t)cs[g] & (C->m2 - 1)) * (C->m2/512) * 64;
    S.base0 = _mm512_loadu_si512(base); S.base1 = _mm512_loadu_si512(base+8);
    S.inc0 = _mm512_loadu_si512(inc); S.inc1 = _mm512_loadu_si512(inc+8);
    S.mod0 = _mm512_loadu_si512(mod); S.mod1 = _mm512_loadu_si512(mod+8);
    S.rowbytes = C->rowbytes;
    int K1 = C->K1, K2 = C->K2;
    stage1_fn s1 = stage1_table[K1];
    R->t_setup += __rdtsc() - ts0;

    int64_t prev_hi_g[G]; for (int g = 0; g < G; g++) prev_hi_g[g] = INT64_MIN;
    for (int k = 0; k < C->nint; k++) {
        /* per-c a-ranges (disjoint across intervals for each c), and their union over the group */
        int64_t alo_g[G], ahi_g[G], alo = INT64_MAX, ahi = INT64_MIN;
        for (int g = 0; g < G; g++) {
            int64_t c = cs[g];
            double xlo = C->ilo[k] * (double)c, xhi = C->ihi[k] * (double)c;
            int64_t lo = xlo <= -(double)C->N ? -C->N : (int64_t)floor(xlo) - 1;
            int64_t hi = xhi >= (double)C->N ? C->N : (int64_t)ceil(xhi) + 1;
            if (lo < -C->N) lo = -C->N;
            if (hi > C->N) hi = C->N;
            if (lo <= prev_hi_g[g]) lo = prev_hi_g[g] + 1;
            alo_g[g] = lo; ahi_g[g] = hi;
            if (lo > hi) continue;
            prev_hi_g[g] = hi;
            if (lo < alo) alo = lo;
            if (hi > ahi) ahi = hi;
        }
        if (alo > ahi) continue;
        int64_t i0 = (alo - C->a_base) >> 9, i1 = (ahi - C->a_base) >> 9;
        for (int64_t ib = i0; ib <= i1; ib += BLOCK_CHUNKS) {
            int64_t ie = ib + BLOCK_CHUNKS - 1; if (ie > i1) ie = i1;
            uint64_t ph[16] __attribute__((aligned(64)));
            int64_t a0 = C->a_base + 512 * ib;
            for (int j = 0; j < 16; j++) ph[j] = 0;
            for (int j = 0; j < C->nvec; j++) { int64_t r = a0 % (int64_t)C->vec[j].M; if (r < 0) r += C->vec[j].M; ph[j] = (uint64_t)r; }
            __m512i phi0 = _mm512_load_si512(ph), phi1 = _mm512_load_si512(ph+8);
            uint64_t tt0 = __rdtsc();
            size_t n = s1(prow, pmask, &S, phi0, phi1, ib, ie, B);
            uint64_t tt1 = __rdtsc();
            R->chunks += (uint64_t)(ie - ib + 1) * G;
            R->chunks_s2 += n;
            size_t n2 = 0;
            __m512i *b1 = B->acc1, *b2 = B->acc2; uint32_t *t1 = B->tag1, *t2 = B->tag2;
            if (K2 > K1) {
                for (size_t q = 0; q < n; q++) {
                    __m512i acc = b1[q]; uint32_t tag = t1[q]; const phi_t *P = &B->phitab[tag >> 2];
                    apply_moduli(&acc, P->phi0, P->phi1, &S, (int)(tag & 3), K1, K2);
                    _mm512_store_si512(&b2[n2], acc); t2[n2] = tag;
                    n2 += (_mm512_test_epi64_mask(acc, acc) != 0);
                }
            } else { b2 = b1; t2 = t1; n2 = n; }
            uint64_t tt2 = __rdtsc();
            R->chunks_s3 += n2;
            uint64_t tcand = 0;
            for (size_t q = 0; q < n2; q++) {
                __m512i acc = b2[q];
                uint32_t tag = t2[q]; const phi_t *P = &B->phitab[tag >> 2];
                int g = (int)(tag & 3);
                int64_t i = ib + (int64_t)(tag >> 2);
                int j;
                for (j = K2; j < C->nvec; j++) {
                    apply_moduli(&acc, P->phi0, P->phi1, &S, g, j, j+1);
                    if (_mm512_test_epi64_mask(acc, acc) == 0) break;
                }
                if (j < C->nvec) continue;
                int64_t c = cs[g];
                if (c > C->N || c <= 0) continue;
                uint64_t w[8] __attribute__((aligned(64)));
                _mm512_store_si512(w, acc);
                uint64_t tc0 = __rdtsc();
                for (int j2 = 0; j2 < C->nsca; j2++) X->cmod[j2] = fastmod_u32((uint32_t)c, C->sca[j2]->magic, C->sca[j2]->q);
                for (int l = 0; l < 8; l++) {
                    uint64_t x = w[l];
                    while (x) {
                        int b = __builtin_ctzll(x); x &= x - 1;
                        int64_t a = C->a_base + 512 * i + 64 * l + b;
                        if (a < alo_g[g] || a > ahi_g[g]) continue;
                        R->cand_vec++;
                        check_candidate(C, R, X, a, c);
                    }
                }
                tcand += __rdtsc() - tc0;
            }
            uint64_t tt3 = __rdtsc();
            R->t_s1 += tt1 - tt0; R->t_s2 += tt2 - tt1; R->t_s3 += tt3 - tt2 - tcand; R->t_cand += tcand;
        }
    }
}

/* ------------------------------------------------------------------ */
/* modulus selection: knapsack over prime powers (weight log q, value -log density, capacity log Mmax, at most maxfac
 * items); the powers of one prime form a group (items consecutive in pp) of which at most one may be taken, and
 * none once the prime is in an earlier modulus (gused) */
static double knapsack_pick(const pptab_t *pp, const int *grp, int npp, const int *gused, double logMmax, int *chosen, int *nchosen, int maxfac) {
    /* DP with the capacity discretised to 0.02 in log2 units */
    const double unit = 0.02;
    int cap = (int)(logMmax / unit);
    int ngrp = npp ? grp[npp - 1] + 1 : 0;
    double *best = calloc(cap + 1, sizeof(double)), *bold = malloc((cap + 1) * sizeof(double));
    int *cnt = calloc(cap + 1, sizeof(int)), *cold = malloc((cap + 1) * sizeof(int));
    int *choice = malloc((size_t)(ngrp + 1) * (cap + 1) * sizeof(int));   /* choice[g][w]: item of group g taken at weight w, or -1 */
    for (size_t k = 0; k < (size_t)(ngrp + 1) * (cap + 1); k++) choice[k] = -1;
    for (int i = 0; i < npp; ) {
        int g = grp[i], i0 = i; while (i < npp && grp[i] == g) i++;
        if (gused[g]) continue;
        memcpy(bold, best, (cap + 1) * sizeof(double)); memcpy(cold, cnt, (cap + 1) * sizeof(int));
        int *ch = choice + (size_t)g * (cap + 1);
        for (int k = i0; k < i; k++) {   /* every item of the group extends the state before the group */
            double val = -log2(pp[k].density); int wt = (int)ceil(log2((double)pp[k].q) / unit);
            for (int w = cap; w >= wt; w--) if (cold[w - wt] < maxfac && bold[w - wt] + val > best[w] + 1e-12) { best[w] = bold[w - wt] + val; cnt[w] = cold[w - wt] + 1; ch[w] = k; }
        }
    }
    int wbest = 0; for (int w = 0; w <= cap; w++) if (best[w] > best[wbest]) wbest = w;
    *nchosen = 0;
    for (int g = ngrp - 1, w = wbest; g >= 0; g--) { int k = choice[(size_t)g * (cap + 1) + w]; if (k >= 0) { chosen[(*nchosen)++] = k; w -= (int)ceil(log2((double)pp[k].q) / unit); } }
    double v = best[wbest];
    free(best); free(bold); free(cnt); free(cold); free(choice);
    return v;
}

/* Codex round 34: value-taking options must have an argument, integers must parse completely */
static int opt_int(int argc, char **argv, int *i, int64_t *out) {
    if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", argv[*i]); return 0; }
    char *end = NULL; errno = 0; int64_t v = strtoll(argv[*i + 1], &end, 10);
    if (errno == ERANGE || end == argv[*i + 1] || *end != 0) { fprintf(stderr, "option %s: invalid integer '%s'\n", argv[*i], argv[*i + 1]); return 0; }
    (*i)++; *out = v; return 1;
}
static int opt_double(int argc, char **argv, int *i, double *out) {
    if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", argv[*i]); return 0; }
    char *end = NULL; errno = 0; double v = strtod(argv[*i + 1], &end);
    if (errno == ERANGE || end == argv[*i + 1] || *end != 0) { fprintf(stderr, "option %s: invalid number '%s'\n", argv[*i], argv[*i + 1]); return 0; }
    (*i)++; *out = v; return 1;
}
/* returns 0 on success (point count in *npts_out), 1 on a usage/option error (Codex round 21: status separated from the count) */
static int run_main(int argc, char **argv, mpz_t *fover, int64_t Nover, int quiet_override, uint64_t *npts_out) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s 'f0 f1 ... f6' N [-t threads] [-v] [-K1 k] [-K2 k] [-nvec n] [-pw2 m] [-Mmax m] [-Pmax p] [-q]\n       %s -batch N [-batchpts] [options] < records   (records \"tag... f0 .. f6\": prints \"npts record\", or with -batchpts \"CURVE record\" + the points)\n", argv[0], argv[0]);
        return 1;
    }
    config_t C; memset(&C, 0, sizeof C);
    for (int i = 0; i <= MAXDEG; i++) mpz_init(C.f[i]);   /* missing trailing coefficients stay 0 */
    if (fover) { for (int i = 0; i <= MAXDEG; i++) mpz_set(C.f[i], fover[i]); C.N = Nover; }
    else {
        char *s = strdup(argv[1]); int i = 0; const char *bad = NULL;
        char *tok = strtok(s, " ,");
        while (tok && i <= MAXDEG) { if (!parse_coeff(C.f[i++], tok)) { bad = tok; break; } tok = strtok(NULL, " ,"); }
        if (bad) { fprintf(stderr, "g2rp: invalid coefficient '%s'\n", bad); free(s); clear_coeffs(C.f); return 1; }
        if (tok) { fprintf(stderr, "g2rp: more than %d coefficients given (extra token '%s')\n", MAXDEG + 1, tok); free(s); clear_coeffs(C.f); return 1; }   /* Codex round 35 */
        if (i == 0) { fprintf(stderr, "g2rp: no coefficients given\n"); free(s); clear_coeffs(C.f); return 1; }
        free(s);
        { char *end = NULL; errno = 0; C.N = strtoll(argv[2], &end, 10); if (errno == ERANGE || end == argv[2] || *end != 0 || C.N <= 0) { fprintf(stderr, "N must be a positive integer, not '%s'\n", argv[2]); clear_coeffs(C.f); return 1; } }   /* Codex round 34 */
    }
    int nthreads = omp_get_max_threads();
    int K1 = -1, K2 = -1, nvec = 12; uint32_t pw2 = 0; uint32_t Mmax = 0, Pmax = 0; int quiet = 0; int maxfac = 4;
    double target = 0.02;   /* target candidates per chunk after vector stage */
    int noodd = 0;
    for (int i = 3; i < argc; i++) {
        int64_t iv = 0; double dv = 0;
        #define OPTI(name, var) else if (!strcmp(argv[i], name)) { if (!opt_int(argc, argv, &i, &iv)) { clear_coeffs(C.f); return 1; } if (iv < -1 || iv > 100000000) { fprintf(stderr, "option %s: value %ld out of range\n", name, (long)iv); clear_coeffs(C.f); return 1; } var = (int)iv; }
        if (0) {}
        OPTI("-t", nthreads)
        else if (!strcmp(argv[i], "-v")) verbose++;
        else if (!strcmp(argv[i], "-q")) quiet = 1;
        OPTI("-K1", K1)
        OPTI("-K2", K2)
        OPTI("-nvec", nvec)
        OPTI("-pw2", pw2)
        OPTI("-Mmax", Mmax)
        OPTI("-Pmax", Pmax)
        OPTI("-maxfac", maxfac)
        else if (!strcmp(argv[i], "-target")) { if (!opt_double(argc, argv, &i, &dv)) { clear_coeffs(C.f); return 1; } target = dv; }
        else if (!strcmp(argv[i], "-noodd")) noodd = 1;
        else if (!strcmp(argv[i], "-batch")) { i++; }
        else if (!strcmp(argv[i], "-batchpts")) { }   /* handled in main: batch mode printing "CURVE <record>" + the points of each curve */
        #undef OPTI
        else { fprintf(stderr, "unknown option %s\n", argv[i]); clear_coeffs(C.f); return 1; }
    }
    /* Codex round 35: option domains (-1 / 0 mean automatic where the defaults say so) */
    if (nthreads < 1 || nthreads > 4096) { fprintf(stderr, "-t must be in [1, 4096]\n"); clear_coeffs(C.f); return 1; }
    if (K1 < -1 || K2 < -1) { fprintf(stderr, "-K1/-K2 must be >= -1\n"); clear_coeffs(C.f); return 1; }
    if (nvec < 1 || nvec > MAXVEC - 1) { fprintf(stderr, "-nvec must be in [1, %d]\n", MAXVEC - 1); clear_coeffs(C.f); return 1; }
    if (pw2 > (1u << 24) || (pw2 > 0 && (pw2 & (pw2 - 1)) != 0)) { fprintf(stderr, "-pw2 must be 0 (automatic) or a power of two <= 2^24\n"); clear_coeffs(C.f); return 1; }
    if (Mmax > 100000000u || Pmax > 10000000u) { fprintf(stderr, "-Mmax must be <= 1e8 and -Pmax <= 1e7 (0 = automatic)\n"); clear_coeffs(C.f); return 1; }
    if (maxfac < 1 || maxfac > 12) { fprintf(stderr, "-maxfac must be in [1, 12]\n"); clear_coeffs(C.f); return 1; }
    if (!(target > 0) || !isfinite(target)) { fprintf(stderr, "-target must be a positive finite number\n"); clear_coeffs(C.f); return 1; }
    if (quiet_override) quiet = 1;
    if (nvec > MAXVEC - 1) nvec = MAXVEC - 1;
    /* automatic setup size: at small heights the search is tiny and the setup (tables ~ sum q^2) must be light;
     * the sieving power still has to bring the candidate count down to a few thousand (N=1e4: ~1e-5) */
    int light = C.N <= SMALL_N;
    if (Pmax == 0) Pmax = C.N <= 3000 ? 100 : C.N <= 10000 ? 150 : light ? 200 : 1100;
    if (Mmax == 0) Mmax = C.N <= 3000 ? 600 : C.N <= 10000 ? 1000 : light ? 1500 : C.N <= 200000 ? 3000 : 12000;
    uint32_t pw2max = pw2 ? pw2 : light ? 512 : 8192;
    omp_set_num_threads(nthreads);
    double t0 = omp_get_wtime();

    /* real intervals where f >= 0 (from the scaled coefficients: exact to 64 bits whatever the coefficient sizes) */
    {
        long double cf[MAXDEG+1]; scaled_coeffs(C.f, cf);
        /* Codex round 33: the scaled coefficients keep only the top 64 bits; with larger coefficients cancellation can be lost
         * and the pruning would be unsound (points in a region where the truncated f is negative but f > 0), so prune only
         * when every coefficient is exactly representable; otherwise search the full range (sieve + exact test stay correct). */
        int exact_cf = 1; for (int i = 0; i <= MAXDEG; i++) if (mpz_sizeinbase(C.f[i], 2) > 63) exact_cf = 0;
        long double roots[MAXDEG+1]; int nr = exact_cf ? real_roots(cf, MAXDEG, roots) : -1;
        if (nr < 0) { if (verbose) fprintf(stderr, "real roots: overflow (coefficient sizes too spread), searching without the real-root restriction\n"); nr = 0; }
        for (int i = 0; i < nr; i++) for (int j = i+1; j < nr; j++) if (roots[j] < roots[i]) { long double t = roots[i]; roots[i] = roots[j]; roots[j] = t; }
        long double pts[MAXDEG+3]; int np = 0; pts[np++] = -INFINITY;
        for (int i = 0; i < nr; i++) pts[np++] = roots[i];
        pts[np++] = INFINITY;
        C.nint = 0;
        for (int i = 0; i + 1 < np; i++) {
            long double lo = pts[i], hi = pts[i+1];
            long double mid = isinf(lo) ? (isinf(hi) ? 0 : hi - 1) : (isinf(hi) ? lo + 1 : 0.5L*(lo+hi));
            if (horner(cf, MAXDEG, mid) >= 0 || (isinf(lo) && isinf(hi))) {
                long double eps = 1e-9L;
                double l = (double)(isinf(lo) ? lo : lo - eps*(fabsl(lo)+1)), h = (double)(isinf(hi) ? hi : hi + eps*(fabsl(hi)+1));
                if (C.nint > 0 && C.ihi[C.nint-1] >= l) C.ihi[C.nint-1] = h;
                else { C.ilo[C.nint] = l; C.ihi[C.nint] = h; C.nint++; }
            }
        }
        if (verbose) {
            fprintf(stderr, "real roots (%d):", nr); for (int i = 0; i < nr; i++) fprintf(stderr, " %.6f", (double)roots[i]);
            fprintf(stderr, "\nintervals (%d):", C.nint); for (int i = 0; i < C.nint; i++) fprintf(stderr, " [%.6g, %.6g]", C.ilo[i], C.ihi[i]);
            fprintf(stderr, "\n");
        }
    }
    /* C.a_base (origin of the chunk index) is set after the power of 2 is chosen: it must be a multiple of m2 */

    double ts_a = omp_get_wtime();
    /* ---- prime powers: exact O(q) density of every q = p^k (p odd prime <= Pmax; q = p, or q <= Mmax); the powers
     *      of a prime form a group for the knapsack; only the tables actually used are built (below) ---- */
    int npp = 0; pptab_t *pp = calloc(1024, sizeof(pptab_t)); int *grp = malloc(1024 * sizeof(int)); int ngrp = 0;
    for (uint32_t p = 3; p <= Pmax; p += 2) {
        if (!is_prime(p)) continue;
        int any = 0;
        for (uint64_t q = p; q == p || q <= Mmax; q *= p) {
            pprows_t R; pp_rows(C.f, (uint32_t)q, &R); double d = pp_density(&R, p); pp_rows_free(&R);
            if (verbose > 1) fprintf(stderr, "  density q=%lu : %.4f\n", q, d);
            if (d < 0.999 && npp < 1024) { pp[npp].q = (uint32_t)q; pp[npp].p = p; pp[npp].density = d; grp[npp] = ngrp; npp++; any = 1; }
        }
        ngrp += any;
    }
    double ts_b = omp_get_wtime();
    /* power of 2 */
    {
        double bestd = 2; uint32_t bestm = 512;
        for (uint32_t m = 512; m <= pw2max; m *= 2) {
            if (pw2 && m != pw2) continue;
            pprows_t R; pp_rows(C.f, m, &R); double d = pp_density(&R, 2); pp_rows_free(&R);
            if (verbose > 1) fprintf(stderr, "  density m=%u : %.4f\n", m, d);
            if (d < bestd - 0.01) { bestd = d; bestm = m; }
        }
        pptab_t T; build_pptab(&T, bestm, 2, C.f);
        C.m2 = bestm; C.d2 = T.density;
        /* the pw2 rows are read at word (chunk index mod m2/512): the chunk origin must be a multiple of m2
         * (the old 512-aligned origin missed points whenever m2 > 512 did not divide ceil(N/512)*512) */
        C.a_base = -(((C.N + bestm - 1) / bestm) * (int64_t)bestm);
        uint32_t nw = bestm / 512;
        C.W2 = aligned_alloc(64, (size_t)64 * bestm * nw);
        C.rowempty2 = calloc(bestm, 1);
        for (uint32_t cb = 0; cb < bestm; cb++) {
            uint64_t *w = (uint64_t*)(C.W2 + (size_t)cb * nw * 64);
            pack_bits(T.ok + (size_t)cb * bestm, bestm, w);
            uint64_t any = 0; for (uint32_t k = 0; k < nw * 8; k++) any |= w[k];
            C.rowempty2[cb] = (any == 0);
        }
        free(T.ok);
    }
    double ts_c = omp_get_wtime();
    /* ---- choose vector moduli by repeated grouped knapsacks ---- */
    int *gused = calloc(ngrp + 1, sizeof(int));
    C.nvec = 0;
    double prod = C.d2;
    uint32_t selM[MAXVEC]; int selnf[MAXVEC]; pptab_t *selfac[MAXVEC][8];
    while (C.nvec < nvec && prod * 512 > target) {
        int chosen[16], nch;
        double v = knapsack_pick(pp, grp, npp, gused, log2((double)Mmax), chosen, &nch, maxfac);
        if (nch == 0 || v <= 0) break;
        uint32_t M = 1; double d = 1;
        for (int k = 0; k < nch; k++) { M *= pp[chosen[k]].q; selfac[C.nvec][k] = &pp[chosen[k]]; gused[grp[chosen[k]]] = 1; d *= pp[chosen[k]].density; }
        selM[C.nvec] = M; selnf[C.nvec] = nch;
        prod *= d;
        C.nvec++;
    }
    /* scalar tables: for each remaining prime its densest power <= 512; the 24 densest of those */
    C.sca = calloc(npp + 1, sizeof(pptab_t*)); C.nsca = 0;
    for (int i = 0; i < npp; i++) {
        if (gused[grp[i]] || pp[i].q > 512) continue;
        if (C.nsca && C.sca[C.nsca-1]->p == pp[i].p) { if (pp[i].density < C.sca[C.nsca-1]->density) C.sca[C.nsca-1] = &pp[i]; }
        else C.sca[C.nsca++] = &pp[i];
    }
    for (int i = 0; i < C.nsca; i++) for (int j = i+1; j < C.nsca; j++) if (C.sca[j]->density < C.sca[i]->density) { pptab_t *t = C.sca[i]; C.sca[i] = C.sca[j]; C.sca[j] = t; }
    if (C.nsca > 24) C.nsca = 24;
    free(gused); free(grp);
    /* build the tables in use in one parallel region, large first (those above 4M entries internally parallel) */
    {
        pptab_t *need[MAXVEC * 8 + 24]; int nneed = 0;
        for (int j = 0; j < C.nvec; j++) for (int k = 0; k < selnf[j]; k++) need[nneed++] = selfac[j][k];
        for (int i = 0; i < C.nsca; i++) need[nneed++] = C.sca[i];
        for (int i = 0; i < nneed; i++) for (int j = i + 1; j < nneed; j++) if (need[j]->q > need[i]->q) { pptab_t *t = need[i]; need[i] = need[j]; need[j] = t; }
        int k0 = 0;
        for (; k0 < nneed && (uint64_t)need[k0]->q * need[k0]->q > (1u << 22); k0++) build_pptab(need[k0], need[k0]->q, need[k0]->p, C.f);
        size_t rest = 0; for (int k = k0; k < nneed; k++) rest += (size_t)need[k]->q * need[k]->q;
        #pragma omp parallel for schedule(dynamic, 1) if(rest > (1u << 22))
        for (int k = k0; k < nneed; k++) build_pptab(need[k], need[k]->q, need[k]->p, C.f);
    }
    double ts_c2 = omp_get_wtime();
    size_t rowbytes = 0;
    for (int j = 0; j < C.nvec; j++) { size_t w = (selM[j] + 512 + 64 + 63) / 64 + 1; if (w * 8 > rowbytes) rowbytes = w * 8; }
    for (int j = 0; j < C.nvec; j++) build_vmod(&C.vec[j], selM[j], selnf[j], selfac[j], rowbytes);
    C.rowbytes = rowbytes;
    C.odd = (mpz_sgn(C.f[6]) == 0 && mpz_sgn(C.f[5]) != 0 && !noodd);
    if (C.odd) {
        /* squarefree divisors e <= N of rad(|f5|): only the primes p <= N dividing f5 can occur in such an e, so f5 is
         * trial-divided by p <= min(N, sqrt(cofactor)) -- with GMP while the cofactor exceeds 64 bits, then in uint64 */
        int64_t primes[64]; int np = 0;
        mpz_t L; mpz_init(L); mpz_abs(L, C.f[5]);
        int64_t p = 2;
        for (; p <= C.N && !mpz_fits_ulong_p(L); p = p == 2 ? 3 : p + 2)
            if (mpz_divisible_ui_p(L, (unsigned long)p)) { if (np < 64) primes[np++] = p; do mpz_divexact_ui(L, L, (unsigned long)p); while (mpz_divisible_ui_p(L, (unsigned long)p)); }
        if (mpz_fits_ulong_p(L)) {
            uint64_t Lu = mpz_get_ui(L), q = (uint64_t)p;
            for (; q <= (uint64_t)C.N && Lu / q >= q; q = q == 2 ? 3 : q + 2) if (Lu % q == 0) { if (np < 64) primes[np++] = (int64_t)q; do Lu /= q; while (Lu % q == 0); }
            if (Lu > 1 && Lu / q < q && Lu <= (uint64_t)C.N && np < 64) primes[np++] = (int64_t)Lu;   /* the cofactor is prime */
        }
        mpz_clear(L);
        if (np > 12) { fprintf(stderr, "leading coefficient has too many prime factors for the quintic path; use -noodd\n"); exit(1); }
        C.ne = 0;
        for (int mask = 0; mask < (1 << np); mask++) {
            int64_t e = 1; int ok = 1;
            for (int i = 0; i < np && ok; i++) if (mask & (1 << i)) { if (e > C.N / primes[i]) ok = 0; else e *= primes[i]; }   /* e <= N without overflow */
            if (ok) C.elist[C.ne++] = e;
        }
        if (verbose) { fprintf(stderr, "quintic model: denominators c = e*d^2, e in {"); for (int i = 0; i < C.ne; i++) fprintf(stderr, "%s%ld", i ? "," : "", C.elist[i]); fprintf(stderr, "}\n"); }
    }
    double ts_d = omp_get_wtime();
    /* scalar tables: bit rows and the phase of a_base */
    for (int i = 0; i < C.nsca; i++) { int64_t r = C.a_base % (int64_t)C.sca[i]->q; if (r < 0) r += C.sca[i]->q; C.sca[i]->abase_mod = (uint32_t)r; build_bits(C.sca[i]); }
    /* stage boundaries */
    {
        double pr = C.d2; int j;
        if (K1 < 0) { for (j = 0; j < C.nvec && pr * 512 > 0.7; j++) pr *= C.vec[j].density; K1 = j; }
        else { for (j = 0; j < K1 && j < C.nvec; j++) pr *= C.vec[j].density; }
        if (K2 < 0) { for (j = K1; j < C.nvec && pr * 512 > 0.03; j++) pr *= C.vec[j].density; K2 = j; }
    }
    if (K1 > C.nvec) K1 = C.nvec;
    if (K2 > C.nvec) K2 = C.nvec;
    if (K2 < K1) K2 = K1;
    C.K1 = K1; C.K2 = K2;
    if (verbose) {
        fprintf(stderr, "pw2 %u density %.4f\n", C.m2, C.d2);
        fprintf(stderr, "vector moduli (%d):", C.nvec); double pr = C.d2;
        for (int j = 0; j < C.nvec; j++) {
            pr *= C.vec[j].density; fprintf(stderr, " %u=", C.vec[j].M);
            for (int k = 0; k < C.vec[j].nf; k++) fprintf(stderr, "%s%u", k ? "*" : "", C.vec[j].fac[k]->q);
            fprintf(stderr, "(%.3f)", C.vec[j].density);
            if (j+1 == K1) fprintf(stderr, " |"); if (j+1 == K2) fprintf(stderr, " ||");
        }
        fprintf(stderr, "\n  cumulative density after vector stage: %.3g (%.4f cand/chunk); K1=%d K2=%d\n", pr, pr*512, K1, K2);
        fprintf(stderr, "scalar moduli (%d):", C.nsca); for (int j = 0; j < C.nsca && j < 40; j++) { fprintf(stderr, " %u", C.sca[j]->q); }
        for (int j = 0; j < C.nsca; j++) pr *= C.sca[j]->density;
        fprintf(stderr, "%s\n  cumulative density after scalar stage: %.3g\n", C.nsca > 40 ? " ..." : "", pr);
        size_t tb = 0; for (int j = 0; j < C.nvec; j++) tb += (size_t)(C.vec[j].M + G) * C.vec[j].rowbytes;
        fprintf(stderr, "vector tables %.1f MB; setup time %.2fs (pp scan %.2f, pw2 %.2f, tables %.2f, vmod build %.2f, rest %.2f)\n", tb / 1048576.0, omp_get_wtime() - t0, ts_b - ts_a, ts_c - ts_b, ts_c2 - ts_c, ts_d - ts_c2, omp_get_wtime() - ts_d);
    }
    omp_set_num_threads(nthreads);
    double t1 = omp_get_wtime();

    thread_result_t *res = calloc(nthreads, sizeof(thread_result_t));
    /* a tiny search (below ~1e5 chunks, i.e. a few ms of work) runs on the calling thread: forking a team costs more,
     * especially on a loaded machine */
    int par = nthreads > 1 && (double)C.N * ((2.0 * C.N) / 512 + 1) > 1e5;
    #pragma omp parallel if(par)
    {
        int tid = omp_get_thread_num();
        thread_result_t *R = &res[tid];
        bufs_t B;
        B.acc1 = aligned_alloc(64, 64 * (BLOCK_CHUNKS * G + 1)); B.acc2 = aligned_alloc(64, 64 * (BLOCK_CHUNKS * G + 1));
        B.tag1 = malloc(sizeof(uint32_t) * (BLOCK_CHUNKS * G + 1)); B.tag2 = malloc(sizeof(uint32_t) * (BLOCK_CHUNKS * G + 1));
        B.phitab = aligned_alloc(64, sizeof(phi_t) * BLOCK_CHUNKS);
        cand_ctx_t X; X.cmod = malloc(sizeof(uint32_t) * (C.nsca + 1)); mpz_init(X.F); mpz_init(X.t1); mpz_init(X.t2);
        X.scratch = aligned_alloc(64, (size_t)(C.nvec + 1) * G * C.rowbytes + 64);
        if (!C.odd) {
            #pragma omp for schedule(dynamic, 4)
            for (int64_t c0 = 1; c0 <= C.N; c0 += G) {
                int64_t cs[G]; for (int g = 0; g < G; g++) cs[g] = c0 + g;
                search_group(&C, R, &X, c0, cs, &B);
            }
        } else {
            /* quintic model: c = e * d^2, e squarefree with all primes dividing f5; rows are d-indexed (built per e) */
            for (int ie = 0; ie < C.ne; ie++) {
                int64_t e = C.elist[ie];
                int64_t dmax = (int64_t)floor(sqrt((double)C.N / (double)e) + 1); while (e * dmax * dmax > C.N) dmax--;
                #pragma omp for schedule(dynamic, 1) nowait
                for (int64_t d0 = 1; d0 <= dmax; d0 += G) {
                    int64_t cs[G]; for (int g = 0; g < G; g++) { int64_t d = d0 + g; cs[g] = e * d * d; }
                    search_group(&C, R, &X, d0, cs, &B);
                }
            }
        }
        mpz_clear(X.F); mpz_clear(X.t1); mpz_clear(X.t2); free(X.cmod); free(X.scratch);
        free(B.acc1); free(B.acc2); free(B.tag1); free(B.tag2); free(B.phitab);
    }
    double t2 = omp_get_wtime();

    uint64_t npts = 0, nx = 0;
    {
        mpz_t f6; mpz_init_set(f6, C.f[MAXDEG]);
        if (mpz_sgn(f6) == 0) { npts += 1; if (!quiet) printf("(1 : 0 : 0)\n"); }
        else if (mpz_sgn(f6) > 0 && mpz_perfect_square_p(f6)) { mpz_sqrt(f6, f6); npts += 2; if (!quiet) gmp_printf("(1 : %Zd : 0)\n(1 : -%Zd : 0)\n", f6, f6); }
        mpz_clear(f6);
    }
    thread_result_t tot; memset(&tot, 0, sizeof tot);
    for (int t = 0; t < nthreads; t++) {
        tot.cand_vec += res[t].cand_vec; tot.cand_scalar += res[t].cand_scalar; tot.cand_gmp += res[t].cand_gmp;
        tot.chunks += res[t].chunks; tot.chunks_s2 += res[t].chunks_s2; tot.chunks_s3 += res[t].chunks_s3; tot.c_skipped += res[t].c_skipped;
        tot.t_s1 += res[t].t_s1; tot.t_s2 += res[t].t_s2; tot.t_s3 += res[t].t_s3; tot.t_cand += res[t].t_cand; tot.t_setup += res[t].t_setup;
        for (size_t i = 0; i < res[t].n; i++) {
            point_t *p = &res[t].pts[i]; nx++;
            if (mpz_sgn(p->b) == 0) { npts += 1; if (!quiet) printf("(%ld : 0 : %ld)\n", p->a, p->c); }
            else { npts += 2; if (!quiet) gmp_printf("(%ld : %Zd : %ld)\n(%ld : -%Zd : %ld)\n", p->a, p->b, p->c, p->a, p->b, p->c); }
        }
    }
    fprintf(stderr, "N=%ld: %lu points (%lu finite x-values), search %.3fs (setup %.2fs), %d threads\n", C.N, npts, nx, t2 - t1, t1 - t0, nthreads);
    if (verbose) {
        fprintf(stderr, "chunks %lu (s2 %lu = %.3f, s3 %lu = %.4f), c skipped %lu, cand vec %lu scalar %lu gmp %lu\n",
                tot.chunks, tot.chunks_s2, (double)tot.chunks_s2/tot.chunks, tot.chunks_s3, (double)tot.chunks_s3/tot.chunks, tot.c_skipped, tot.cand_vec, tot.cand_scalar, tot.cand_gmp);
        double ch = (double)tot.chunks;
        fprintf(stderr, "tsc/chunk: setup %.2f s1 %.2f, s2 %.2f, s3 %.2f, cand %.2f, total %.2f  (%.2f ns/chunk/thread)\n", tot.t_setup/ch, tot.t_s1/ch, tot.t_s2/ch, tot.t_s3/ch, tot.t_cand/ch, (tot.t_setup+tot.t_s1+tot.t_s2+tot.t_s3+tot.t_cand)/ch, (t2-t1)*1e9/ch*nthreads);
    }
    /* free */
    for (int t = 0; t < nthreads; t++) { for (size_t i = 0; i < res[t].n; i++) mpz_clear(res[t].pts[i].b); free(res[t].pts); }
    free(res);
    for (int j = 0; j < C.nvec; j++) { free(C.vec[j].rows); free(C.vec[j].rowempty); }
    free(C.W2); free(C.rowempty2);
    for (int i = 0; i < npp; i++) { free(pp[i].ok); free(pp[i].bits); }
    free(pp); free(C.sca); clear_coeffs(C.f);
    *npts_out = npts; return 0;
}

int main(int argc, char **argv) {
    int batch = 0, batchpts = 0; int64_t Nb = 0;
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "-batchpts")) batchpts = 1;
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "-batch")) { if (i + 1 >= argc) { fprintf(stderr, "-batch needs a height bound N\n"); return 1; } batch = 1; char *end = NULL; errno = 0; Nb = strtoll(argv[i+1], &end, 10); if (errno == ERANGE || end == argv[i+1] || *end != 0 || Nb <= 0) { fprintf(stderr, "-batch: N must be a positive integer, not '%s'\n", argv[i+1]); return 1; } }   /* Codex rounds 22-23: full validation */
    if (!batch) { uint64_t np = 0; return run_main(argc, argv, NULL, 0, 0, &np); }
    /* batch: read "tag... f0 .. f6" lines: the LAST 7 fields are the coefficients; print "npts <line>" */
    char line[8192];
    mpz_t f[7]; for (int i = 0; i < 7; i++) mpz_init(f[i]);
    while (fgets(line, sizeof line, stdin)) {
        char tmp[8192]; strcpy(tmp, line);
        size_t LL = strlen(line); if (LL > 0 && line[LL-1] != '\n' && !feof(stdin)) { fprintf(stderr, "g2rp: batch line too long (> %zu bytes)\n", sizeof line - 2); return 1; }
        char *toks[256]; int nt = 0; char *tok = strtok(tmp, " \t\n,");
        while (tok && nt < 256) { toks[nt++] = tok; tok = strtok(NULL, " \t\n,"); }
        if (tok) { fprintf(stderr, "g2rp: batch record with more than 256 fields: %s", line); return 1; }   /* Codex round 25: never truncate a record */
        if (nt == 0) continue;                                                     /* blank line */
        if (nt < 7) { fprintf(stderr, "g2rp: batch record with fewer than 7 coefficients: %s", line); return 1; }   /* Codex round 24 */
        for (int i = 0; i < 7; i++)   /* whole token, any size (GMP) */
            if (!parse_coeff(f[i], toks[nt - 7 + i])) { fprintf(stderr, "g2rp: invalid coefficient '%s' in batch record: %s", toks[nt - 7 + i], line); return 1; }
        size_t L = strlen(line); if (L && line[L-1] == '\n') line[L-1] = 0;
        if (batchpts) {   /* -batchpts: the record as a "CURVE <record>" header followed by its points (g2rank's stream format) */
            printf("CURVE %s\n", line);
            uint64_t np = 0; if (run_main(argc, argv, f, Nb, 0, &np) != 0) { fprintf(stderr, "g2rp: batch failed on line: %s\n", line); return 1; }
            fflush(stdout); continue;
        }
        uint64_t np = 0; if (run_main(argc, argv, f, Nb, 1, &np) != 0) { fprintf(stderr, "g2rp: batch failed on line: %s\n", line); return 1; }
        printf("%lu %s\n", np, line); fflush(stdout);
    }
    for (int i = 0; i < 7; i++) mpz_clear(f[i]);
    return 0;
}

/* k3pts.c -- rational points of the K3 double plane w^2 = S(x0,x1,x2) (S a plane sextic) up to a height bound.
 *
 * Enumerates primitive (x0, x1, x2) with |xi| <= H (one representative per projective point: x2 > 0, or x2 = 0 and
 * x1 > 0, or (1:0:0)) and prints those with S(x0,x1,x2) a perfect square (or zero: a point on the branch curve).
 * Points on a given list of lines (tritangents / node lines, whose lifts are the "forced" points on every fibre) are
 * skipped.
 *
 * Method: for each (x1, x2) the values S(x0, x1, x2), x0 = x0lo .. H, are sieved with quadratic-residue bit rows
 * modulo a set of moduli chosen automatically: for every prime power q <= Qmax the exact fraction of residue
 * triples with S a square mod q is measured on the q^3 grid, and moduli (products of prime powers, <= Mcap) are
 * packed by a knapsack on -log2(density) until the survivors are cheaper to test exactly than to sieve further.
 * Rows are periodic (period m, stored to cover m + 576 bits) so that a 512-bit window at any phase is two unaligned
 * loads and a VBMI2 funnel shift; four consecutive x1 rows share the per-chunk phases (which depend only on x0).
 * If S is even in x0 only x0 >= 0 is sieved and both signs are reported.  Rows with gcd(x1, x2) divisible by
 * 2, 3, 5 or 7 are ANDed with a period-210 coprimality mask.  Survivors are tested exactly: int128 while
 * sum |S_n| H^6 < 2^125, otherwise through a mod-p quadratic-residue prefilter (four 16-bit primes) and GMP.
 *
 *   k3pts sextic.txt H [-t threads] [-skip lines.txt] [-q] [-v] [-Q qmax] [-M mcap] [-nmod n] [-range x2lo x2hi] [-nosym]
 * output lines: "x0 x1 x2" (a trailing "W" marks S = 0).  Sextic file: 28 coefficients in linesearch order
 * (i = 6..0, j = 6-i..0, k = 6-i-j for x0^i x1^j x2^k).  -range restricts x2 to [x2lo, x2hi] (splitting a run).
 */
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "simd.h"
#include <gmp.h>
#include <omp.h>

typedef __int128 i128;
typedef unsigned __int128 u128;
#define G 4                      /* x1 rows processed together */
#define MAXMOD 64
#define MAXFAC 3
#define NPRE 4                   /* mod-p prefilter primes (GMP regime) */

static int64_t S[28];            /* coefficients */
static int SI[28], SJ[28], SK[28];
static int nskip = 0; static int64_t skipL[256][3];
static int verbose = 0;

/* ------------------------------------------------------------------ */
static inline uint64_t fastmod_magic(uint32_t m) { return UINT64_C(0xFFFFFFFFFFFFFFFF) / m + 1; }
static inline uint32_t fastmod_u32(uint32_t a, uint64_t M, uint32_t m) { uint64_t low = M * a; return (uint32_t)(((u128)low * m) >> 64); }
static inline uint32_t modp(int64_t a, uint32_t p) { int64_t r = a % (int64_t)p; return (uint32_t)(r < 0 ? r + p : r); }
static inline int64_t gcd64(int64_t sa, int64_t sb) {   /* binary gcd */
    uint64_t a = sa < 0 ? -(uint64_t)sa : (uint64_t)sa, b = sb < 0 ? -(uint64_t)sb : (uint64_t)sb;
    if (a == 0) return (int64_t)b;
    if (b == 0) return (int64_t)a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do { b >>= __builtin_ctzll(b); if (a > b) { uint64_t t = a; a = b; b = t; } b -= a; } while (b);
    return (int64_t)(a << sh);
}

/* pack a byte row (values 0/1) into bits, then extend it periodically (period q) over words*64 bits */
static inline void pack_bits(const uint8_t *ok, uint32_t n, uint64_t *w) {
    uint32_t j = 0;
    for (; j + 64 <= n; j += 64) w[j >> 6] = _mm512_cmpneq_epi8_mask(_mm512_loadu_si512(ok + j), _mm512_setzero_si512());
    if (j < n) { __mmask64 k = (__mmask64)((UINT64_C(1) << (n - j)) - 1); w[j >> 6] = _mm512_mask_cmpneq_epi8_mask(k, _mm512_maskz_loadu_epi8(k, ok + j), _mm512_setzero_si512()); }
}
static void extend_periodic(uint64_t *w, uint32_t q, size_t words) {
    size_t nb = words * 64; uint32_t Q = q;
    if (Q < 64) {
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
static inline uint32_t row_words(uint32_t m) { return (m + 63) / 64 + 10; }   /* phase < m, then 8 words + 1 for the funnel shift */

/* ------------------------------------------------------------------ */
/* prime-power tiles: rows (x2q, x1q) of W words, bit x0 = [S(x0,x1q,x2q) is a square mod q, and not p | x0,x1,x2] */
typedef struct { uint32_t q, p; uint32_t W; uint64_t *rows; double density; } tile_t;

static void build_tile(tile_t *T, uint32_t q, uint32_t p, uint32_t W) {
    T->q = q; T->p = p; T->W = W;
    T->rows = aligned_alloc(64, (size_t)q * q * W * 8);
    uint8_t *qr = calloc(q, 1); for (uint64_t v = 0; v < q; v++) qr[(v * v) % q] = 1;
    uint32_t sm[28]; for (int n = 0; n < 28; n++) sm[n] = modp(S[n], q);
    uint64_t magic = fastmod_magic(q);
    uint64_t cnt = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:cnt)
    for (uint32_t x2 = 0; x2 < q; x2++) {
        uint32_t pw2[7]; pw2[0] = 1; for (int e = 1; e <= 6; e++) pw2[e] = fastmod_u32(pw2[e-1] * x2, magic, q);
        uint8_t *ok = malloc(q);
        for (uint32_t x1 = 0; x1 < q; x1++) {
            uint32_t pw1[7]; pw1[0] = 1; for (int e = 1; e <= 6; e++) pw1[e] = fastmod_u32(pw1[e-1] * x1, magic, q);
            uint32_t c[7] = {0,0,0,0,0,0,0};
            for (int n = 0; n < 28; n++) { uint32_t t = fastmod_u32(sm[n] * pw1[SJ[n]], magic, q); t = fastmod_u32(t * pw2[SK[n]], magic, q); c[SI[n]] += t; }
            for (int i = 0; i <= 6; i++) c[i] = fastmod_u32(c[i], magic, q);
            int z12 = (x1 % p == 0 && x2 % p == 0);
            for (uint32_t x0 = 0; x0 < q; x0++) {
                uint32_t v = c[6];
                for (int i = 5; i >= 0; i--) v = fastmod_u32(v * x0 + c[i], magic, q);
                int excl = z12 && x0 % p == 0;   /* p | x0, x1, x2: never primitive */
                ok[x0] = qr[v] && !excl;
                cnt += ok[x0];
            }
            uint64_t *w = T->rows + ((size_t)x2 * q + x1) * W;
            memset(w, 0, (size_t)W * 8); pack_bits(ok, q, w); extend_periodic(w, q, W);
        }
        free(ok);
    }
    free(qr);
    /* density over the residue triples that primitive points can have (not all divisible by p) */
    double tot = (double)q * q * q - (double)(q / p) * (q / p) * (q / p);
    T->density = (double)cnt / tot;
}

/* ------------------------------------------------------------------ */
typedef struct { uint32_t m; int nf; tile_t *fac[MAXFAC]; uint32_t fq[MAXFAC]; double density; uint32_t W; uint64_t *T; } mod_t;
static mod_t MODS[MAXMOD]; static int nmod = 0;

/* rows (x2r, x1r') for x1r' in [0, m + G - 1) (the last G-1 rows repeat rows 0..G-2, so G consecutive x1 are contiguous) */
static void build_modulus(mod_t *V) {
    uint32_t m = V->m, W = row_words(m); V->W = W;
    size_t nrow = (size_t)m * (m + G - 1);
    V->T = aligned_alloc(64, nrow * W * 8);
    if (!V->T) { fprintf(stderr, "allocation failed for modulus %u\n", m); exit(1); }
    #pragma omp parallel for schedule(static)
    for (uint32_t x2r = 0; x2r < m; x2r++) {
        for (uint32_t x1i = 0; x1i < m + G - 1; x1i++) {
            uint32_t x1r = x1i % m;
            uint64_t *w = V->T + ((size_t)x2r * (m + G - 1) + x1i) * W;
            for (int f = 0; f < V->nf; f++) {
                tile_t *F = V->fac[f]; uint32_t q = F->q;
                const uint64_t *src = F->rows + ((size_t)(x2r % q) * q + (x1r % q)) * F->W;
                if (f == 0) memcpy(w, src, (size_t)W * 8); else for (uint32_t k = 0; k < W; k++) w[k] &= src[k];
            }
        }
    }
}

/* knapsack over the candidate tiles (weight log2 q, value -log2 density, capacity log2 Mcap, at most MAXFAC items) */
static double knapsack_pick(tile_t **items, int *used, int nitems, double logMcap, int *chosen, int *nchosen) {
    const double unit = 0.02;
    int cap = (int)(logMcap / unit);
    double *best = calloc(cap + 1, sizeof(double)); int *cnt = calloc(cap + 1, sizeof(int));
    char *from = calloc((size_t)(cap + 1) * nitems, 1);
    for (int i = 0; i < nitems; i++) {
        if (used[i]) continue;
        double val = -log2(items[i]->density); int wt = (int)ceil(log2((double)items[i]->q) / unit);
        char *fr = from + (size_t)i * (cap + 1);
        for (int w = cap; w >= wt; w--) if (cnt[w - wt] < MAXFAC && best[w - wt] + val > best[w] + 1e-12) { best[w] = best[w - wt] + val; cnt[w] = cnt[w - wt] + 1; fr[w] = 1; }
    }
    int wbest = 0; for (int w = 0; w <= cap; w++) if (best[w] > best[wbest]) wbest = w;
    *nchosen = 0; int w = wbest;
    for (int i = nitems - 1; i >= 0; i--) if (from[(size_t)i * (cap + 1) + w]) { chosen[(*nchosen)++] = i; w -= (int)ceil(log2((double)items[i]->q) / unit); }
    double v = best[wbest]; free(best); free(cnt); free(from); return v;
}

static int is_prime(uint32_t n) { if (n < 2) return 0; for (uint32_t d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }

/* choose the moduli: per prime the best power <= Qmax (exact densities), then repeated knapsacks */
static void choose_moduli(uint32_t Qmax, uint32_t Mcap, int nmodmax, double stop_density) {
    tile_t *cand[512]; int ncand = 0;
    uint32_t Wt = row_words(Mcap > Qmax ? Mcap : Qmax);   /* tiles are as wide as the widest possible modulus row */
    for (uint32_t p = 2; p <= Qmax; p++) {
        if (!is_prime(p)) continue;
        tile_t *bestT = NULL;
        for (uint64_t q = p; q <= Qmax; q *= p) {
            tile_t *T = calloc(1, sizeof(tile_t));
            build_tile(T, (uint32_t)q, p, Wt);
            if (verbose > 1) fprintf(stderr, "  q=%lu density %.4f\n", q, T->density);
            /* a higher power must gain at least 0.02 in density over the best lower one */
            if (!bestT || T->density < bestT->density - 0.02) { if (bestT) { free(bestT->rows); free(bestT); } bestT = T; }
            else { free(T->rows); free(T); }
        }
        if (bestT->density < 0.995) cand[ncand++] = bestT; else { free(bestT->rows); free(bestT); }
    }
    int used[512] = {0};
    double prod = 1;
    while (nmod < nmodmax && prod > stop_density) {
        int chosen[MAXFAC + 1], nch;
        double v = knapsack_pick(cand, used, ncand, log2((double)Mcap), chosen, &nch);
        if (nch == 0 || v <= 0) break;
        mod_t *V = &MODS[nmod]; V->m = 1; V->nf = nch; V->density = 1;
        for (int k = 0; k < nch; k++) { V->m *= cand[chosen[k]]->q; V->fac[k] = cand[chosen[k]]; V->fq[k] = cand[chosen[k]]->q; V->density *= cand[chosen[k]]->density; used[chosen[k]] = 1; }
        prod *= V->density; nmod++;
    }
    /* best moduli first (the kernel applies them in order; only the AND order changes) */
    for (int i = 0; i < nmod; i++) for (int j = i + 1; j < nmod; j++) if (MODS[j].density < MODS[i].density) { mod_t t = MODS[i]; MODS[i] = MODS[j]; MODS[j] = t; }
    for (int i = 0; i < nmod; i++) build_modulus(&MODS[i]);
    for (int i = 0; i < ncand; i++) { free(cand[i]->rows); free(cand[i]); }
    if (verbose) {
        fprintf(stderr, "moduli (%d):", nmod); double pr = 1;
        for (int i = 0; i < nmod; i++) { pr *= MODS[i].density; fprintf(stderr, " %u=", MODS[i].m); for (int k = 0; k < MODS[i].nf; k++) fprintf(stderr, "%s%u", k ? "*" : "", MODS[i].fq[k]); fprintf(stderr, "(%.3f)", MODS[i].density); }
        size_t tb = 0; for (int i = 0; i < nmod; i++) tb += (size_t)MODS[i].m * (MODS[i].m + G - 1) * MODS[i].W * 8;
        fprintf(stderr, "\n  cumulative density %.3g (%.3f survivors per 512), tables %.1f MB\n", pr, pr * 512, tb / 1048576.0);
    }
}

/* ------------------------------------------------------------------ */
/* coprimality mask rows of period 210: pattern bit b set for prime {2,3,5,7}[b] dividing gcd(x1,x2) */
static uint64_t *MASK210; static uint32_t W210; static uint8_t TAB210[210];   /* TAB210[x] = pattern of {2,3,5,7} dividing x */
static void build_mask210(void) {
    W210 = row_words(210);
    MASK210 = aligned_alloc(64, (size_t)16 * W210 * 8);
    static const int pr[4] = {2, 3, 5, 7};
    for (int x = 0; x < 210; x++) { int p = 0; for (int b = 0; b < 4; b++) if (x % pr[b] == 0) p |= 1 << b; TAB210[x] = p; }
    uint8_t ok[210];
    for (int pat = 0; pat < 16; pat++) {
        for (int x = 0; x < 210; x++) { int o = 1; for (int b = 0; b < 4; b++) if ((pat >> b & 1) && x % pr[b] == 0) o = 0; ok[x] = o; }
        uint64_t *w = MASK210 + (size_t)pat * W210; memset(w, 0, (size_t)W210 * 8); pack_bits(ok, 210, w); extend_periodic(w, 210, W210);
    }
}

static int is_square_i128(i128 v, i128 *root) {
    if (v < 0) return 0;
    if (v == 0) { *root = 0; return 1; }
    i128 s = (i128)sqrtl((long double)v);
    while (s * s > v) s--;
    while ((s + 1) * (s + 1) <= v) s++;
    *root = s; return s * s == v;
}

/* mod-p prefilter for the GMP regime: quadratic-residue tables for four 16-bit primes */
static const uint32_t PRE[NPRE] = {65521, 65519, 65497, 65479};
static uint8_t *PREQR[NPRE]; static uint64_t PREMAGIC[NPRE]; static uint32_t SPRE[NPRE][28];
static void build_prefilter(void) {
    for (int j = 0; j < NPRE; j++) {
        uint32_t p = PRE[j]; PREQR[j] = calloc(p, 1); for (uint64_t v = 0; v < p; v++) PREQR[j][(v * v) % p] = 1;
        PREMAGIC[j] = fastmod_magic(p); for (int n = 0; n < 28; n++) SPRE[j][n] = modp(S[n], p);
    }
}

/* per-row lazily computed data */
typedef struct {
    int have_coef, have_gcd, have_pre, have_skip; int64_t g12;
    i128 c[7];                      /* int128 regime */
    mpz_t cz[7];                    /* GMP regime */
    uint32_t cm[NPRE][7];           /* coefficients mod the prefilter primes */
    int rowskip; int32_t E[256] __attribute__((aligned(64)));   /* skip lines: whole row on a line (A = 0); E_s = B_s x1 + C_s x2 */
} rowctx_t;

/* Points on the skip lines survive every modulus (S is a square on them), and lines with |A| = 1 cross almost every
 * row, so they are rejected before the gcd and the exact test: per row E_s = B_s x1 + C_s x2 (lazily, 16 lines per
 * zmm), per survivor "A_s x0 + E_s = 0 for some s".  32-bit lanes are exact while 3 * max|coef| * H < 2^31
 * (otherwise a scalar int64 loop is used). */
static int32_t LA[256] __attribute__((aligned(64))), LB[256] __attribute__((aligned(64))), LC[256] __attribute__((aligned(64)));
static int nlb = 0, lines_vec = 0;   /* 16-line blocks; vector path usable */
static void setup_lines(int64_t H) {
    nlb = (nskip + 15) / 16;
    int64_t maxc = 0;
    for (int s = 0; s < nlb * 16; s++) {
        if (s < nskip) { LA[s] = (int32_t)skipL[s][0]; LB[s] = (int32_t)skipL[s][1]; LC[s] = (int32_t)skipL[s][2]; for (int k = 0; k < 3; k++) { int64_t v = skipL[s][k] < 0 ? -skipL[s][k] : skipL[s][k]; if (v > maxc) maxc = v; } }
        else { LA[s] = 0; LB[s] = 0; LC[s] = 1; }   /* padding: never 0 as long as x2 != 0, handled below */
    }
    lines_vec = nskip > 0 && 3 * maxc * H < ((int64_t)1 << 31) && maxc < ((int64_t)1 << 31);
}
static void row_lines(rowctx_t *X, int64_t x1, int64_t x2) {
    X->rowskip = 0;
    if (lines_vec) {
        __m512i v1 = _mm512_set1_epi32((int32_t)x1), v2 = _mm512_set1_epi32((int32_t)x2), zero = _mm512_setzero_si512();
        for (int k = 0; k < nlb; k++) {
            __m512i e = _mm512_add_epi32(_mm512_mullo_epi32(_mm512_load_si512(LB + 16 * k), v1), _mm512_mullo_epi32(_mm512_load_si512(LC + 16 * k), v2));
            __mmask16 pad = (k == nlb - 1 && (nskip & 15)) ? (__mmask16)~((1u << (nskip & 15)) - 1) : 0;
            e = _mm512_mask_mov_epi32(e, pad, _mm512_set1_epi32(1));
            _mm512_store_si512(X->E + 16 * k, e);
            if (_mm512_mask_cmpeq_epi32_mask(_mm512_cmpeq_epi32_mask(_mm512_load_si512(LA + 16 * k), zero), e, zero)) X->rowskip = 1;
        }
    } else {
        for (int s = 0; s < nskip; s++) if (skipL[s][0] == 0 && skipL[s][1] * x1 + skipL[s][2] * x2 == 0) X->rowskip = 1;
    }
    X->have_skip = 1;
}
static inline int on_skipline(const rowctx_t *X, int64_t x0, int64_t x1, int64_t x2) {
    if (lines_vec) {
        __m512i v0 = _mm512_set1_epi32((int32_t)x0), zero = _mm512_setzero_si512(); __mmask16 hit = 0;
        for (int k = 0; k < nlb; k++) hit |= _mm512_cmpeq_epi32_mask(_mm512_add_epi32(_mm512_mullo_epi32(_mm512_load_si512(LA + 16 * k), v0), _mm512_load_si512(X->E + 16 * k)), zero);
        return hit != 0;
    }
    for (int s = 0; s < nskip; s++) if (skipL[s][0] * x0 + skipL[s][1] * x1 + skipL[s][2] * x2 == 0) return 1;
    return 0;
}

static void mpz_set_i128(mpz_t r, i128 v) {
    int neg = v < 0; u128 a = neg ? -(u128)v : (u128)v;
    mpz_set_ui(r, (unsigned long)(a >> 64)); mpz_mul_2exp(r, r, 64); mpz_add_ui(r, r, (unsigned long)a);
    if (neg) mpz_neg(r, r);
}
/* coefficients c_i(x1,x2) of the sextic in x0 */
static void row_coefficients(rowctx_t *X, int64_t x1, int64_t x2, int use128, mpz_t tmp) {
    i128 p1[7], p2[7]; p1[0] = 1; p2[0] = 1; for (int e = 1; e <= 6; e++) { p1[e] = p1[e-1] * x1; p2[e] = p2[e-1] * x2; }
    if (use128) {
        for (int i = 0; i <= 6; i++) X->c[i] = 0;
        for (int n = 0; n < 28; n++) X->c[SI[n]] += (i128)S[n] * p1[SJ[n]] * p2[SK[n]];
    } else {   /* x1^j x2^k fits in int128 for H < 2^21; the products with S_n may not */
        for (int i = 0; i <= 6; i++) mpz_set_ui(X->cz[i], 0);
        for (int n = 0; n < 28; n++) { mpz_set_i128(tmp, p1[SJ[n]] * p2[SK[n]]); mpz_mul_si(tmp, tmp, S[n]); mpz_add(X->cz[SI[n]], X->cz[SI[n]], tmp); }
    }
    X->have_coef = 1;
}
static void row_prefilter(rowctx_t *X, int64_t x1, int64_t x2) {
    for (int j = 0; j < NPRE; j++) {
        uint32_t p = PRE[j]; uint64_t M = PREMAGIC[j];
        uint32_t pw1[7], pw2[7]; pw1[0] = pw2[0] = 1; uint32_t x1p = modp(x1, p), x2p = modp(x2, p);
        for (int e = 1; e <= 6; e++) { pw1[e] = fastmod_u32(pw1[e-1] * x1p, M, p); pw2[e] = fastmod_u32(pw2[e-1] * x2p, M, p); }
        uint32_t c[7] = {0,0,0,0,0,0,0};
        for (int n = 0; n < 28; n++) { uint32_t t = fastmod_u32(SPRE[j][n] * pw1[SJ[n]], M, p); c[SI[n]] += fastmod_u32(t * pw2[SK[n]], M, p); }
        for (int i = 0; i <= 6; i++) X->cm[j][i] = fastmod_u32(c[i], M, p);
    }
    X->have_pre = 1;
}
static inline int prefilter_pass(const rowctx_t *X, int64_t x0) {
    for (int j = 0; j < NPRE; j++) {
        uint32_t p = PRE[j]; uint64_t M = PREMAGIC[j]; uint32_t x = modp(x0, p), v = X->cm[j][6];
        for (int i = 5; i >= 0; i--) v = fastmod_u32(v * x + X->cm[j][i], M, p);
        if (!PREQR[j][v]) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
typedef struct { char *out; size_t olen, ocap; uint64_t ncand, nsurv, npts; mpz_t mval, mtmp; } thr_t;

static void flush_out(thr_t *T) {
    #pragma omp critical
    { fwrite(T->out, 1, T->olen, stdout); fflush(stdout); }
    T->olen = 0;
}
static inline void emit(thr_t *T, int64_t x0, int64_t x1, int64_t x2, int iszero, int quiet) {
    T->npts++;
    if (quiet) return;
    if (T->olen + 80 > T->ocap) { if (T->ocap >= (1u << 22)) flush_out(T); else { T->ocap *= 2; T->out = realloc(T->out, T->ocap); } }
    T->olen += sprintf(T->out + T->olen, "%ld %ld %ld%s\n", x0, x1, x2, iszero ? " W" : "");
}


/* value-taking options must have an argument and numeric values must parse completely */
static char *arg_str(int argc, char **argv, int *i) { if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", argv[*i]); exit(1); } return argv[++*i]; }
static __attribute__((unused)) long long arg_ll(int argc, char **argv, int *i) { const char *opt = argv[*i]; char *s = arg_str(argc, argv, i), *end = NULL; errno = 0; long long v = strtoll(s, &end, 10); if (errno == ERANGE || end == s || *end != 0) { fprintf(stderr, "option %s: invalid integer '%s'\n", opt, s); exit(1); } return v; }
static __attribute__((unused)) double arg_d(int argc, char **argv, int *i) { const char *opt = argv[*i]; char *s = arg_str(argc, argv, i), *end = NULL; errno = 0; double v = strtod(s, &end); if (errno == ERANGE || end == s || *end != 0) { fprintf(stderr, "option %s: invalid number '%s'\n", opt, s); exit(1); } return v; }

/* integer options/arguments are validated as long long against their range BEFORE any narrowing cast */
static long long chk_ll(const char *what, const char *s, long long lo, long long hi) { char *end = NULL; errno = 0; long long v = strtoll(s, &end, 10); if (errno != 0 || end == s || *end != 0) { fprintf(stderr, "%s: invalid integer '%s'\n", what, s); exit(1); } if (v < lo || v > hi) { fprintf(stderr, "%s: %lld is outside [%lld, %lld]\n", what, v, lo, hi); exit(1); } return v; }
static int arg_int(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return (int)chk_ll(opt, argv[++*i], lo, hi); }
static __attribute__((unused)) long long arg_llr(int argc, char **argv, int *i, long long lo, long long hi) { const char *opt = argv[*i]; if (*i + 1 >= argc) { fprintf(stderr, "option %s needs an argument\n", opt); exit(1); } return chk_ll(opt, argv[++*i], lo, hi); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: k3pts sextic.txt H [-t threads] [-skip lines.txt] [-q] [-v] [-Q qmax] [-M mcap] [-nmod n] [-range x2lo x2hi] [-nosym]\n"); return 1; }
    int64_t H = chk_ll("H", argv[2], 1, 1000000000); int nthreads = omp_get_max_threads(); int quiet = 0; char *skipfile = NULL;
    uint32_t Qmax = 256, Mcap = 512; int nmodmax = 40; int64_t x2lo = 0, x2hi = H; int nosym = 0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "-t")) nthreads = arg_int(argc, argv, &i, 1, 4096);
        else if (!strcmp(argv[i], "-skip")) skipfile = arg_str(argc, argv, &i);
        else if (!strcmp(argv[i], "-q")) quiet = 1;
        else if (!strcmp(argv[i], "-v")) verbose++;
        else if (!strcmp(argv[i], "-Q")) Qmax = (uint32_t)arg_int(argc, argv, &i, 2, 1000000);
        else if (!strcmp(argv[i], "-M")) Mcap = (uint32_t)arg_int(argc, argv, &i, 1, 1000000);
        else if (!strcmp(argv[i], "-nmod")) nmodmax = arg_int(argc, argv, &i, 1, 1000);
        else if (!strcmp(argv[i], "-range")) { x2lo = arg_llr(argc, argv, &i, 0, 1000000000); x2hi = arg_llr(argc, argv, &i, 0, 1000000000); }
        else if (!strcmp(argv[i], "-nosym")) nosym = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (nmodmax > MAXMOD) nmodmax = MAXMOD;
    if (x2lo < 0) x2lo = 0; if (x2hi > H) x2hi = H;
    omp_set_num_threads(nthreads);
    { FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; } char tok[128]; int n = 0;   /* exactly 28 complete int64 coefficients, nothing else */
      while (fscanf(f, "%127s", tok) == 1) { if (n == 28) { fprintf(stderr, "%s: more than 28 coefficients\n", argv[1]); return 1; } char *end = NULL; errno = 0; long long v = strtoll(tok, &end, 10); if (errno != 0 || end == tok || *end != 0) { fprintf(stderr, "%s: coefficient %d is not an int64 integer: '%s'\n", argv[1], n + 1, tok); return 1; } S[n++] = v; }
      if (n != 28) { fprintf(stderr, "need 28 coefficients, got %d\n", n); return 1; } fclose(f); }
    { int n = 0; for (int i = 6; i >= 0; i--) for (int j = 6 - i; j >= 0; j--) { SI[n] = i; SJ[n] = j; SK[n] = 6 - i - j; n++; } }
    if (skipfile) {   /* a skip record is exactly three complete int64 integers (comment/blank lines skipped); anything else aborts */
      FILE *f = fopen(skipfile, "r"); if (!f) { perror(skipfile); return 1; } char *buf = NULL; size_t bcap = 0; ssize_t blen; long rec = 0;
      while ((blen = getline(&buf, &bcap, f)) > 0) { rec++; if (buf[0] == '#') continue; char *save = NULL, *tok = strtok_r(buf, " \t\r\n", &save); if (!tok) continue; int64_t v[3]; int nt = 0;
        for (; tok; tok = strtok_r(NULL, " \t\r\n", &save)) { if (nt == 3) { fprintf(stderr, "skip file %s: record %ld has more than 3 fields\n", skipfile, rec); return 1; } char *e = NULL; errno = 0; long long x = strtoll(tok, &e, 10); if (errno != 0 || e == tok || *e != 0) { fprintf(stderr, "skip file %s: record %ld, field %d is not an int64 integer: '%s'\n", skipfile, rec, nt + 1, tok); return 1; } v[nt++] = x; }
        if (nt != 3) { fprintf(stderr, "skip file %s: record %ld has %d fields, need 3 (A B C)\n", skipfile, rec, nt); return 1; }
        if (nskip == 256) { fprintf(stderr, "skip file %s has more than 256 lines (the skip table is fixed at 256)\n", skipfile); return 1; }
        skipL[nskip][0] = v[0]; skipL[nskip][1] = v[1]; skipL[nskip][2] = v[2]; nskip++; }
      free(buf); fclose(f); }
    setup_lines(H);
    /* symmetry x0 -> -x0: no odd powers of x0 */
    int sym = !nosym; for (int n = 0; n < 28; n++) if (SI[n] & 1 && S[n] != 0) sym = 0;
    /* exact-test regime: all intermediate values are bounded by sum |S_n| H^6 */
    long double bound = 0; for (int n = 0; n < 28; n++) bound += fabsl((long double)S[n]);
    bound *= powl((long double)H, 6);
    int use128 = bound < ldexpl(1.0L, 125) && H < (1 << 21);
    if (getenv("K3PTS_GMP")) use128 = 0;   /* testing knob: force the GMP regime */
    /* survivors cost ~300 cycles (int128) or ~2000 (GMP, after the prefilter) against ~3 per 512 candidates per modulus */
    double stop_density = (use128 ? 3.0 / 300 : 3.0 / 2000) / 512;

    double t0 = omp_get_wtime();
    choose_moduli(Qmax, Mcap, nmodmax, stop_density);
    build_mask210(); build_prefilter();
    double t1 = omp_get_wtime();
    if (!quiet) fprintf(stderr, "tables built in %.2fs; H=%ld, %d threads, %d skip lines%s%s\n", t1 - t0, H, nthreads, nskip, sym ? ", x0 -> -x0 symmetry" : "", use128 ? "" : ", GMP regime");

    /* x0 range and per-chunk phases (they depend only on the chunk, not on the row) */
    int64_t x0lo = sym ? 0 : -H; int64_t nx0 = H - x0lo + 1; int nchunks = (int)((nx0 + 511) / 512);
    uint32_t *OFF = malloc((size_t)nchunks * MAXMOD * sizeof(uint32_t)), *SH = malloc((size_t)nchunks * MAXMOD * sizeof(uint32_t)), *OFF210 = malloc(nchunks * sizeof(uint32_t)), *SH210 = malloc(nchunks * sizeof(uint32_t));
    for (int c = 0; c < nchunks; c++) {
        int64_t x0b = x0lo + 512 * (int64_t)c;
        for (int mi = 0; mi < nmod; mi++) { uint32_t ph = modp(x0b, MODS[mi].m); OFF[(size_t)c * MAXMOD + mi] = (ph >> 6) * 8; SH[(size_t)c * MAXMOD + mi] = ph & 63; }
        uint32_t ph = modp(x0b, 210); OFF210[c] = (ph >> 6) * 8; SH210[c] = ph & 63;
    }
    uint64_t ncand = 0, nsurv = 0, npts = 0;
    #pragma omp parallel reduction(+:ncand, nsurv, npts)
    {
        thr_t T; T.ocap = 1 << 20; T.out = malloc(T.ocap); T.olen = 0; T.ncand = T.nsurv = T.npts = 0; mpz_init(T.mval); mpz_init(T.mtmp);
        rowctx_t X[G]; for (int g = 0; g < G; g++) for (int i = 0; i <= 6; i++) mpz_init(X[g].cz[i]);
        const uint8_t *blk[MAXMOD], *base[MAXMOD]; size_t rb[MAXMOD]; uint32_t x1r[MAXMOD];
        for (int mi = 0; mi < nmod; mi++) rb[mi] = (size_t)MODS[mi].W * 8;
        #pragma omp for schedule(dynamic, 1)
        for (int64_t x2 = x2lo; x2 <= x2hi; x2++) {
            int64_t x1lo = x2 == 0 ? 1 : -H;
            int pat2 = TAB210[modp(x2, 210)]; uint32_t x1m = modp(x1lo, 210);
            for (int mi = 0; mi < nmod; mi++) { mod_t *V = &MODS[mi]; blk[mi] = (const uint8_t*)(V->T + (size_t)modp(x2, V->m) * (V->m + G - 1) * V->W); x1r[mi] = modp(x1lo, V->m); }
            for (int64_t x1g = x1lo; x1g <= H; x1g += G) {
                int ng = H - x1g + 1 < G ? (int)(H - x1g + 1) : G;
                int pat[G];
                for (int g = 0; g < G; g++) {
                    uint32_t xm = x1m + g; if (xm >= 210) xm -= 210;
                    pat[g] = TAB210[xm] & pat2; X[g].have_coef = X[g].have_gcd = X[g].have_pre = X[g].have_skip = 0;
                }
                x1m += G; if (x1m >= 210) x1m -= 210;
                for (int mi = 0; mi < nmod; mi++) {
                    base[mi] = blk[mi] + (size_t)x1r[mi] * rb[mi];
                    x1r[mi] += G; while (x1r[mi] >= MODS[mi].m) x1r[mi] -= MODS[mi].m;
                }
                for (int c = 0; c < nchunks; c++) {
                    __m512i acc[G];
                    for (int g = 0; g < G; g++) {
                        if (pat[g]) { const uint8_t *r = (const uint8_t*)(MASK210 + (size_t)pat[g] * W210) + OFF210[c]; acc[g] = _mm512_shrdv_epi64(_mm512_loadu_si512(r), _mm512_loadu_si512(r + 8), _mm512_set1_epi64(SH210[c])); }
                        else acc[g] = _mm512_set1_epi64(-1);
                    }
                    const uint32_t *off = OFF + (size_t)c * MAXMOD, *sh = SH + (size_t)c * MAXMOD;
                    for (int mi = 0; mi < nmod; mi++) {
                        const uint8_t *r = base[mi] + off[mi]; size_t RB = rb[mi]; __m512i cnt = _mm512_set1_epi64(sh[mi]);
                        for (int g = 0; g < G; g++) acc[g] = _mm512_and_si512(acc[g], _mm512_shrdv_epi64(_mm512_loadu_si512(r + g * RB), _mm512_loadu_si512(r + g * RB + 8), cnt));
                    }
                    int64_t x0b = x0lo + 512 * (int64_t)c;
                    int64_t rem = H - x0b + 1; if (rem > 512) rem = 512;
                    T.ncand += (uint64_t)rem * ng;
                    for (int g = 0; g < ng; g++) {
                        if (_mm512_test_epi64_mask(acc[g], acc[g]) == 0) continue;
                        uint64_t w[8] __attribute__((aligned(64))); _mm512_store_si512(w, acc[g]);
                        int64_t x1 = x1g + g; rowctx_t *R = &X[g];
                        for (int l = 0; l < 8; l++) {
                            uint64_t x = w[l]; int64_t xb = x0b + 64 * l;
                            if (xb > H) break;
                            if (xb + 63 > H) x &= (UINT64_C(1) << (H - xb + 1)) - 1;
                            while (x) {
                                int t = __builtin_ctzll(x); x &= x - 1; int64_t x0 = xb + t; T.nsurv++;
                                if (!R->have_skip) row_lines(R, x1, x2);
                                if (R->rowskip) continue;
                                int skip0 = on_skipline(R, x0, x1, x2), skip1 = (sym && x0 != 0) ? on_skipline(R, -x0, x1, x2) : 1;
                                if (skip0 && skip1) continue;
                                if (!R->have_gcd) { R->g12 = gcd64(x1, x2); R->have_gcd = 1; }
                                if (R->g12 > 1 && gcd64(R->g12, x0) != 1) continue;
                                int issq = 0, iszero = 0;
                                if (use128) {
                                    if (!R->have_coef) row_coefficients(R, x1, x2, 1, T.mtmp);
                                    i128 v = R->c[6]; for (int i = 5; i >= 0; i--) v = v * x0 + R->c[i];
                                    i128 root; if (v == 0) { iszero = 1; issq = 1; } else issq = is_square_i128(v, &root);
                                } else {
                                    if (!R->have_pre) row_prefilter(R, x1, x2);
                                    if (!prefilter_pass(R, x0)) continue;
                                    if (!R->have_coef) row_coefficients(R, x1, x2, 0, T.mtmp);
                                    mpz_set(T.mval, R->cz[6]); for (int i = 5; i >= 0; i--) { mpz_mul_si(T.mval, T.mval, x0); mpz_add(T.mval, T.mval, R->cz[i]); }
                                    if (mpz_sgn(T.mval) == 0) { iszero = 1; issq = 1; } else if (mpz_sgn(T.mval) > 0) issq = mpz_perfect_square_p(T.mval);
                                }
                                if (!issq) continue;
                                if (!skip0) emit(&T, x0, x1, x2, iszero, quiet);
                                if (!skip1) emit(&T, -x0, x1, x2, iszero, quiet);
                            }
                        }
                    }
                }
            }
        }
        flush_out(&T);
        ncand += T.ncand; nsurv += T.nsurv; npts += T.npts;
        free(T.out); mpz_clear(T.mval); mpz_clear(T.mtmp); for (int g = 0; g < G; g++) for (int i = 0; i <= 6; i++) mpz_clear(X[g].cz[i]);
    }
    double t2 = omp_get_wtime();
    fprintf(stderr, "H=%ld: %lu candidates, %lu sieve survivors, %lu points (excluding %d lines); %.2fs\n", H, ncand, nsurv, npts, nskip, t2 - t1);
    return 0;
}

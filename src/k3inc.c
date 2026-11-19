/* k3inc.c -- incidence counting between the rational points of a K3 double plane (from k3pts) and all lines
 * A x0 + B x1 + C x2 = 0 with max(|A|,|B|,|C|) <= L.
 *
 * Stage A: for every point of height <= HA, enumerate the lattice of lines through it inside the box and count
 *          incidences in a hash table; keep lines with >= kA points.
 * Stage B: for every point of height <= HB (all points in the file), enumerate its lines again and count only the
 *          stage-A candidates; print candidates with >= kB points, sorted by the stage-B count.
 *
 *   k3inc points.txt L HA kA kB [-t threads] [-slots 2^n] [-A Alo Ahi]
 * points file: "x0 x1 x2 [W]" per line (k3pts output).  Output: "countB countA A B C".
 * -A Alo Ahi restricts both stages to the lines with Alo <= A <= Ahi (0 <= Alo <= Ahi <= L; the sign representative has
 * A >= 0, so the default 0..L is the whole box): a box too large for one table is processed in slabs, each slab a run whose
 * table holds only its own lines, and the outputs concatenated (the lines of the slabs are disjoint; k3inc_select.py does
 * not need a sorted input).  The lattice enumeration narrows its range to the slab, so a slab costs about its share of
 * the box plus a per-point overhead of one interval computation per lattice row.
 *
 * Hash table: open addressing with linear probing over 2^n slots.  For L < 2^14 the key (A, B+L, C+L, 3 x 15 bits at
 * most) and the count share one 64-bit word per slot ("packed": one cache line per increment, half the memory); for
 * larger L the table falls back to a 64-bit key array plus a 32-bit count array.  Increments are batched per thread
 * (32 keys): the slots are prefetched first, then the atomic updates are applied, which hides most of the DRAM latency
 * of the random accesses.  The output is sorted by (countB, countA) descending, then (A, B, C) ascending, so it is
 * reproducible (the slot order used before depended on probing races between threads).
 */
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/mman.h>
#include "simd.h"

/* zero-filled table on (transparent) huge pages when available: fewer TLB misses on the random accesses */
static void *table_alloc(size_t bytes) {
    void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) return NULL;
    madvise(p, bytes, MADV_HUGEPAGE);
    return p;
}

static int64_t gcd64(int64_t a, int64_t b) {   /* binary gcd */
    uint64_t x = a < 0 ? -(uint64_t)a : (uint64_t)a, y = b < 0 ? -(uint64_t)b : (uint64_t)b;
    if (!x) return (int64_t)y; if (!y) return (int64_t)x;
    int sh = __builtin_ctzll(x | y); x >>= __builtin_ctzll(x);
    do { y >>= __builtin_ctzll(y); if (x > y) { uint64_t t = x; x = y; y = t; } y -= x; } while (y);
    return (int64_t)(x << sh);
}
static int64_t egcd(int64_t a, int64_t b, int64_t *x, int64_t *y) { if (b == 0) { *x = 1; *y = 0; return a; } int64_t x1, y1; int64_t g = egcd(b, a % b, &x1, &y1); *x = y1; *y = x1 - (a / b) * y1; return g; }

/* integer basis U, V of {(A,B,C): A p0 + B p1 + C p2 = 0} for primitive (p0,p1,p2) */
static void dual_basis(int64_t p0, int64_t p1, int64_t p2, int64_t U[3], int64_t V[3]) {
    if (p1 == 0 && p2 == 0) { U[0] = 0; U[1] = 1; U[2] = 0; V[0] = 0; V[1] = 0; V[2] = 1; return; }
    int64_t x, y; int64_t g = egcd(p1 < 0 ? -p1 : p1, p2 < 0 ? -p2 : p2, &x, &y); if (g < 0) g = -g;
    V[0] = 0; V[1] = p2 / g; V[2] = -p1 / g;
    int64_t u, v; int64_t g2 = egcd((p1 / g) < 0 ? -(p1 / g) : p1 / g, (p2 / g) < 0 ? -(p2 / g) : p2 / g, &u, &v);
    if (g2 < 0) { u = -u; v = -v; }
    if (p1 < 0) u = -u;
    if (p2 < 0) v = -v;
    U[0] = g; U[1] = -p0 * u; U[2] = -p0 * v;
}

/* Gauss reduction of a 2-D lattice basis (Euclidean norm) */
static void reduce2(int64_t U[3], int64_t V[3]) {
    for (int it = 0; it < 64; it++) {
        double nu = (double)U[0]*U[0] + (double)U[1]*U[1] + (double)U[2]*U[2], nv = (double)V[0]*V[0] + (double)V[1]*V[1] + (double)V[2]*V[2];
        if (nu > nv) { for (int i = 0; i < 3; i++) { int64_t t = U[i]; U[i] = V[i]; V[i] = t; } double t = nu; nu = nv; nv = t; }
        double d = (double)U[0]*V[0] + (double)U[1]*V[1] + (double)U[2]*V[2]; int64_t q = (int64_t)llround(d / nu);
        if (q == 0) break;
        for (int i = 0; i < 3; i++) V[i] -= q * U[i];
    }
}

#define KEY(A,B,C,L) ((((uint64_t)((A) + (L))) << 42) | (((uint64_t)((B) + (L))) << 21) | (uint64_t)((C) + (L)))
static inline uint64_t hash64(uint64_t k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k; }

/* ---- hash table ---- */
static uint64_t *hkeys; static uint32_t *hcnt; static uint64_t hmask;   /* two-array layout: key+1 (0 = empty), count */
static int packed = 0; static int sbits = 0; static uint64_t kmask, CNT1, cntmax; static int64_t GL;   /* packed layout */

/* packed key from the 63-bit KEY: (A, B+L, C+L) in 3 x sbits bits, plus 1 (0 = empty slot) */
static inline uint64_t compact(uint64_t key) {
    uint64_t a = (key >> 42) - (uint64_t)GL, b = (key >> 21) & 0x1FFFFF, c = key & 0x1FFFFF;
    return ((a << (2 * sbits)) | (b << sbits) | c) + 1;
}
static inline void decode(uint64_t ck, int64_t *A, int64_t *B, int64_t *C) {
    ck -= 1; uint64_t sm = ((uint64_t)1 << sbits) - 1;
    *A = (int64_t)(ck >> (2 * sbits)); *B = (int64_t)((ck >> sbits) & sm) - GL; *C = (int64_t)(ck & sm) - GL;
}

static uint64_t hprobe_limit = 0;   /* = number of slots; a probe sequence longer than this means the table is full */
static void htable_full(void) { fprintf(stderr, "k3inc: hash table full (more distinct lines than slots) -- rerun with a larger -slots\n"); exit(2); }
/* stage A: insert-or-increment at the probe start h */
static inline void hinc(uint64_t key, uint64_t h) {
    if (packed) {
        uint64_t want = compact(key); uint64_t probes = 0;
        for (;;) {
            if (++probes > hprobe_limit) htable_full();
            uint64_t cur = __atomic_load_n(&hkeys[h], __ATOMIC_RELAXED);
            if ((cur & kmask) == want) { if ((cur >> (3 * sbits)) < cntmax) __atomic_fetch_add(&hkeys[h], CNT1, __ATOMIC_RELAXED); return; }
            if (cur == 0) {
                uint64_t expect = 0;
                if (__atomic_compare_exchange_n(&hkeys[h], &expect, want | CNT1, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return;
                if ((expect & kmask) == want) { __atomic_fetch_add(&hkeys[h], CNT1, __ATOMIC_RELAXED); return; }
            }
            h = (h + 1) & hmask;
        }
    }
    uint64_t probes = 0;
    for (;;) {
        if (++probes > hprobe_limit) htable_full();
        uint64_t cur = __atomic_load_n(&hkeys[h], __ATOMIC_RELAXED);
        if (cur == key + 1) { __atomic_fetch_add(&hcnt[h], 1, __ATOMIC_RELAXED); return; }
        if (cur == 0) { uint64_t expect = 0; if (__atomic_compare_exchange_n(&hkeys[h], &expect, key + 1, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) { __atomic_fetch_add(&hcnt[h], 1, __ATOMIC_RELAXED); return; } if (expect == key + 1) { __atomic_fetch_add(&hcnt[h], 1, __ATOMIC_RELAXED); return; } }
        h = (h + 1) & hmask;
    }
}
/* stage B: increment only if present */
static inline void hinc_if(uint64_t key, uint64_t h) {
    if (packed) {
        uint64_t want = compact(key); uint64_t probes = 0;
        for (;;) { if (++probes > hprobe_limit) return; uint64_t cur = hkeys[h]; if ((cur & kmask) == want) { if ((cur >> (3 * sbits)) < cntmax) __atomic_fetch_add(&hkeys[h], CNT1, __ATOMIC_RELAXED); return; } if (cur == 0) return; h = (h + 1) & hmask; }
    }
    { uint64_t probes = 0; for (;;) { if (++probes > hprobe_limit) return; uint64_t cur = hkeys[h]; if (cur == key + 1) { __atomic_fetch_add(&hcnt[h], 1, __ATOMIC_RELAXED); return; } if (cur == 0) return; h = (h + 1) & hmask; } }
}

/* per-thread batching of increments: prefetch the slots of a batch, then apply */
#define KB 32
typedef struct { uint64_t key[KB], h[KB]; int n; } kbuf_t;
static __thread kbuf_t KBUF;
static int stage = 0;
static int64_t ALO_ = 0;   /* copy of the slab bound ALO for the direct table (ALO is declared below) */
static void kflush(void) {
    kbuf_t *b = &KBUF;
    for (int i = 0; i < b->n; i++) { b->h[i] = hash64(b->key[i]) & hmask; __builtin_prefetch(&hkeys[b->h[i]], 1, 0); if (!packed) __builtin_prefetch(&hcnt[b->h[i]], 1, 0); }
    if (stage == 0) for (int i = 0; i < b->n; i++) hinc(b->key[i], b->h[i]);
    else for (int i = 0; i < b->n; i++) hinc_if(b->key[i], b->h[i]);
    b->n = 0;
}
static void kadd(uint64_t key) { kbuf_t *b = &KBUF; b->key[b->n++] = key; if (b->n == KB) kflush(); }

/* ---- direct table ----
 * When every point is counted in stage A and the box is small, the lines of the box index a table of one byte per line
 * directly: line (A, B, C) is entry ((A - ALO) (2L+1) + B + L) (2L+1) + C + L.  Each thread counts in a table of its own
 * (plain increments, no probing, an eighth of the memory of the hash table per thread), and the tables are added at the
 * end.  A byte stops at 255; what a line receives beyond that is counted in a small table of the thread (lines with so
 * many points are few: the lines through a singular point of the sextic), so the counts are exact.  If that small table
 * fills up, the direct count is abandoned and the hash table is used. */
#define OVN 4096                 /* slots of the overflow table of a thread */
typedef struct { uint64_t key[OVN]; uint32_t cnt[OVN]; int used; } ovt_t;   /* key = entry + 1 (0 = empty) */
static __thread uint8_t *DTAB; static __thread ovt_t *DOV; static int doverflow = 0; static int64_t DW;   /* DW = 2L+1 */
static void ov_inc(ovt_t *o, uint64_t h) {
    uint64_t g = hash64(h) & (OVN - 1);
    for (;;) { if (o->key[g] == h + 1) { o->cnt[g]++; return; } if (o->key[g] == 0) { if (o->used >= OVN / 2) { doverflow = 1; return; } o->key[g] = h + 1; o->cnt[g] = 1; o->used++; return; } g = (g + 1) & (OVN - 1); }
}
static uint32_t ov_get(const ovt_t *o, uint64_t h) {
    uint64_t g = hash64(h) & (OVN - 1);
    for (;;) { if (o->key[g] == h + 1) return o->cnt[g]; if (o->key[g] == 0) return 0; g = (g + 1) & (OVN - 1); }
}
static inline uint64_t didx(uint64_t key) { int64_t A = (int64_t)(key >> 42) - GL, B = (int64_t)((key >> 21) & 0x1FFFFF), C = (int64_t)(key & 0x1FFFFF); return (uint64_t)(((A - ALO_) * DW + B) * DW + C); }
static void dflush(void) {
    kbuf_t *b = &KBUF; uint8_t *t = DTAB;
    for (int i = 0; i < b->n; i++) { b->h[i] = didx(b->key[i]); __builtin_prefetch(&t[b->h[i]], 1, 0); }
    for (int i = 0; i < b->n; i++) { uint8_t v = t[b->h[i]]; if (v == 255) ov_inc(DOV, b->h[i]); else t[b->h[i]] = v + 1; }
    b->n = 0;
}
static void dadd(uint64_t key) { kbuf_t *b = &KBUF; b->key[b->n++] = key; if (b->n == KB) dflush(); }

/* enumerate the primitive lines through P in the box, calling f(key) */
typedef void (*linecb)(uint64_t);
static int64_t ALO = 0, AHI = -1;   /* slab: lines with ALO <= A <= AHI (set in main; default 0..L) */
static uint8_t COP[256][256];   /* COP[a][b] = 1 iff gcd(a, b) = 1 (filled in main) */
static inline int coprime(int64_t s, int64_t t) { uint64_t a = s < 0 ? -(uint64_t)s : (uint64_t)s, b = t < 0 ? -(uint64_t)t : (uint64_t)t; return (a | b) < 256 ? COP[a][b] : gcd64(s, t) == 1; }
static void lines_through(int64_t p0, int64_t p1, int64_t p2, int64_t L, linecb f) {
    int64_t U[3], V[3]; dual_basis(p0, p1, p2, U, V); reduce2(U, V);
    double nu = sqrt((double)U[0]*U[0] + (double)U[1]*U[1] + (double)U[2]*U[2]), nv = sqrt((double)V[0]*V[0] + (double)V[1]*V[1] + (double)V[2]*V[2]);
    int64_t T = (int64_t)(2.0 * sqrt(3.0) * L / nv) + 1; (void)nu;
    for (int64_t t = -T; t <= T; t++) {
        double lo = -1e300, hi = 1e300; int ok = 1;
        for (int i = 0; i < 3; i++) {
            double tv = (double)t * V[i], bl = i == 0 ? (double)ALO : -(double)L, bh = i == 0 ? (double)AHI : (double)L;   /* coordinate bounds: the slab for A, the box for B, C */
            if (U[i] == 0) { if (tv < bl || tv > bh) { ok = 0; break; } continue; }
            double a = (bl - tv) / U[i], b = (bh - tv) / U[i]; if (a > b) { double x = a; a = b; b = x; }
            if (a > lo) lo = a; if (b < hi) hi = b;
        }
        if (!ok) continue;   /* no test of lo > hi here: see dpoint; an empty row has s0 > s1 */
        int64_t s0 = (int64_t)ceil(lo - 1e-9), s1 = (int64_t)floor(hi + 1e-9);
        for (int64_t s = s0; s <= s1; s++) {
            int64_t A = s * U[0] + t * V[0], B = s * U[1] + t * V[1], C = s * U[2] + t * V[2];
            if (A < ALO || A > AHI || B < -L || B > L || C < -L || C > L) continue;
            if (A < 0 || (A == 0 && B < 0) || (A == 0 && B == 0 && C <= 0)) continue;   /* each line once: canonical sign only (and not the zero vector) */
            if (!coprime(s, t)) continue;   /* U, V span the saturated lattice of lines through P, so gcd(A, B, C) = gcd(s, t) */
            f(KEY(A, B, C, L));
        }
    }
}


/* The lines through P for the direct table, whole box (ALO = 0, AHI = L).  With D = 2L+1 the entry of the line (A, B, C) in
 * canonical sign is A D^2 + B D + C + L D + L, and for either sign of (A, B, C) it is |A D^2 + B D + C| + L D + L (the
 * digits B, C are at most L in absolute value, so the sign of the number is the sign of its first nonzero digit).  The
 * number is linear in (s, t) for the line s U + t V, so a row of the lattice is an arithmetic progression of entries.
 * Each line is met once: t > 0, or t = 0 and s = 1.  The values of s prime to t are read from a table of bit masks. */
static uint64_t CM[256][8];   /* CM[t]: bit s + 256 is set iff gcd(s, t) = 1, for -256 <= s <= 255 (filled in main) */
static int64_t DOFF;          /* L D + L */
static void dflush2(void) {
    kbuf_t *b = &KBUF; uint8_t *t = DTAB;
    for (int i = 0; i < b->n; i++) __builtin_prefetch(&t[b->h[i]], 1, 0);
    for (int i = 0; i < b->n; i++) { uint8_t v = t[b->h[i]]; if (v == 255) ov_inc(DOV, b->h[i]); else t[b->h[i]] = v + 1; }
    b->n = 0;
}
static inline void dput(int64_t lin) { kbuf_t *b = &KBUF; b->h[b->n++] = (uint64_t)(lin < 0 ? -lin : lin) + (uint64_t)DOFF; if (b->n == KB) dflush2(); }
static inline int inbox(const int64_t *U, const int64_t *V, int64_t s, int64_t t, int64_t L) {
    for (int i = 0; i < 3; i++) { int64_t x = s * U[i] + t * V[i]; if (x < -L || x > L) return 0; }
    return 1;
}
static void dpoint(int64_t p0, int64_t p1, int64_t p2, int64_t L) {
    int64_t U[3], V[3]; dual_basis(p0, p1, p2, U, V); reduce2(U, V);
    double nv = sqrt((double)V[0]*V[0] + (double)V[1]*V[1] + (double)V[2]*V[2]);
    int64_t T = (int64_t)(2.0 * sqrt(3.0) * L / nv) + 1;
    int64_t stepU = (U[0] * DW + U[1]) * DW + U[2], stepV = (V[0] * DW + V[1]) * DW + V[2];
    double iu[3]; for (int i = 0; i < 3; i++) iu[i] = U[i] ? 1.0 / (double)U[i] : 0;
    if (inbox(U, V, 1, 0, L)) dput(stepU);
    for (int64_t t = 1; t <= T; t++) {
        double lo = -1e300, hi = 1e300; int ok = 1;
        for (int i = 0; i < 3; i++) {
            double tv = (double)t * V[i];
            if (U[i] == 0) { if (tv < -(double)L || tv > (double)L) { ok = 0; break; } continue; }
            double a = (-(double)L - tv) * iu[i], b = ((double)L - tv) * iu[i]; if (a > b) { double x = a; a = b; b = x; }
            if (a > lo) lo = a; if (b < hi) hi = b;
        }
        /* lo > hi is not tested before the tolerance: when the row meets the box in one point, a corner where two coordinates
         * of the line are +-L, lo and hi are the same integer up to the rounding of the products with iu[], and they can come
         * out as 1 and 1 - 2^-53.  An empty row has s0 > s1 below. */
        if (!ok) continue;
        int64_t s0 = (int64_t)ceil(lo - 1e-6), s1 = (int64_t)floor(hi + 1e-6);
        while (s0 <= s1 && !inbox(U, V, s0, t, L)) s0++;        /* the box is convex: exact ends give the exact row */
        while (s1 >= s0 && !inbox(U, V, s1, t, L)) s1--;
        if (s0 > s1) continue;
        int64_t base = t * stepV;
        if (t < 256 && s0 >= -256 && s1 <= 255) {
            const uint64_t *cm = CM[t]; int b0 = (int)(s0 + 256), b1 = (int)(s1 + 256);
            for (int w = b0 >> 6; w <= b1 >> 6; w++) {
                uint64_t m = cm[w];
                if (w == b0 >> 6) m &= ~UINT64_C(0) << (b0 & 63);
                if (w == b1 >> 6) m &= ~UINT64_C(0) >> (63 - (b1 & 63));
                while (m) { int64_t sv = 64 * (int64_t)w - 256 + __builtin_ctzll(m); m &= m - 1; dput(base + sv * stepU); }
            }
        } else {
            for (int64_t sv = s0; sv <= s1; sv++) if (gcd64(sv, t) == 1) dput(base + sv * stepU);
        }
    }
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
    if (argc < 6) { fprintf(stderr, "usage: k3inc points.txt L HA kA kB [-t threads] [-slots log2] [-A Alo Ahi]\n"); return 1; }
    int64_t L = chk_ll("L", argv[2], 1, (1 << 20) - 1), HA = chk_ll("HA", argv[3], 1, 1000000000); int kA = (int)chk_ll("kA", argv[4], 0, 1000000), kB = (int)chk_ll("kB", argv[5], 0, 1000000); int nthreads = omp_get_max_threads(), lg = 30;
    AHI = L;
    for (int i = 6; i < argc; i++) { if (!strcmp(argv[i], "-t")) nthreads = arg_int(argc, argv, &i, 1, 4096); else if (!strcmp(argv[i], "-slots")) lg = arg_int(argc, argv, &i, 4, 40); else if (!strcmp(argv[i], "-A")) { ALO = arg_llr(argc, argv, &i, 0, L); AHI = arg_llr(argc, argv, &i, ALO, L); } else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; } }
    if (L < 1 || L >= (1 << 20)) { fprintf(stderr, "L must be in [1, 2^20)\n"); return 1; }
    if (lg < 4 || lg > 40) { fprintf(stderr, "-slots must be in [4, 40]\n"); return 1; }
    omp_set_num_threads(nthreads);
    GL = L; while (((int64_t)1 << sbits) <= 2 * L) sbits++;   /* B+L, C+L in [0, 2L] */
    packed = 3 * sbits <= 45;                                  /* >= 19 count bits */
    if (packed) { kmask = ((uint64_t)1 << (3 * sbits)) - 1; CNT1 = (uint64_t)1 << (3 * sbits); cntmax = ((uint64_t)1 << (64 - 3 * sbits)) - 1; }
    /* read points */
    /* a point record is "x0 x1 x2" or "x0 x1 x2 W" (k3pts output), every coordinate a complete int64; comment
     * ('#') and blank lines are skipped; any other record aborts instead of being dropped from the incidence counts */
    size_t cap = 1 << 20, n = 0; int64_t (*P)[3] = malloc(cap * sizeof *P); char *buf = NULL; size_t bcap = 0; ssize_t blen; long recno = 0; FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    while ((blen = getline(&buf, &bcap, f)) > 0) { recno++; if (buf[0] == '#') continue; char *save = NULL, *tok = strtok_r(buf, " \t\r\n", &save); if (!tok) continue; int64_t v[3]; int nt = 0;
        for (; tok; tok = strtok_r(NULL, " \t\r\n", &save)) { if (nt == 3) { if (!strcmp(tok, "W") && strtok_r(NULL, " \t\r\n", &save) == NULL) break; fprintf(stderr, "%s: record %ld has more than 3 coordinates (only a trailing 'W' marker is allowed)\n", argv[1], recno); return 1; } char *e = NULL; errno = 0; long long x = strtoll(tok, &e, 10); if (errno != 0 || e == tok || *e != 0) { fprintf(stderr, "%s: record %ld, coordinate %d is not an int64 integer: '%s'\n", argv[1], recno, nt + 1, tok); return 1; } v[nt++] = x; }
        if (nt != 3) { fprintf(stderr, "%s: record %ld has %d coordinates, need 3\n", argv[1], recno, nt); return 1; }
        if (n == cap) { cap *= 2; P = realloc(P, cap * sizeof *P); } P[n][0] = v[0]; P[n][1] = v[1]; P[n][2] = v[2]; n++; }
    free(buf); fclose(f);
    size_t nA = 0; for (size_t i = 0; i < n; i++) { int64_t h = llabs(P[i][0]); if (llabs(P[i][1]) > h) h = llabs(P[i][1]); if (llabs(P[i][2]) > h) h = llabs(P[i][2]); if (h <= HA) nA++; }
    fprintf(stderr, "%zu points read, %zu with height <= %ld for stage A; L=%ld, %d threads, 2^%d slots (%s layout)%s", n, nA, HA, L, nthreads, lg, packed ? "packed" : "two-array", ALO == 0 && AHI == L ? "\n" : "");
    if (!(ALO == 0 && AHI == L)) fprintf(stderr, ", slab %ld <= A <= %ld\n", ALO, AHI);
    for (int a = 0; a < 256; a++) for (int b = 0; b < 256; b++) COP[a][b] = gcd64(a, b) == 1;
    DW = 2 * L + 1; ALO_ = ALO; DOFF = L * DW + L;
    for (int t = 1; t < 256; t++) for (int sv = -256; sv <= 255; sv++) if (gcd64(sv, t) == 1) CM[t][(sv + 256) >> 6] |= UINT64_C(1) << ((sv + 256) & 63);
    int whole = ALO == 0 && AHI == L;
    size_t dsize = (size_t)(AHI - ALO + 1) * (size_t)DW * (size_t)DW;
    if (nA == n && kA <= 254 && kB <= 254 && dsize <= ((size_t)1 << 30) && dsize * (size_t)nthreads <= ((size_t)1 << 33) && !getenv("K3INC_HASH")) {
        double td0 = omp_get_wtime();
        uint8_t **tabs = calloc(nthreads, sizeof *tabs); ovt_t **ovs = calloc(nthreads, sizeof *ovs); int nt = 0, bad = 0;
        #pragma omp parallel
        {
            uint8_t *t = table_alloc(dsize); ovt_t *o = calloc(1, sizeof *o);
            #pragma omp critical
            { if (!t || !o) bad = 1; else { tabs[nt] = t; ovs[nt] = o; nt++; } }
            DTAB = t; DOV = o;
            if (t && o) {
                if (whole) {
                    #pragma omp for schedule(dynamic, 64)
                    for (size_t i = 0; i < n; i++) dpoint(P[i][0], P[i][1], P[i][2], L);
                    dflush2();
                } else {
                    #pragma omp for schedule(dynamic, 64)
                    for (size_t i = 0; i < n; i++) lines_through(P[i][0], P[i][1], P[i][2], L, dadd);
                    dflush();
                }
            }
        }
        if (!bad && !doverflow) {
            /* the lines with at least max(kA, kB) points, in the order (A, B, C) */
            int kmin = kA > kB ? kA : kB; if (kmin < 1) kmin = 1;
            typedef struct { uint32_t c; int32_t A, B, C; } drec_t;
            int nchunk = 256; size_t chunk = (dsize + nchunk - 1) / nchunk; drec_t **crec = calloc(nchunk, sizeof *crec); size_t *cn = calloc(nchunk, sizeof *cn);
            #pragma omp parallel for schedule(dynamic, 1)
            for (int c = 0; c < nchunk; c++) {
                size_t h0 = (size_t)c * chunk, h1 = h0 + chunk < dsize ? h0 + chunk : dsize, cap = 1024, m = 0; drec_t *rr = malloc(cap * sizeof *rr);
                __m512i vk = _mm512_set1_epi8((char)kmin);
                for (size_t h = h0; h < h1; h++) {
                    if ((h & 63) == 0 && h + 64 <= h1) {   /* 64 entries at once: the saturating sum of the bytes is below kmin almost everywhere */
                        __m512i sum = _mm512_loadu_si512(tabs[0] + h); for (int q = 1; q < nt; q++) sum = _mm512_adds_epu8(sum, _mm512_loadu_si512(tabs[q] + h));
                        if (_mm512_cmpge_epu8_mask(sum, vk) == 0) { h += 63; continue; }
                    }
                    uint32_t v = 0; for (int q = 0; q < nt; q++) { uint8_t b = tabs[q][h]; v += b; if (b == 255) v += ov_get(ovs[q], h); }
                    if (v < (uint32_t)kmin) continue;
                    if (m == cap) { cap *= 2; rr = realloc(rr, cap * sizeof *rr); }
                    rr[m].c = v; rr[m].C = (int32_t)((int64_t)(h % (size_t)DW) - L); rr[m].B = (int32_t)((int64_t)((h / (size_t)DW) % (size_t)DW) - L); rr[m].A = (int32_t)((int64_t)(h / ((size_t)DW * (size_t)DW)) + ALO); m++;
                }
                crec[c] = rr; cn[c] = m;
            }
            size_t m = 0; for (int c = 0; c < nchunk; c++) m += cn[c];
            drec_t *rec = malloc((m + 1) * sizeof(drec_t)); size_t j = 0;
            for (int c = 0; c < nchunk; c++) { memcpy(rec + j, crec[c], cn[c] * sizeof(drec_t)); j += cn[c]; free(crec[c]); }
            /* by count, descending; within a count the order (A, B, C) is kept (counting sort) */
            uint32_t maxc = 0; for (size_t i = 0; i < m; i++) if (rec[i].c > maxc) maxc = rec[i].c;
            size_t *start = calloc((size_t)maxc + 2, sizeof *start); for (size_t i = 0; i < m; i++) start[maxc - rec[i].c + 1]++;
            for (uint32_t v = 1; v <= maxc + 1; v++) start[v] += start[v - 1];
            drec_t *srt = malloc((m + 1) * sizeof(drec_t)); for (size_t i = 0; i < m; i++) srt[start[maxc - rec[i].c]++] = rec[i];
            for (size_t i = 0; i < m; i++) printf("%u %u %d %d %d\n", srt[i].c, srt[i].c, srt[i].A, srt[i].B, srt[i].C);
            if (fflush(stdout) != 0 || ferror(stdout)) { fprintf(stderr, "k3inc: write error on stdout (%s) -- output incomplete\n", strerror(errno)); return 2; }
            fprintf(stderr, "direct table (%zu MB per thread): %.1fs, %zu lines with >= %d points printed\n", dsize >> 20, omp_get_wtime() - td0, m, kmin);
            return 0;
        }
        for (int q = 0; q < nt; q++) munmap(tabs[q], dsize);
        free(tabs);
        fprintf(stderr, "direct table abandoned (%s): using the hash table\n", bad ? "allocation failed" : "too many lines with more than 255 points");
    }
    size_t nslots = (size_t)1 << lg;
    hmask = nslots - 1; hprobe_limit = nslots; hkeys = table_alloc(nslots * sizeof(uint64_t)); hcnt = packed ? NULL : table_alloc(nslots * sizeof(uint32_t));
    if (!hkeys || (!packed && !hcnt)) { fprintf(stderr, "hash allocation failed\n"); return 1; }
    double t0 = omp_get_wtime();
    stage = 0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < n; i++) { int64_t h = llabs(P[i][0]); if (llabs(P[i][1]) > h) h = llabs(P[i][1]); if (llabs(P[i][2]) > h) h = llabs(P[i][2]); if (h <= HA) lines_through(P[i][0], P[i][1], P[i][2], L, kadd); }
        kflush();
    }
    double t1 = omp_get_wtime();
    /* compact: keep candidates with count >= kA (rebuilt into a fresh table of the same size, in slot order, so that
       deletions leave no tombstones and the slot order is reproducible), store countA, counts restart at 0 for stage B */
    size_t used = 0, ncand = 0;
    uint64_t *k2 = table_alloc(nslots * sizeof(uint64_t)); uint32_t *cntA = table_alloc(nslots * sizeof(uint32_t));
    if (!k2 || !cntA) { fprintf(stderr, "hash allocation failed\n"); return 1; }
    int saturated = 0;
    /* parallel scan of the stage-A table in chunks (candidates collected per chunk, in slot order), then sequential
       insertion in slot order (the insertion order decides the placement of colliding candidates) */
    int nchunk = 256; size_t chunk = (nslots + nchunk - 1) / nchunk;
    uint64_t **ckey = calloc(nchunk, sizeof *ckey); uint32_t **ccnt = calloc(nchunk, sizeof *ccnt); size_t *cn = calloc(nchunk, sizeof *cn);
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:used, saturated)
    for (int c = 0; c < nchunk; c++) {
        size_t h0 = (size_t)c * chunk, h1 = h0 + chunk < nslots ? h0 + chunk : nslots, cap = 1024, m = 0;
        uint64_t *kk = malloc(cap * sizeof *kk); uint32_t *cc = malloc(cap * sizeof *cc);
        for (size_t h = h0; h < h1; h++) if (hkeys[h]) {
            used++;
            uint64_t key, cnt;
            if (packed) { cnt = hkeys[h] >> (3 * sbits); if (cnt == cntmax) saturated++; int64_t A, B, C; decode(hkeys[h] & kmask, &A, &B, &C); key = KEY(A, B, C, L); }
            else { cnt = hcnt[h]; key = hkeys[h] - 1; }
            if (cnt < (uint64_t)kA) continue;
            if (m == cap) { cap *= 2; kk = realloc(kk, cap * sizeof *kk); cc = realloc(cc, cap * sizeof *cc); }
            kk[m] = key; cc[m] = (uint32_t)cnt; m++;
        }
        ckey[c] = kk; ccnt[c] = cc; cn[c] = m;
    }
    for (int c = 0; c < nchunk; c++) {
        for (size_t i = 0; i < cn[c]; i++) {
            uint64_t key = ckey[c][i], g = hash64(key) & hmask; while (k2[g]) g = (g + 1) & hmask;
            k2[g] = packed ? compact(key) : key + 1; cntA[g] = ccnt[c][i]; ncand++;
        }
        free(ckey[c]); free(ccnt[c]);
    }
    free(ckey); free(ccnt); free(cn);
    if (saturated) fprintf(stderr, "warning: %d stage-A counts saturated at %lu\n", saturated, cntmax);
    fprintf(stderr, "stage A: %.1fs, %zu distinct lines, %zu candidates with >= %d points (load %.2f)\n", t1 - t0, used, ncand, kA, (double)used / nslots);
    munmap(hkeys, nslots * sizeof(uint64_t)); hkeys = k2;
    if (!packed) memset(hcnt, 0, nslots * sizeof(uint32_t));
    stage = 1;
    if (nA == n) {
        /* every point was counted in stage A: the second pass would recount the same incidences, so the stage-A counts
           (saturated at cntmax in the packed layout, as stage B's would be) are the final counts */
        for (size_t g = 0; g < nslots; g++) if (hkeys[g]) { if (packed) hkeys[g] |= (uint64_t)cntA[g] << (3 * sbits); else hcnt[g] = cntA[g]; }
    } else {
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 64)
        for (size_t i = 0; i < n; i++) lines_through(P[i][0], P[i][1], P[i][2], L, kadd);
        kflush();
    }
    }
    double t2 = omp_get_wtime();
    /* output sorted by stage-B count, then stage-A count, then (A, B, C): a deterministic order (the slot order
       used before depended on which thread won a probing race in stage A, so ties came out in a different order on
       every run) */
    typedef struct { uint32_t cB, cA; int64_t A, B, C; } orec_t;
    size_t m = 0; for (size_t h = 0; h < nslots; h++) if (hkeys[h]) { uint64_t c = packed ? hkeys[h] >> (3 * sbits) : hcnt[h]; if (c >= (uint64_t)kB) m++; }
    orec_t *rec = malloc((m + 1) * sizeof(orec_t)); size_t j = 0;
    for (size_t h = 0; h < nslots; h++) if (hkeys[h]) {
        uint64_t c = packed ? hkeys[h] >> (3 * sbits) : hcnt[h]; if (c < (uint64_t)kB) continue;
        int64_t A, B, C;
        if (packed) decode(hkeys[h] & kmask, &A, &B, &C);
        else { uint64_t key = hkeys[h] - 1; A = (int64_t)(key >> 42) - L; B = (int64_t)((key >> 21) & 0x1FFFFF) - L; C = (int64_t)(key & 0x1FFFFF) - L; }
        rec[j].cB = (uint32_t)c; rec[j].cA = cntA[h]; rec[j].A = A; rec[j].B = B; rec[j].C = C; j++;
    }
    int cmpd(const void *a, const void *b) { const orec_t *x = a, *y = b;
        if (x->cB != y->cB) return x->cB < y->cB ? 1 : -1; if (x->cA != y->cA) return x->cA < y->cA ? 1 : -1;
        if (x->A != y->A) return x->A < y->A ? -1 : 1; if (x->B != y->B) return x->B < y->B ? -1 : 1; return (x->C > y->C) - (x->C < y->C); }
    qsort(rec, m, sizeof(orec_t), cmpd);
    for (size_t i = 0; i < m; i++) printf("%u %u %ld %ld %ld\n", rec[i].cB, rec[i].cA, rec[i].A, rec[i].B, rec[i].C);
    /* a full disk makes the redirected printf fail silently: on Sep 6, 2026 the 106A <= 10^4 run lost 36 of 48 slabs (98 % of its
       lines) this way with exit status 0.  Flush, check the stream, fail loudly (exit 2, like a full table). */
    if (fflush(stdout) != 0 || ferror(stdout)) { fprintf(stderr, "k3inc: write error on stdout (%s) -- output incomplete\n", strerror(errno)); return 2; }
    fprintf(stderr, "stage B: %.1fs, %zu lines with >= %d points printed\n", t2 - t1, m, kB);
    return 0;
}

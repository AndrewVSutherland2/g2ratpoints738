/* simd.h -- the 512-bit integer operations used by k3pts.c, k3inc.c, g2rp.c and cubic.c.
 *
 * On a processor with AVX-512 (the subsets F, BW, VL, VBMI2 and IFMA) these are the intrinsics of <immintrin.h>, and the
 * programs compile to the code that produced the stored run.  On any other little-endian 64-bit machine, or when
 * NO_AVX512 is defined, each operation that the programs use is defined below in plain C with the same result bit for
 * bit, so the programs compute the same thing, more slowly (see the README for the measured difference).  A masked load
 * reads only the elements selected by its mask, as the instruction does: the programs rely on that at the ends of
 * their tables.
 */
#ifndef SIMD_H
#define SIMD_H
#include <stdint.h>
#include <string.h>

#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512VL__) && defined(__AVX512VBMI2__) && defined(__AVX512IFMA__) && !defined(NO_AVX512)

#include <immintrin.h>
#include <x86intrin.h>
#define SIMD_NAME "AVX-512"

#else

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "the programs of src/ need a little-endian machine"
#endif
#include <time.h>
#define SIMD_NAME "plain C"
#define SIMD_F static inline __attribute__((always_inline, unused))

typedef union { uint64_t q[8]; uint32_t d[16]; uint8_t b[64]; } __attribute__((aligned(64))) __m512i;
typedef union { uint32_t d[8]; } __attribute__((aligned(32))) __m256i;
typedef struct { uint64_t q[2]; } __m128i;
typedef uint8_t __mmask8;
typedef uint16_t __mmask16;
typedef uint64_t __mmask64;
#define SIMD_M52 ((UINT64_C(1) << 52) - 1)

/* constants, loads, stores */
SIMD_F __m512i _mm512_setzero_si512(void) { __m512i r; memset(&r, 0, 64); return r; }
SIMD_F __m512i _mm512_set1_epi64(long long x) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = (uint64_t)x; return r; }
SIMD_F __m512i _mm512_set1_epi32(int x) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = (uint32_t)x; return r; }
SIMD_F __m512i _mm512_set1_epi8(char x) { __m512i r; memset(&r, (unsigned char)x, 64); return r; }
SIMD_F __m512i _mm512_loadu_si512(const void *p) { __m512i r; memcpy(&r, p, 64); return r; }
SIMD_F __m512i _mm512_load_si512(const void *p) { __m512i r; memcpy(&r, p, 64); return r; }
SIMD_F void _mm512_store_si512(void *p, __m512i a) { memcpy(p, &a, 64); }
SIMD_F __m512i _mm512_maskz_loadu_epi8(__mmask64 k, const void *p) {
    __m512i r; const uint8_t *s = (const uint8_t *)p; for (int i = 0; i < 64; i++) { if ((k >> i) & 1) r.b[i] = s[i]; else r.b[i] = 0; } return r; }
SIMD_F __m512i _mm512_maskz_loadu_epi64(__mmask8 k, const void *p) {
    __m512i r; const uint8_t *s = (const uint8_t *)p; for (int i = 0; i < 8; i++) { if ((k >> i) & 1) memcpy(&r.q[i], s + 8 * i, 8); else r.q[i] = 0; } return r; }
SIMD_F __m256i _mm256_setzero_si256(void) { __m256i r; memset(&r, 0, 32); return r; }
SIMD_F __m256i _mm256_set1_epi32(int x) { __m256i r; for (int i = 0; i < 8; i++) r.d[i] = (uint32_t)x; return r; }
SIMD_F __m128i _mm_cvtsi32_si128(int x) { __m128i r; r.q[0] = (uint32_t)x; r.q[1] = 0; return r; }

/* arithmetic on 64-bit, 32-bit and 8-bit elements */
SIMD_F __m512i _mm512_add_epi64(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = a.q[i] + b.q[i]; return r; }
SIMD_F __m512i _mm512_sub_epi64(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = a.q[i] - b.q[i]; return r; }
SIMD_F __m512i _mm512_and_si512(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = a.q[i] & b.q[i]; return r; }
SIMD_F __m512i _mm512_srli_epi64(__m512i a, unsigned n) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = n > 63 ? 0 : a.q[i] >> n; return r; }
SIMD_F __m512i _mm512_slli_epi64(__m512i a, unsigned n) { __m512i r; for (int i = 0; i < 8; i++) r.q[i] = n > 63 ? 0 : a.q[i] << n; return r; }
/* the 128-bit number b:a shifted right by c mod 64, low half */
SIMD_F __m512i _mm512_shrdv_epi64(__m512i a, __m512i b, __m512i c) {
    __m512i r; for (int i = 0; i < 8; i++) { unsigned s = (unsigned)(c.q[i] & 63); r.q[i] = s ? (a.q[i] >> s) | (b.q[i] << (64 - s)) : a.q[i]; } return r; }
/* a + the low (high) 52 bits of the 104-bit product of the low 52 bits of b and of c */
SIMD_F __m512i _mm512_madd52lo_epu64(__m512i a, __m512i b, __m512i c) {
    __m512i r; for (int i = 0; i < 8; i++) r.q[i] = a.q[i] + (((b.q[i] & SIMD_M52) * (c.q[i] & SIMD_M52)) & SIMD_M52); return r; }
SIMD_F __m512i _mm512_madd52hi_epu64(__m512i a, __m512i b, __m512i c) {
    __m512i r; for (int i = 0; i < 8; i++) r.q[i] = a.q[i] + (uint64_t)(((unsigned __int128)(b.q[i] & SIMD_M52) * (c.q[i] & SIMD_M52)) >> 52); return r; }
SIMD_F __m512i _mm512_add_epi32(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = a.d[i] + b.d[i]; return r; }
SIMD_F __m512i _mm512_sub_epi32(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = a.d[i] - b.d[i]; return r; }
SIMD_F __m512i _mm512_mullo_epi32(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = a.d[i] * b.d[i]; return r; }
SIMD_F __m512i _mm512_min_epu32(__m512i a, __m512i b) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = a.d[i] < b.d[i] ? a.d[i] : b.d[i]; return r; }
SIMD_F __m512i _mm512_sll_epi32(__m512i a, __m128i c) { __m512i r; for (int i = 0; i < 16; i++) r.d[i] = c.q[0] > 31 ? 0 : a.d[i] << c.q[0]; return r; }
SIMD_F __m512i _mm512_adds_epu8(__m512i a, __m512i b) {
    __m512i r; for (int i = 0; i < 64; i++) { uint8_t s = (uint8_t)(a.b[i] + b.b[i]); r.b[i] = s < a.b[i] ? 255 : s; } return r; }

/* operations under a mask: element i is computed where bit i of k is set and copied from src elsewhere */
SIMD_F __m512i _mm512_mask_sub_epi64(__m512i src, __mmask8 k, __m512i a, __m512i b) {
    __m512i r; for (int i = 0; i < 8; i++) r.q[i] = (k >> i) & 1 ? a.q[i] - b.q[i] : src.q[i]; return r; }
SIMD_F __m512i _mm512_mask_sub_epi32(__m512i src, __mmask16 k, __m512i a, __m512i b) {
    __m512i r; for (int i = 0; i < 16; i++) r.d[i] = (k >> i) & 1 ? a.d[i] - b.d[i] : src.d[i]; return r; }
SIMD_F __m512i _mm512_mask_add_epi32(__m512i src, __mmask16 k, __m512i a, __m512i b) {
    __m512i r; for (int i = 0; i < 16; i++) r.d[i] = (k >> i) & 1 ? a.d[i] + b.d[i] : src.d[i]; return r; }
SIMD_F __m512i _mm512_mask_mov_epi32(__m512i src, __mmask16 k, __m512i a) {
    __m512i r; for (int i = 0; i < 16; i++) r.d[i] = (k >> i) & 1 ? a.d[i] : src.d[i]; return r; }

/* comparisons: bit i of the result is the outcome for element i */
SIMD_F __mmask8 _mm512_test_epi64_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 8; i++) m |= (unsigned)((a.q[i] & b.q[i]) != 0) << i; return (__mmask8)m; }
SIMD_F __mmask8 _mm512_cmpeq_epi64_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 8; i++) m |= (unsigned)(a.q[i] == b.q[i]) << i; return (__mmask8)m; }
SIMD_F __mmask8 _mm512_cmpge_epu64_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 8; i++) m |= (unsigned)(a.q[i] >= b.q[i]) << i; return (__mmask8)m; }
SIMD_F __mmask16 _mm512_test_epi32_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 16; i++) m |= (unsigned)((a.d[i] & b.d[i]) != 0) << i; return (__mmask16)m; }
SIMD_F __mmask16 _mm512_cmpeq_epi32_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 16; i++) m |= (unsigned)(a.d[i] == b.d[i]) << i; return (__mmask16)m; }
SIMD_F __mmask16 _mm512_cmpge_epu32_mask(__m512i a, __m512i b) { unsigned m = 0; for (int i = 0; i < 16; i++) m |= (unsigned)(a.d[i] >= b.d[i]) << i; return (__mmask16)m; }
SIMD_F __mmask16 _mm512_mask_cmpeq_epi32_mask(__mmask16 k, __m512i a, __m512i b) { return (__mmask16)(k & _mm512_cmpeq_epi32_mask(a, b)); }
SIMD_F __mmask8 _mm256_test_epi32_mask(__m256i a, __m256i b) { unsigned m = 0; for (int i = 0; i < 8; i++) m |= (unsigned)((a.d[i] & b.d[i]) != 0) << i; return (__mmask8)m; }
/* on bytes the result is zero almost always where the programs use it, so that case is decided first, in a loop that the
   compiler can vectorise */
SIMD_F __mmask64 _mm512_cmpneq_epi8_mask(__m512i a, __m512i b) {
    uint64_t any = 0; for (int i = 0; i < 8; i++) any |= a.q[i] ^ b.q[i]; if (!any) return 0;
    uint64_t m = 0; for (int i = 0; i < 64; i++) m |= (uint64_t)(a.b[i] != b.b[i]) << i; return m; }
SIMD_F __mmask64 _mm512_mask_cmpneq_epi8_mask(__mmask64 k, __m512i a, __m512i b) { return k & _mm512_cmpneq_epi8_mask(a, b); }
SIMD_F __mmask64 _mm512_cmpge_epu8_mask(__m512i a, __m512i b) {
    uint8_t any = 0; for (int i = 0; i < 64; i++) any |= (uint8_t)(a.b[i] >= b.b[i]); if (!any) return 0;
    uint64_t m = 0; for (int i = 0; i < 64; i++) m |= (uint64_t)(a.b[i] >= b.b[i]) << i; return m; }

/* gathers: 32-bit elements read at base + index * scale (bytes) */
SIMD_F __m512i _mm512_i32gather_epi32(__m512i idx, const void *base, int scale) {
    __m512i r; for (int i = 0; i < 16; i++) memcpy(&r.d[i], (const char *)base + (int64_t)(int32_t)idx.d[i] * scale, 4); return r; }
SIMD_F __m256i _mm512_mask_i64gather_epi32(__m256i src, __mmask8 k, __m512i idx, const void *base, int scale) {
    __m256i r; for (int i = 0; i < 8; i++) { if ((k >> i) & 1) memcpy(&r.d[i], (const char *)base + (int64_t)idx.q[i] * scale, 4); else r.d[i] = src.d[i]; } return r; }

#define _MM_HINT_T0 3
#define _mm_prefetch(p, hint) __builtin_prefetch((p), 0, 3)

/* a clock for the profile that g2rp prints (it does not enter any result) */
SIMD_F uint64_t __rdtsc(void) {
#if defined(__x86_64__) || defined(__i386__)
    return __builtin_ia32_rdtsc();
#else
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (uint64_t)t.tv_sec * UINT64_C(1000000000) + (uint64_t)t.tv_nsec;
#endif
}

#endif
#endif

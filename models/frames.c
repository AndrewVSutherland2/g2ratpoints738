/* frames.c -- the frame search of reduce.gp (step 1) in C.
 *
 *   frames NKEEP < input > output
 *
 * input:  28 integers (the coefficients c of the raw sextic S, linesearch order), n, then n lines "A B C" (the known
 *         lines, in the order of reduce.gp's list LNEW).
 * For every triple a < b < d of lines with M = (l_a; l_b; l_d) invertible, reduce.gp computes the coefficient vector of
 * S(M^-1 y), normalised by a rational square to an integral vector with squarefree content (normvec), and its size score
 * (scorev: largest log10 |coefficient| + the sum of the log10 |coefficient| over 28 x 10); it keeps the NKEEP frames of
 * smallest score.  Here: S(A y) with A = adj(M) (it differs from S(M^-1 y) by det^6, a square) by evaluation at the 28
 * points (i, j, k), i + j + k = 6, and interpolation; division by the largest square dividing the content; the score in
 * double precision.
 * output: the number of invertible frames, then the frames "a b d" (1-based) with the NKEEP + 6 smallest scores (more if
 *         the next ones are within 1e-9), in increasing order of score, ties in the order of enumeration; reduce.gp
 *         recomputes these exactly and keeps NKEEP.  A frame whose content has a cofactor that cannot be decided here
 *         is printed as "F a b d" for reduce.gp to treat itself.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#define NM 28
static int32_t EI[NM], EJ[NM], EK[NM];           /* exponents of x0, x1, x2 of the n-th monomial */
static mpz_t C[NM], WN[NM][NM], WD;              /* S; inverse of the evaluation matrix = WN / WD */
static int32_t PT[NM][3];
static mpz_t L[512][3];
typedef struct { double s; int32_t a, b, d; int64_t idx; } cand_t;
static int cmpc(const void *x, const void *y) { cand_t *u = (cand_t *)x, *v = (cand_t *)y; if (u->s != v->s) return u->s < v->s ? -1 : 1; return u->idx < v->idx ? -1 : (u->idx > v->idx); }
static double mlog10(mpz_t z) { long e; double d = mpz_get_d_2exp(&e, z); return log10(fabs(d)) + (double)e * 0.30102999566398119521; }
int main(int argc, char **argv) {
  int32_t nkeep = argc > 1 ? atoi(argv[1]) : 6, n = 0, t = 0;
  for (int32_t i = 6; i >= 0; i--) for (int32_t j = 6 - i; j >= 0; j--) { EI[t] = i; EJ[t] = j; EK[t] = 6 - i - j; t++; }
  t = 0; for (int32_t i = 0; i <= 6; i++) for (int32_t j = 0; j <= 6 - i; j++) { PT[t][0] = i; PT[t][1] = j; PT[t][2] = 6 - i - j; t++; }
  for (int32_t i = 0; i < NM; i++) { mpz_init(C[i]); if (gmp_scanf("%Zd", C[i]) != 1) { fprintf(stderr, "frames: coefficients\n"); return 1; } }
  if (scanf("%d", &n) != 1 || n < 3 || n > 512) { fprintf(stderr, "frames: number of lines\n"); return 1; }
  for (int32_t i = 0; i < n; i++) for (int32_t k = 0; k < 3; k++) { mpz_init(L[i][k]); if (gmp_scanf("%Zd", L[i][k]) != 1) { fprintf(stderr, "frames: lines\n"); return 1; } }
  /* inverse of E[t][m] = monomial m at the point t, by Gauss-Jordan over Q */
  { static mpq_t E[NM][2 * NM]; mpq_t f, g; mpq_init(f); mpq_init(g);
    for (int32_t r = 0; r < NM; r++) for (int32_t c = 0; c < 2 * NM; c++) { mpq_init(E[r][c]);
      if (c < NM) { int64_t v = 1; for (int32_t e = 0; e < EI[c]; e++) v *= PT[r][0]; for (int32_t e = 0; e < EJ[c]; e++) v *= PT[r][1]; for (int32_t e = 0; e < EK[c]; e++) v *= PT[r][2]; mpq_set_si(E[r][c], v, 1); }
      else mpq_set_si(E[r][c], c - NM == r, 1); }
    for (int32_t c = 0; c < NM; c++) { int32_t p = c; while (p < NM && mpq_sgn(E[p][c]) == 0) p++;
      if (p == NM) { fprintf(stderr, "frames: singular evaluation matrix\n"); return 1; }
      if (p != c) for (int32_t k = 0; k < 2 * NM; k++) mpq_swap(E[p][k], E[c][k]);
      mpq_inv(f, E[c][c]); for (int32_t k = 0; k < 2 * NM; k++) mpq_mul(E[c][k], E[c][k], f);
      for (int32_t r = 0; r < NM; r++) if (r != c && mpq_sgn(E[r][c])) { mpq_set(f, E[r][c]); for (int32_t k = 0; k < 2 * NM; k++) { mpq_mul(g, f, E[c][k]); mpq_sub(E[r][k], E[r][k], g); } } }
    mpz_init_set_ui(WD, 1);
    for (int32_t r = 0; r < NM; r++) for (int32_t c = 0; c < NM; c++) mpz_lcm(WD, WD, mpq_denref(E[r][NM + c]));
    for (int32_t r = 0; r < NM; r++) for (int32_t c = 0; c < NM; c++) { mpz_init(WN[r][c]); mpz_divexact(WN[r][c], WD, mpq_denref(E[r][NM + c])); mpz_mul(WN[r][c], WN[r][c], mpq_numref(E[r][NM + c])); } }
  /* small primes */
  static uint32_t pr[7000]; int32_t np = 0; { static uint8_t sv[65536]; for (uint32_t p = 2; p < 65536; p++) if (!sv[p]) { pr[np++] = p; for (uint32_t q = p * p; q < 65536; q += p) sv[q] = 1; } }
  int64_t cap = (int64_t)n * (n - 1) * (n - 2) / 6, nc = 0, cnt = 0, idx = 0; cand_t *cd = malloc(cap * sizeof(cand_t));
  mpz_t A[3][3], det, x[3], pw[3][7], val[NM], v[NM], g, q, r, tt, in; mpz_inits(det, g, q, r, tt, in, NULL);
  for (int32_t i = 0; i < 3; i++) { mpz_init(x[i]); for (int32_t k = 0; k < 3; k++) mpz_init(A[i][k]); for (int32_t k = 0; k < 7; k++) mpz_init(pw[i][k]); }
  for (int32_t i = 0; i < NM; i++) { mpz_init(val[i]); mpz_init(v[i]); }
  for (int32_t a = 0; a < n; a++) for (int32_t b = a + 1; b < n; b++) for (int32_t d = b + 1; d < n; d++, idx++) {
    mpz_t *M[3] = { L[a], L[b], L[d] };
    /* adjugate: A[c][r] = cofactor (r, c) */
    for (int32_t rr = 0; rr < 3; rr++) for (int32_t cc = 0; cc < 3; cc++) { int32_t r1 = (rr + 1) % 3, r2 = (rr + 2) % 3, c1 = (cc + 1) % 3, c2 = (cc + 2) % 3;
      mpz_mul(A[cc][rr], M[r1][c1], M[r2][c2]); mpz_submul(A[cc][rr], M[r1][c2], M[r2][c1]); }
    mpz_mul(det, M[0][0], A[0][0]); mpz_addmul(det, M[0][1], A[1][0]); mpz_addmul(det, M[0][2], A[2][0]);
    if (mpz_sgn(det) == 0) continue;
    cnt++;
    for (int32_t p = 0; p < NM; p++) {
      for (int32_t i = 0; i < 3; i++) { mpz_mul_si(x[i], A[i][0], PT[p][0]); mpz_addmul_ui(x[i], A[i][1], PT[p][1]); mpz_addmul_ui(x[i], A[i][2], PT[p][2]);
        mpz_set_ui(pw[i][0], 1); for (int32_t e = 1; e <= 6; e++) mpz_mul(pw[i][e], pw[i][e - 1], x[i]); }
      mpz_set_ui(val[p], 0);
      int32_t m = 0;
      for (int32_t i = 6; i >= 0; i--) { mpz_set_ui(in, 0);
        for (int32_t j = 6 - i; j >= 0; j--, m++) { if (mpz_sgn(C[m]) == 0) continue; mpz_mul(tt, pw[1][j], pw[2][6 - i - j]); mpz_addmul(in, tt, C[m]); }
        mpz_addmul(val[p], in, pw[0][i]); } }
    mpz_set_ui(g, 0);
    for (int32_t m = 0; m < NM; m++) { mpz_set_ui(v[m], 0); for (int32_t p = 0; p < NM; p++) if (mpz_sgn(WN[m][p])) mpz_addmul(v[m], WN[m][p], val[p]);
      if (!mpz_divisible_p(v[m], WD)) { fprintf(stderr, "frames: interpolation not integral\n"); return 1; }
      mpz_divexact(v[m], v[m], WD); mpz_gcd(g, g, v[m]); }
    if (mpz_sgn(g) == 0) { fprintf(stderr, "frames: zero form\n"); return 1; }
    /* q = the largest integer with q^2 | g */
    mpz_set_ui(q, 1); mpz_set(r, g); int hard = 0;
    for (int32_t k = 0; k < np && mpz_cmp_ui(r, 1) > 0; k++) { if (mpz_cmp_ui(r, (unsigned long)pr[k] * pr[k]) < 0) break; int32_t e = 0; while (mpz_divisible_ui_p(r, pr[k])) { mpz_divexact_ui(r, r, pr[k]); e++; } for (int32_t h = 0; h < e / 2; h++) mpz_mul_ui(q, q, pr[k]); }
    if (mpz_cmp_ui(r, 4294836225UL) > 0) {      /* r > 65535^2 has no prime factor below 65536; a smaller r is 1 or a prime */
      if (mpz_perfect_square_p(r)) { mpz_sqrt(tt, r); mpz_mul(q, q, tt); }
      else if (!mpz_probab_prime_p(r, 25)) hard = 1; }
    if (hard) { printf("F %d %d %d\n", a + 1, b + 1, d + 1); continue; }
    double lq = 2.0 * mlog10(q), mx = -1e300, sm = 0;
    for (int32_t m = 0; m < NM; m++) if (mpz_sgn(v[m])) { double l = mlog10(v[m]) - lq; if (l > mx) mx = l; sm += l; }
    cd[nc].s = mx + sm / 28.0 / 10.0; cd[nc].a = a + 1; cd[nc].b = b + 1; cd[nc].d = d + 1; cd[nc].idx = idx; nc++;
  }
  qsort(cd, nc, sizeof(cand_t), cmpc);
  int64_t K = nkeep + 6; if (K > nc) K = nc; while (K < nc && K < nkeep + 60 && cd[K].s - cd[K - 1].s < 1e-9) K++;
  printf("N %ld\n", (long)cnt);
  for (int64_t k = 0; k < K; k++) printf("C %d %d %d %.12f\n", cd[k].a, cd[k].b, cd[k].d, cd[k].s);
  return 0;
}

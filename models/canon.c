/* canon.c -- canonical form of the configuration of the classes of degree 0 and 1 of a polarization D' of X_163.
 *
 *   canon LAT < rows of config.gp > "id n2 index nroots nlines hash"
 *
 * A row gives the contracted curves and the line classes of D' (see config.gp) as vectors of NS.  The
 * configuration is the matrix of their intersection numbers, up to a simultaneous permutation of rows and columns that
 * keeps roots and classes of degree 1 apart.  Canonical form: colour refinement (a class is recoloured by its colour
 * and the multiset of (colour, intersection number) with all classes, through a 64-bit hash of that multiset, which is
 * an invariant whatever its collisions) until the partition is stable; while some colour class has more than one member,
 * every member of the first such class is singled out in turn and the search continues; the smallest matrix reached is
 * the canonical form.  It is a matrix of the configuration, so equal canonical forms mean isomorphic configurations
 * (when the classes generate NS, index 1, this is an isometry of NS taking one polarization to the other, hence an
 * isomorphism of the double planes), and isomorphic configurations have equal canonical forms because every step is
 * invariant.  The output is a 128-bit hash of the canonical matrix.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXM 600
static int64_t GZ[20][20];
static int32_t m, nr;
static uint8_t W[MAXM][MAXM];
static uint64_t T[2 * MAXM * 16 + 16];
static uint8_t *cur, *best; static int32_t have;
static uint64_t leaves;
static uint64_t smix(uint64_t *s) { uint64_t z = (*s += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
typedef struct { int32_t c; uint64_t h; int32_t i; } ckey_t;
static int cmpkey(const void *a, const void *b) { ckey_t *x = (ckey_t *)a, *y = (ckey_t *)b; if (x->c != y->c) return x->c < y->c ? -1 : 1; if (x->h != y->h) return x->h < y->h ? -1 : 1; return 0; }
/* refine: col (any non-negative values) is replaced by the stable refinement, numbered 0..nc-1 invariantly */
static int32_t refine(int32_t *col) {
  static ckey_t k[MAXM]; int32_t nc = -1, n;
  for (int32_t i = 0; i < m; i++) { k[i].c = col[i]; k[i].h = 0; k[i].i = i; }
  for (;;) {
    qsort(k, m, sizeof(ckey_t), cmpkey);
    n = 0; for (int32_t a = 0; a < m; a++) { if (a && (k[a].c != k[a-1].c || k[a].h != k[a-1].h)) n++; col[k[a].i] = n; } n++;
    if (n == nc || n == m) return n;
    nc = n;
    for (int32_t a = 0; a < m; a++) { int32_t i = k[a].i; uint64_t s = 0; uint8_t *w = W[i];
      for (int32_t j = 0; j < m; j++) s += T[col[j] * 16 + w[j]];
      k[a].c = col[i]; k[a].h = s; }
  }
}
static void canon(int32_t *col0) {
  int32_t col[MAXM], c2[MAXM], cnt[MAXM], p[MAXM]; memcpy(col, col0, m * sizeof(int32_t));
  int32_t nc = refine(col);
  if (nc == m) { leaves++;
    for (int32_t i = 0; i < m; i++) p[col[i]] = i;
    for (int32_t a = 0; a < m; a++) for (int32_t b = 0; b < m; b++) cur[a * m + b] = W[p[a]][p[b]];
    if (!have || memcmp(cur, best, (size_t)m * m) < 0) { memcpy(best, cur, (size_t)m * m); have = 1; }
    return; }
  memset(cnt, 0, nc * sizeof(int32_t)); for (int32_t i = 0; i < m; i++) cnt[col[i]]++;
  int32_t c = 0; while (cnt[c] == 1) c++;
  for (int32_t v = 0; v < m; v++) if (col[v] == c) {
    for (int32_t i = 0; i < m; i++) c2[i] = 2 * col[i] + 1;
    c2[v]--; canon(c2); }
}
int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: canon LAT < rows\n"); return 2; }
  FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
  size_t cap = 1 << 24; char *line = malloc(cap); int ok = 0;
  while (getline(&line, &cap, f) > 0) if (!strncmp(line, "GZ ", 3)) { char *q = line + 3; for (int a = 0; a < 20; a++) for (int b = 0; b < 20; b++) GZ[a][b] = strtoll(q, &q, 10); ok = 1; break; }
  fclose(f); if (!ok) { fprintf(stderr, "no GZ row in %s\n", argv[1]); return 1; }
  uint64_t seed = 163; for (size_t i = 0; i < sizeof(T) / sizeof(T[0]); i++) T[i] = smix(&seed);
  cur = malloc(MAXM * MAXM); best = malloc(MAXM * MAXM);
  static int64_t V[MAXM][20], GV[MAXM][20];
  while (getline(&line, &cap, stdin) > 0) {
    char id[128]; int n2, nl; long idx; int pos;
    if (sscanf(line, "%127s %d %ld %d %d | %n", id, &n2, &idx, &nr, &nl, &pos) < 5) continue;
    m = nr + nl; if (m > MAXM) { fprintf(stderr, "%s: %d classes, more than MAXM\n", id, m); return 1; }
    char *q = line + pos;
    for (int32_t i = 0; i < m; i++) for (int a = 0; a < 20; a++) V[i][a] = strtoll(q, &q, 10);
    for (int32_t i = 0; i < m; i++) for (int a = 0; a < 20; a++) { int64_t s = 0; for (int b = 0; b < 20; b++) s += GZ[a][b] * V[i][b]; GV[i][a] = s; }
    for (int32_t i = 0; i < m; i++) for (int32_t j = 0; j < m; j++) { int64_t s = 0; for (int a = 0; a < 20; a++) s += V[i][a] * GV[j][a]; s += 4;
      if (s < 0 || s > 15) { fprintf(stderr, "%s: intersection number %ld out of range\n", id, (long)(s - 4)); return 1; } W[i][j] = (uint8_t)s; }
    int32_t col[MAXM]; for (int32_t i = 0; i < m; i++) col[i] = i < nr ? 0 : 1;
    have = 0; leaves = 0; canon(col);
    uint64_t h1 = 0xcbf29ce484222325ULL ^ (uint64_t)m, h2 = 0x84222325cbf29ce4ULL ^ (uint64_t)nr;
    for (size_t t = 0; t < (size_t)m * m; t++) { h1 = (h1 ^ best[t]) * 0x100000001b3ULL; h2 = (h2 + best[t] + 1) * 0x9e3779b97f4a7c15ULL; h2 ^= h2 >> 29; }
    printf("%s %d %ld %d %d %016lx%016lx %lu\n", id, n2, idx, nr, nl, (unsigned long)h1, (unsigned long)h2, (unsigned long)leaves);
  }
  return 0;
}

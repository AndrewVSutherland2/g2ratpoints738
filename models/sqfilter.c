/* sqfilter.c : reads lines "A B C f0 f1 ... f6" (output of src/linesearch <sextic> 1 -exactbatch: the binary sextic
 * f(x) = sum f_d x^d, x = t/s, of the restriction of the sextic to the line A x0 + B x1 + C x2 = 0) and prints "A B C"
 * for the lines on which f is the square of a binary cubic form (exact test with GMP: a square over Q is the square of an
 * integral form since f is integral).  Build: gcc -O2 -o sqfilter sqfilter.c -lgmp */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

/* is the binary form F(x, 1) = sum_{d=lo}^{hi} f[d] x^d (f[lo], f[hi] != 0, degree 6 as a binary form) a perfect square? */
static int is_square_form(mpz_t f[7]) {
    int lo = 0, hi = 6;
    while (lo <= 6 && mpz_sgn(f[lo]) == 0) lo++;
    if (lo > 6) return 1;   /* the zero form (degenerate line): reported as a square, sorted out by the caller */
    while (mpz_sgn(f[hi]) == 0) hi--;
    if ((lo & 1) || ((6 - hi) & 1)) return 0;      /* the factors x^lo and s^(6-hi) must have even exponents */
    int m = hi - lo;                                 /* degree of the middle part, even */
    if (m & 1) return 0;
    int h = m / 2;
    if (!mpz_perfect_square_p(f[hi])) return 0;
    mpz_t g[4], t1, t2; for (int k = 0; k <= 3; k++) mpz_init(g[k]); mpz_init(t1); mpz_init(t2);
    mpz_sqrt(g[h], f[hi]);
    int ok = 1;
    /* g_k for k = h-1 .. 0: coefficient of x^(lo + h + k) of g^2 is 2 g_h g_k + sum_{i+j = h+k, k < i,j < h} g_i g_j */
    for (int k = h - 1; k >= 0 && ok; k--) {
        mpz_set(t1, f[lo + h + k]);
        for (int i = k + 1; i <= h - 1; i++) { int j = h + k - i; if (j <= k || j >= h || j < i) continue; if (i == j) { mpz_mul(t2, g[i], g[i]); mpz_sub(t1, t1, t2); } else { mpz_mul(t2, g[i], g[j]); mpz_submul_ui(t1, t2, 2); } }
        mpz_mul_ui(t2, g[h], 2);
        if (!mpz_divisible_p(t1, t2)) { ok = 0; break; }
        mpz_divexact(g[k], t1, t2);
    }
    if (ok) {   /* full check g^2 == middle part */
        for (int d = 0; d <= m && ok; d++) {
            mpz_set_ui(t1, 0);
            for (int i = 0; i <= h; i++) { int j = d - i; if (j < 0 || j > h) continue; mpz_mul(t2, g[i], g[j]); mpz_add(t1, t1, t2); }
            if (mpz_cmp(t1, f[lo + d]) != 0) ok = 0;
        }
    }
    for (int k = 0; k <= 3; k++) mpz_clear(g[k]);
    mpz_clear(t1); mpz_clear(t2);
    return ok;
}

int main(void) {
    char *buf = NULL; size_t cap = 0; ssize_t len;
    mpz_t f[7]; for (int d = 0; d < 7; d++) mpz_init(f[d]);
    long nread = 0, nsq = 0;
    while ((len = getline(&buf, &cap, stdin)) > 0) {
        char *save = NULL, *tok = strtok_r(buf, " \t\r\n", &save);
        char *abc[3]; int n = 0, bad = 0;
        for (n = 0; n < 3 && tok; n++) { abc[n] = tok; tok = strtok_r(NULL, " \t\r\n", &save); }
        if (n < 3) { fprintf(stderr, "sqfilter: short record: %s\n", buf); return 1; }
        for (int d = 0; d < 7; d++) { if (!tok || mpz_set_str(f[d], tok, 10) != 0) { bad = 1; break; } tok = strtok_r(NULL, " \t\r\n", &save); }
        if (bad || tok) { fprintf(stderr, "sqfilter: malformed record (expected A B C and 7 integers)\n"); return 1; }
        nread++;
        if (is_square_form(f)) { nsq++; printf("%s %s %s\n", abc[0], abc[1], abc[2]); }
    }
    fprintf(stderr, "sqfilter: %ld lines read, %ld squares\n", nread, nsq);
    for (int d = 0; d < 7; d++) mpz_clear(f[d]);
    free(buf);
    return 0;
}

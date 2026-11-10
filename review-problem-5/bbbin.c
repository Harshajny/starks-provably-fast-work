// BabyBear (p = 2013265921), mu32 = <zeta>, zeta = 760005850; RS[mu32, deg <= 15] (k = 16), agreement 18.
// A 14-set C (complement of an 18-set S) passes for a word with syndrome V (V_m = 32 a_{32-m}, m = 1..16)
// iff l_e = sum_{j=0..14} c_j(C) V_{e+1+j} = 0 for e = 0, 1 (GRS dual; prod_{t != s}(s - t) = 32/s on mu32).
// Binomial V_m1 = 1, V_m2 = lam: rows (c_{m1-1-e}, c_{m2-1-e}), e = 0,1. Rank 1: one lam; rank 0: every lam.
// Monomial V_m: c_{m-1} = c_{m-2} = 0 (index 14, where c_14 = 1, never vanishes).
// Output: "B m1 m2 lam mask" | "U m1 m2 0 mask" | "M m m 0 mask".
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef uint64_t u64;
#define P 2013265921ULL
static inline u64 mm(u64 a, u64 b) { return a * b % P; }
static inline u64 ms(u64 a, u64 b) { return (a + P - b) % P; }
static u64 pw(u64 b, u64 e) { u64 r = 1; while (e) { if (e & 1) r = mm(r, b); b = mm(b, b); e >>= 1; } return r; }
static u64 Z[32]; static int tj0[600], tj1[600], ntask, nt; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
typedef struct { char *b; size_t l, c; } buf_t;
static void em(buf_t *B, const char *s, int n) { if (B->l + n > B->c) { B->c = (B->c + n) * 2; B->b = realloc(B->b, B->c); } memcpy(B->b + B->l, s, n); B->l += n; }
static inline u64 C(const u64 *c, int j) { return (j < 0 || j > 14) ? 0 : c[j]; }
static void leaf(const u64 *c, uint32_t mask, buf_t *B) {
    char s[80]; int n;
    for (int m = 1; m <= 16; m++) if (m - 1 != 14 && m - 2 != 14 && !C(c, m - 1) && !C(c, m - 2)) { n = snprintf(s, 80, "M %d %d 0 %08x\n", m, m, mask); em(B, s, n); }
    for (int a = 1; a <= 16; a++) for (int b = a + 1; b <= 16; b++) {
        u64 x0 = C(c, a - 1), x1 = C(c, a - 2), y0 = C(c, b - 1), y1 = C(c, b - 2);
        if (ms(mm(x0, y1), mm(x1, y0))) continue;
        int k = (x0 || y0) ? 0 : ((x1 || y1) ? 1 : -1);
        if (k < 0) { n = snprintf(s, 80, "U %d %d 0 %08x\n", a, b, mask); em(B, s, n); continue; }
        u64 x = k ? x1 : x0, y = k ? y1 : y0;
        if (!x || !y) continue;
        u64 lam = mm(P - x, pw(y, P - 2));
        n = snprintf(s, 80, "B %d %d %llu %08x\n", a, b, (unsigned long long)lam, mask); em(B, s, n);
    }
}
static void rec(u64 c[15][15], int d, int st, uint32_t mask, buf_t *B) {
    if (d == 14) { leaf(c[14], mask, B); return; }
    for (int j = st; j <= 32 - (14 - d); j++) {
        for (int t = 0; t <= d + 1; t++) { u64 hi = t ? c[d][t - 1] : 0, lo = t <= d ? mm(c[d][t], Z[j]) : 0; c[d + 1][t] = ms(hi, lo); }
        rec(c, d + 1, j + 1, mask | (1u << j), B);
    }
}
static void *wk(void *a) {
    (void)a; u64 c[15][15]; buf_t B = { 0 };
    for (;;) { pthread_mutex_lock(&mu); int k = nt++; pthread_mutex_unlock(&mu); if (k >= ntask) break;
        int js[2] = { tj0[k], tj1[k] }; c[0][0] = 1;
        for (int s = 0; s < 2; s++) for (int t = 0; t <= s + 1; t++) { u64 hi = t ? c[s][t - 1] : 0, lo = t <= s ? mm(c[s][t], Z[js[s]]) : 0; c[s + 1][t] = ms(hi, lo); }
        rec(c, 2, js[1] + 1, (1u << js[0]) | (1u << js[1]), &B);
        if (B.l > (1u << 20)) { pthread_mutex_lock(&mu); fwrite(B.b, 1, B.l, stdout); pthread_mutex_unlock(&mu); B.l = 0; } }
    pthread_mutex_lock(&mu); fwrite(B.b, 1, B.l, stdout); pthread_mutex_unlock(&mu); free(B.b); return NULL;
}
int main(int argc, char **argv) {
    int th = argc > 1 ? atoi(argv[1]) : 10; u64 z = 760005850ULL;
    if (pw(z, 16) != P - 1 || pw(z, 32) != 1) { fprintf(stderr, "zeta is not of order 32\n"); return 1; }
    Z[0] = 1; for (int j = 1; j < 32; j++) Z[j] = mm(Z[j - 1], z);
    for (int a = 0; a < 32; a++) for (int b = a + 1; b < 32; b++) if (32 - b - 1 >= 12) { tj0[ntask] = a; tj1[ntask] = b; ntask++; }
    pthread_t t[64]; for (int k = 0; k < th; k++) pthread_create(&t[k], 0, wk, 0); for (int k = 0; k < th; k++) pthread_join(t[k], 0);
    return 0;
}

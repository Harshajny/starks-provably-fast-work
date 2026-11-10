// Exhaustive list sizes for every word with at most two high monomials, RS[mu32, 17], agreement 20,
// over F_{p^2} = F_p[i], p = 2^31 - 1, mu32 = <w>, w = 906276279 + 1179735656 i.
//
// A word mod the code is v = sum_{d=17..31} a_d t^d; its syndrome is V_m = 32 a_{32-m}, m = 1..15.
// A 12-set C (complement of a 20-set S) passes iff l_e = sum_j c_j(C) V_{e+1+j} = 0 for e = 0,1,2,
// where prod_{t in C}(x - t) = sum_{j=0}^{12} c_j x^j (c_j = 0 outside 0..12, c_12 = 1).
// Binomial V (V_m1 = alpha, V_m2 = beta, m1 < m2): rows (x_e, y_e) = (c_{m1-1-e}, c_{m2-1-e});
// C passes for (alpha:beta) iff alpha x_e + beta y_e = 0 for all e. Rank 2: never. Rank 1: one ratio
// lambda = beta/alpha = -x/y (a true binomial needs x, y != 0). Rank 0: every lambda ("universal").
// Monomial V_m: C passes iff c_{m-1-e} = 0 for e = 0,1,2.
//
// Output, one line per passing (class, C):  kind m1 m2 lam_re lam_im mask
//   kind B = binomial at lambda, U = binomial universal (every lambda), M = monomial (m1 = m2 = m).
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef uint64_t u64;
#define P 2147483647ULL
typedef struct { u64 r, i; } fe;

static inline u64 red(u64 x) { x = (x & P) + (x >> 31); x = (x & P) + (x >> 31); return x >= P ? x - P : x; }
static inline fe fmul(fe a, fe b) {
    u64 ac = red(a.r * b.r), bd = red(a.i * b.i);
    fe z = { red(ac + P - bd), red(a.r * b.i + a.i * b.r) };
    return z;
}
static inline fe fsub(fe a, fe b) { fe z = { red(a.r + P - b.r), red(a.i + P - b.i) }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }

static fe Z[32];
static int task_j0[600], task_j1[600], ntask;
static int next_task = 0;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static FILE *out;
static u64 leaves_total = 0;

typedef struct { char *buf; size_t len, cap; u64 leaves; } tctx;
static void emit(tctx *t, const char *s, int n) {
    if (t->len + n + 1 > t->cap) { t->cap = (t->cap + n + 1) * 2; t->buf = realloc(t->buf, t->cap); }
    for (int k = 0; k < n; k++) t->buf[t->len++] = s[k];
}

static inline fe C(const fe *c, int j) { fe z = { 0, 0 }; if (j < 0 || j > 12) return z; return c[j]; }

static void leaf(tctx *t, const fe *c, uint32_t mask) {
    t->leaves++;
    char s[96];
    // monomials
    for (int m = 1; m <= 15; m++) {
        if (fz(C(c, m - 1)) && fz(C(c, m - 2)) && fz(C(c, m - 3)) && m - 1 != 12 && m - 2 != 12 && m - 3 != 12) {
            int n = snprintf(s, sizeof s, "M %d %d 0 0 %08x\n", m, m, mask); emit(t, s, n);
        }
    }
    // binomials
    for (int m1 = 1; m1 <= 15; m1++) for (int m2 = m1 + 1; m2 <= 15; m2++) {
        fe x0 = C(c, m1 - 1), x1 = C(c, m1 - 2), x2 = C(c, m1 - 3);
        fe y0 = C(c, m2 - 1), y1 = C(c, m2 - 2), y2 = C(c, m2 - 3);
        if (!fz(fsub(fmul(x0, y1), fmul(x1, y0)))) continue;
        if (!fz(fsub(fmul(x0, y2), fmul(x2, y0)))) continue;
        if (!fz(fsub(fmul(x1, y2), fmul(x2, y1)))) continue;
        fe xs[3] = { x0, x1, x2 }, ys[3] = { y0, y1, y2 };
        int k = -1;
        for (int e = 0; e < 3; e++) if (!fz(xs[e]) || !fz(ys[e])) { k = e; break; }
        int n;
        if (k < 0) { n = snprintf(s, sizeof s, "U %d %d 0 0 %08x\n", m1, m2, mask); emit(t, s, n); continue; }
        if (fz(xs[k]) || fz(ys[k])) continue; // kernel is a monomial, counted under M
        fe lam = fmul(xs[k], finv(ys[k])); lam.r = (P - lam.r) % P; lam.i = (P - lam.i) % P; // -x/y
        n = snprintf(s, sizeof s, "B %d %d %llu %llu %08x\n", m1, m2, (unsigned long long)lam.r, (unsigned long long)lam.i, mask);
        emit(t, s, n);
    }
}

static void rec(tctx *t, fe c[13][13], int depth, int start, uint32_t mask) {
    if (depth == 12) { leaf(t, c[12], mask); return; }
    for (int j = start; j <= 32 - (12 - depth); j++) {
        fe *o = c[depth], *nw = c[depth + 1];
        // nw = o * (x - Z[j]); o has degree depth
        fe zero = { 0, 0 };
        for (int d = 0; d <= depth + 1; d++) {
            fe hi = d >= 1 ? o[d - 1] : zero;
            fe lo = d <= depth ? fmul(o[d], Z[j]) : zero;
            nw[d] = fsub(hi, lo);
        }
        rec(t, c, depth + 1, j + 1, mask | (1u << j));
    }
}

static void *worker(void *arg) {
    (void)arg;
    tctx t = { 0 };
    fe c[13][13];
    for (;;) {
        pthread_mutex_lock(&mu); int k = next_task++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int j0 = task_j0[k], j1 = task_j1[k];
        c[0][0].r = 1; c[0][0].i = 0;
        uint32_t mask = 0;
        int js[2] = { j0, j1 };
        for (int s = 0; s < 2; s++) {
            fe *o = c[s], *nw = c[s + 1]; fe zero = { 0, 0 };
            for (int d = 0; d <= s + 1; d++) {
                fe hi = d >= 1 ? o[d - 1] : zero;
                fe lo = d <= s ? fmul(o[d], Z[js[s]]) : zero;
                nw[d] = fsub(hi, lo);
            }
            mask |= 1u << js[s];
        }
        rec(&t, c, 2, j1 + 1, mask);
        if (t.len > (1u << 20)) { pthread_mutex_lock(&mu); fwrite(t.buf, 1, t.len, out); pthread_mutex_unlock(&mu); t.len = 0; }
    }
    pthread_mutex_lock(&mu); fwrite(t.buf, 1, t.len, out); leaves_total += t.leaves; pthread_mutex_unlock(&mu);
    free(t.buf);
    return NULL;
}

int main(int argc, char **argv) {
    int nth = argc > 1 ? atoi(argv[1]) : 10;
    out = stdout;
    fe w = { 906276279ULL, 1179735656ULL };
    Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int j0 = 0; j0 < 32; j0++) for (int j1 = j0 + 1; j1 < 32; j1++)
        if (32 - j1 - 1 >= 10) { task_j0[ntask] = j0; task_j1[ntask] = j1; ntask++; }
    pthread_t th[64];
    for (int k = 0; k < nth; k++) pthread_create(&th[k], NULL, worker, NULL);
    for (int k = 0; k < nth; k++) pthread_join(th[k], NULL);
    fprintf(stderr, "leaves=%llu tasks=%d\n", (unsigned long long)leaves_total, ntask);
    return 0;
}

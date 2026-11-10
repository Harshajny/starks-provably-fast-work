// Exhaustive pass for every word with exactly three high monomials, RS[mu32, 17], agreement 20,
// over F_{p^2} = F_p[i], p = 2^31 - 1, mu32 = <w>, w = 906276279 + 1179735656 i.
//
// Syndrome V supported on m1 < m2 < m3 (values alpha, beta, gamma). A 12-set C passes iff
// M_C (alpha, beta, gamma)^T = 0, where M_C[e] = (c_{m1-1-e}, c_{m2-1-e}, c_{m3-1-e}), e = 0, 1, 2,
// and prod_{t in C}(x - t) = sum c_j x^j (c_j = 0 outside 0..12, c_12 = 1). See NOTES.md.
//   det M_C != 0 : passes for nothing.
//   rank 2       : passes for exactly one point (1 : b : g); printed only if genuine (b, g != 0).
//   rank 1       : passes for a projective line n . (alpha, beta, gamma) = 0; printed as its normal n.
//   rank 0       : passes for every point.
// Output lines:  P m1 m2 m3 b_re b_im g_re g_im mask | L m1 m2 m3 n0r n0i n1r n1i n2r n2i mask | U m1 m2 m3 mask
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
static inline fe fadd(fe a, fe b) { fe z = { red(a.r + b.r), red(a.i + b.i) }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }

static fe Z[32];
static int task_j0[600], task_j1[600], ntask, next_task = 0;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static u64 leaves_total = 0;

typedef struct { char *buf; size_t len, cap; u64 leaves; } tctx;
static void emit(tctx *t, const char *s, int n) {
    if (t->len + n + 1 > t->cap) { t->cap = (t->cap + n + 1) * 2; t->buf = realloc(t->buf, t->cap); }
    for (int k = 0; k < n; k++) t->buf[t->len++] = s[k];
}

static inline void cross(const fe a[3], const fe b[3], fe o[3]) {
    o[0] = fsub(fmul(a[1], b[2]), fmul(a[2], b[1]));
    o[1] = fsub(fmul(a[2], b[0]), fmul(a[0], b[2]));
    o[2] = fsub(fmul(a[0], b[1]), fmul(a[1], b[0]));
}

static void leaf(tctx *t, const fe *c, uint32_t mask) {
    t->leaves++;
    fe col[16][3];
    fe zero = { 0, 0 };
    for (int m = 1; m <= 15; m++)
        for (int e = 0; e < 3; e++) { int j = m - 1 - e; col[m][e] = (j < 0 || j > 12) ? zero : c[j]; }
    char s[200];
    for (int m1 = 1; m1 <= 13; m1++) for (int m2 = m1 + 1; m2 <= 14; m2++) {
        fe w[3]; cross(col[m1], col[m2], w);
        for (int m3 = m2 + 1; m3 <= 15; m3++) {
            fe det = fadd(fadd(fmul(w[0], col[m3][0]), fmul(w[1], col[m3][1])), fmul(w[2], col[m3][2]));
            if (!fz(det)) continue;
            fe row[3][3];
            for (int e = 0; e < 3; e++) { row[e][0] = col[m1][e]; row[e][1] = col[m2][e]; row[e][2] = col[m3][e]; }
            fe k[3]; int found = 0;
            int pr[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
            for (int q = 0; q < 3 && !found; q++) {
                cross(row[pr[q][0]], row[pr[q][1]], k);
                if (!fz(k[0]) || !fz(k[1]) || !fz(k[2])) found = 1;
            }
            int n;
            if (found) { // rank 2: unique point
                if (fz(k[0]) || fz(k[1]) || fz(k[2])) continue; // not a genuine trinomial
                fe inv0 = finv(k[0]); fe b = fmul(k[1], inv0), g = fmul(k[2], inv0);
                n = snprintf(s, sizeof s, "P %d %d %d %llu %llu %llu %llu %08x\n", m1, m2, m3,
                             (unsigned long long)b.r, (unsigned long long)b.i, (unsigned long long)g.r, (unsigned long long)g.i, mask);
                emit(t, s, n); continue;
            }
            int r = -1; // rank <= 1: find a nonzero row
            for (int e = 0; e < 3 && r < 0; e++) if (!fz(row[e][0]) || !fz(row[e][1]) || !fz(row[e][2])) r = e;
            if (r < 0) { n = snprintf(s, sizeof s, "U %d %d %d %08x\n", m1, m2, m3, mask); emit(t, s, n); continue; }
            int f = 0; while (fz(row[r][f])) f++;
            fe inv = finv(row[r][f]); fe nn[3];
            for (int q = 0; q < 3; q++) nn[q] = fmul(row[r][q], inv);
            n = snprintf(s, sizeof s, "L %d %d %d %llu %llu %llu %llu %llu %llu %08x\n", m1, m2, m3,
                         (unsigned long long)nn[0].r, (unsigned long long)nn[0].i, (unsigned long long)nn[1].r,
                         (unsigned long long)nn[1].i, (unsigned long long)nn[2].r, (unsigned long long)nn[2].i, mask);
            emit(t, s, n);
        }
    }
}

static void rec(tctx *t, fe c[13][13], int depth, int start, uint32_t mask) {
    if (depth == 12) { leaf(t, c[12], mask); return; }
    for (int j = start; j <= 32 - (12 - depth); j++) {
        fe *o = c[depth], *nw = c[depth + 1]; fe zero = { 0, 0 };
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
    tctx t = { 0 }; fe c[13][13];
    for (;;) {
        pthread_mutex_lock(&mu); int k = next_task++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int js[2] = { task_j0[k], task_j1[k] };
        c[0][0].r = 1; c[0][0].i = 0; uint32_t mask = 0;
        for (int s = 0; s < 2; s++) {
            fe *o = c[s], *nw = c[s + 1]; fe zero = { 0, 0 };
            for (int d = 0; d <= s + 1; d++) {
                fe hi = d >= 1 ? o[d - 1] : zero;
                fe lo = d <= s ? fmul(o[d], Z[js[s]]) : zero;
                nw[d] = fsub(hi, lo);
            }
            mask |= 1u << js[s];
        }
        rec(&t, c, 2, js[1] + 1, mask);
        if (t.len > (1u << 20)) { pthread_mutex_lock(&mu); fwrite(t.buf, 1, t.len, stdout); pthread_mutex_unlock(&mu); t.len = 0; }
    }
    pthread_mutex_lock(&mu); fwrite(t.buf, 1, t.len, stdout); leaves_total += t.leaves; pthread_mutex_unlock(&mu);
    free(t.buf);
    return NULL;
}

int main(int argc, char **argv) {
    int nth = argc > 1 ? atoi(argv[1]) : 10;
    int limit = argc > 2 ? atoi(argv[2]) : 1 << 30; // optional: only the first `limit` tasks (timing test)
    fe w = { 906276279ULL, 1179735656ULL };
    Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int j0 = 0; j0 < 32; j0++) for (int j1 = j0 + 1; j1 < 32; j1++)
        if (32 - j1 - 1 >= 10 && ntask < limit) { task_j0[ntask] = j0; task_j1[ntask] = j1; ntask++; }
    pthread_t th[64];
    for (int k = 0; k < nth; k++) pthread_create(&th[k], NULL, worker, NULL);
    for (int k = 0; k < nth; k++) pthread_join(th[k], NULL);
    fprintf(stderr, "leaves=%llu tasks=%d\n", (unsigned long long)leaves_total, ntask);
    return 0;
}

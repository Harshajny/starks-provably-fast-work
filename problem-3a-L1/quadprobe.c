// Feasibility probe for one four-monomial support (m1..m4): for each 12-set C, the 3x4 system
// M_C (a1..a4)^T = 0 (rows e = 0,1,2: c_{mq-1-e}). Counts C by rank of M_C, and counts how many
// rank-3 C (one word each) are canonical for their member (the 21-set test of trinom2.c).
// Usage: ./quadprobe THREADS TASK_LIMIT m1 m2 m3 m4
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef uint64_t u64;
#define P 2147483647ULL
typedef struct { u64 r, i; } fe;
static inline u64 red(u64 x) { x = (x & P) + (x >> 31); x = (x & P) + (x >> 31); return x >= P ? x - P : x; }
static inline fe fmul(fe a, fe b) { u64 ac = red(a.r * b.r), bd = red(a.i * b.i); fe z = { red(ac + P - bd), red(a.r * b.i + a.i * b.r) }; return z; }
static inline fe fsub(fe a, fe b) { fe z = { red(a.r + P - b.r), red(a.i + P - b.i) }; return z; }
static inline fe fadd(fe a, fe b) { fe z = { red(a.r + b.r), red(a.i + b.i) }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static fe Z[32]; static int M[4];
static int tj0[600], tj1[600], ntask, nexttask; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static u64 g_rank[4], g_canon;

static inline fe det3(fe a[3], fe b[3], fe c[3]) {
    fe x = fsub(fmul(b[1], c[2]), fmul(b[2], c[1])), y = fsub(fmul(b[2], c[0]), fmul(b[0], c[2])), z = fsub(fmul(b[0], c[1]), fmul(b[1], c[0]));
    return fadd(fadd(fmul(a[0], x), fmul(a[1], y)), fmul(a[2], z));
}

static void leaf(const fe *c, uint32_t mask, u64 *rk, u64 *canon) {
    fe col[4][3], zero = { 0, 0 };
    for (int q = 0; q < 4; q++) for (int e = 0; e < 3; e++) { int j = M[q] - 1 - e; col[q][e] = (j < 0 || j > 12) ? zero : c[j]; }
    // kernel of 3x4 via signed 3x3 minors
    fe k[4];
    k[0] = det3(col[1], col[2], col[3]);
    k[1] = det3(col[0], col[2], col[3]); k[1].r = (P - k[1].r) % P; k[1].i = (P - k[1].i) % P;
    k[2] = det3(col[0], col[1], col[3]);
    k[3] = det3(col[0], col[1], col[2]); k[3].r = (P - k[3].r) % P; k[3].i = (P - k[3].i) % P;
    if (fz(k[0]) && fz(k[1]) && fz(k[2]) && fz(k[3])) { rk[2]++; return; } // rank <= 2 (lines/planes), counted only
    rk[3]++;
    if (fz(k[0]) || fz(k[1]) || fz(k[2]) || fz(k[3])) return; // fewer than four monomials: covered elsewhere
    int maxS = 31; while ((mask >> maxS) & 1) maxS--;
    for (int y = 0; y < maxS; y++) {
        if (!((mask >> y) & 1)) continue;
        fe cp[12]; cp[11] = c[12];
        for (int j = 11; j >= 1; j--) cp[j - 1] = fadd(c[j], fmul(Z[y], cp[j]));
        int all = 1;
        for (int e = 0; e < 4 && all; e++) {
            fe l = zero;
            for (int q = 0; q < 4; q++) { int j = M[q] - 1 - e; if (j >= 0 && j <= 11) l = fadd(l, fmul(cp[j], k[q])); }
            if (!fz(l)) all = 0;
        }
        if (all) return; // not canonical
    }
    (*canon)++;
}
static void rec(fe c[13][13], int d, int s, uint32_t mask, u64 *rk, u64 *canon) {
    if (d == 12) { leaf(c[12], mask, rk, canon); return; }
    for (int j = s; j <= 32 - (12 - d); j++) {
        fe *o = c[d], *nw = c[d + 1], zero = { 0, 0 };
        for (int t = 0; t <= d + 1; t++) { fe hi = t >= 1 ? o[t - 1] : zero, lo = t <= d ? fmul(o[t], Z[j]) : zero; nw[t] = fsub(hi, lo); }
        rec(c, d + 1, j + 1, mask | (1u << j), rk, canon);
    }
}
static void *worker(void *a) {
    (void)a; fe c[13][13]; u64 rk[4] = { 0 }, canon = 0;
    for (;;) {
        pthread_mutex_lock(&mu); int k = nexttask++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        c[0][0].r = 1; c[0][0].i = 0; memset(&c[0][1], 0, sizeof c[0][1]);
        rec(c, 0, 0, 0, rk, &canon); // placeholder, replaced below
        break;
    }
    (void)rk; (void)canon; return NULL;
}
// simpler driver: each task fixes the first two points
static void *worker2(void *a) {
    (void)a; fe c[13][13]; u64 rk[4] = { 0 }, canon = 0;
    for (;;) {
        pthread_mutex_lock(&mu); int k = nexttask++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int js[2] = { tj0[k], tj1[k] }; c[0][0].r = 1; c[0][0].i = 0; uint32_t mask = 0;
        for (int s = 0; s < 2; s++) {
            fe *o = c[s], *nw = c[s + 1], zero = { 0, 0 };
            for (int t = 0; t <= s + 1; t++) { fe hi = t >= 1 ? o[t - 1] : zero, lo = t <= s ? fmul(o[t], Z[js[s]]) : zero; nw[t] = fsub(hi, lo); }
            mask |= 1u << js[s];
        }
        rec(c, 2, js[1] + 1, mask, rk, &canon);
    }
    pthread_mutex_lock(&mu); for (int q = 0; q < 4; q++) g_rank[q] += rk[q]; g_canon += canon; pthread_mutex_unlock(&mu);
    return NULL;
}
int main(int argc, char **argv) {
    int nth = atoi(argv[1]), lim = atoi(argv[2]);
    for (int q = 0; q < 4; q++) M[q] = atoi(argv[3 + q]);
    fe w = { 906276279ULL, 1179735656ULL }; Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int a = 0; a < 32; a++) for (int b = a + 1; b < 32; b++) if (32 - b - 1 >= 10 && ntask < lim) { tj0[ntask] = a; tj1[ntask] = b; ntask++; }
    (void)worker;
    pthread_t th[64]; for (int k = 0; k < nth; k++) pthread_create(&th[k], 0, worker2, 0);
    for (int k = 0; k < nth; k++) pthread_join(th[k], 0);
    printf("support m=(%d,%d,%d,%d): rank3=%llu rank<=2=%llu canonical_genuine_rank3=%llu\n", M[0], M[1], M[2], M[3],
           (unsigned long long)g_rank[3], (unsigned long long)g_rank[2], (unsigned long long)g_canon);
    return 0;
}

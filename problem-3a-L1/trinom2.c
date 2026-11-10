// Distinct-member pass for every word with exactly three high monomials, RS[mu32, 17], agreement 20,
// over F_{p^2} = F_p[i], p = 2^31 - 1, mu32 = <w>, w = 906276279 + 1179735656 i.
//
// Same 3x3 test as trinom.c (see NOTES.md), plus a canonical-set filter so that each output record is
// one distinct list member:
//   a passing 20-set S (complement C) is canonical for its member iff no y in C with y < max(S)
//   makes the 21-set S + y pass. S + y passes iff l'_e = sum_j c'_j V_{e+1+j} = 0 for e = 0..3,
//   where prod_{t in C - y}(x - t) = sum c'_j x^j = Q_C(x) / (x - y).
// Every member has exactly one canonical 20-set (the 20 lowest-index points of its agreement set).
//
// Outputs:
//   argv[2] (binary): one record per canonical rank-2 hit at a genuine point (1 : b : g):
//                     u8 m1, m2, m3, pad; u32 b_re, b_im, g_re, g_im   (20 bytes)
//   argv[3] (text):   "L m1 m2 m3 n0r n0i n1r n1i n2r n2i mask" for rank-1 C that are canonical for the
//                     whole line (no y makes S + y pass on the whole line); "U m1 m2 m3 mask" for rank 0.
// For a point Pt of a triple, distinct(Pt) <= #records at Pt + #L lines through Pt + #U (exact on recount).
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
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }
static inline void cross(const fe a[3], const fe b[3], fe o[3]) {
    o[0] = fsub(fmul(a[1], b[2]), fmul(a[2], b[1]));
    o[1] = fsub(fmul(a[2], b[0]), fmul(a[0], b[2]));
    o[2] = fsub(fmul(a[0], b[1]), fmul(a[1], b[0]));
}

static fe Z[32];
static int task_j0[600], task_j1[600], ntask, next_task = 0;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static FILE *fbin, *ftxt;
static u64 leaves_total = 0, recs_total = 0, drop_total = 0;

typedef struct { unsigned char *b; size_t bl, bc; char *t; size_t tl, tc; u64 leaves, recs, drops; } tctx;
static void emitb(tctx *x, const void *s, size_t n) { if (x->bl + n > x->bc) { x->bc = (x->bc + n) * 2; x->b = realloc(x->b, x->bc); } memcpy(x->b + x->bl, s, n); x->bl += n; }
static void emitt(tctx *x, const char *s, size_t n) { if (x->tl + n > x->tc) { x->tc = (x->tc + n) * 2; x->t = realloc(x->t, x->tc); } memcpy(x->t + x->tl, s, n); x->tl += n; }

// rows of the 21-set test for complement C - y, triple (m1, m2, m3): r[e] = (c'_{m1-1-e}, c'_{m2-1-e}, c'_{m3-1-e})
static void rows21(const fe *c, int y, int m1, int m2, int m3, fe r[4][3]) {
    fe cp[12]; fe zero = { 0, 0 };
    cp[11] = c[12];
    for (int j = 11; j >= 1; j--) cp[j - 1] = fadd(c[j], fmul(Z[y], cp[j]));
    int ms[3] = { m1, m2, m3 };
    for (int e = 0; e < 4; e++) for (int q = 0; q < 3; q++) { int j = ms[q] - 1 - e; r[e][q] = (j < 0 || j > 11) ? zero : cp[j]; }
}

static void leaf(tctx *t, const fe *c, uint32_t mask) {
    t->leaves++;
    fe col[16][3]; fe zero = { 0, 0 };
    for (int m = 1; m <= 15; m++) for (int e = 0; e < 3; e++) { int j = m - 1 - e; col[m][e] = (j < 0 || j > 12) ? zero : c[j]; }
    int maxS = 31; while ((mask >> maxS) & 1) maxS--;
    char s[200];
    for (int m1 = 1; m1 <= 13; m1++) for (int m2 = m1 + 1; m2 <= 14; m2++) {
        fe w[3]; cross(col[m1], col[m2], w);
        for (int m3 = m2 + 1; m3 <= 15; m3++) {
            fe det = fadd(fadd(fmul(w[0], col[m3][0]), fmul(w[1], col[m3][1])), fmul(w[2], col[m3][2]));
            if (!fz(det)) continue;
            fe row[3][3];
            for (int e = 0; e < 3; e++) { row[e][0] = col[m1][e]; row[e][1] = col[m2][e]; row[e][2] = col[m3][e]; }
            fe k[3]; int found = 0; int pr[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
            for (int q = 0; q < 3 && !found; q++) { cross(row[pr[q][0]], row[pr[q][1]], k); if (!fz(k[0]) || !fz(k[1]) || !fz(k[2])) found = 1; }
            if (found) {
                if (fz(k[0]) || fz(k[1]) || fz(k[2])) continue;
                fe inv0 = finv(k[0]); fe pt[3] = { { 1, 0 }, fmul(k[1], inv0), fmul(k[2], inv0) };
                int canon = 1;
                for (int y = 0; y < maxS && canon; y++) {
                    if (!((mask >> y) & 1)) continue;
                    fe r[4][3]; rows21(c, y, m1, m2, m3, r);
                    int all = 1;
                    for (int e = 0; e < 4 && all; e++) {
                        fe l = fadd(fadd(r[e][0], fmul(r[e][1], pt[1])), fmul(r[e][2], pt[2]));
                        if (!fz(l)) all = 0;
                    }
                    if (all) canon = 0;
                }
                if (!canon) { t->drops++; continue; }
                unsigned char rec[20] = { (unsigned char)m1, (unsigned char)m2, (unsigned char)m3, 0 };
                uint32_t v[4] = { (uint32_t)pt[1].r, (uint32_t)pt[1].i, (uint32_t)pt[2].r, (uint32_t)pt[2].i };
                memcpy(rec + 4, v, 16); emitb(t, rec, 20); t->recs++;
                continue;
            }
            int r0 = -1;
            for (int e = 0; e < 3 && r0 < 0; e++) if (!fz(row[e][0]) || !fz(row[e][1]) || !fz(row[e][2])) r0 = e;
            int n;
            if (r0 < 0) { n = snprintf(s, sizeof s, "U %d %d %d %08x\n", m1, m2, m3, mask); emitt(t, s, n); continue; }
            int f = 0; while (fz(row[r0][f])) f++;
            fe inv = finv(row[r0][f]); fe nn[3]; for (int q = 0; q < 3; q++) nn[q] = fmul(row[r0][q], inv);
            int canon = 1;
            for (int y = 0; y < maxS && canon; y++) {
                if (!((mask >> y) & 1)) continue;
                fe r[4][3]; rows21(c, y, m1, m2, m3, r);
                int all = 1;
                for (int e = 0; e < 4 && all; e++) { fe x[3]; cross(r[e], nn, x); if (!fz(x[0]) || !fz(x[1]) || !fz(x[2])) all = 0; }
                if (all) canon = 0;
            }
            if (!canon) { t->drops++; continue; }
            n = snprintf(s, sizeof s, "L %d %d %d %llu %llu %llu %llu %llu %llu %08x\n", m1, m2, m3,
                         (unsigned long long)nn[0].r, (unsigned long long)nn[0].i, (unsigned long long)nn[1].r,
                         (unsigned long long)nn[1].i, (unsigned long long)nn[2].r, (unsigned long long)nn[2].i, mask);
            emitt(t, s, n);
        }
    }
}

static void rec(tctx *t, fe c[13][13], int depth, int start, uint32_t mask) {
    if (depth == 12) { leaf(t, c[12], mask); return; }
    for (int j = start; j <= 32 - (12 - depth); j++) {
        fe *o = c[depth], *nw = c[depth + 1]; fe zero = { 0, 0 };
        for (int d = 0; d <= depth + 1; d++) { fe hi = d >= 1 ? o[d - 1] : zero; fe lo = d <= depth ? fmul(o[d], Z[j]) : zero; nw[d] = fsub(hi, lo); }
        rec(t, c, depth + 1, j + 1, mask | (1u << j));
    }
}

static void flush(tctx *t) {
    pthread_mutex_lock(&mu);
    fwrite(t->b, 1, t->bl, fbin); fwrite(t->t, 1, t->tl, ftxt);
    pthread_mutex_unlock(&mu);
    t->bl = t->tl = 0;
}

static void *worker(void *arg) {
    (void)arg;
    tctx t; memset(&t, 0, sizeof t); fe c[13][13];
    for (;;) {
        pthread_mutex_lock(&mu); int k = next_task++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int js[2] = { task_j0[k], task_j1[k] };
        c[0][0].r = 1; c[0][0].i = 0; uint32_t mask = 0;
        for (int s = 0; s < 2; s++) {
            fe *o = c[s], *nw = c[s + 1]; fe zero = { 0, 0 };
            for (int d = 0; d <= s + 1; d++) { fe hi = d >= 1 ? o[d - 1] : zero; fe lo = d <= s ? fmul(o[d], Z[js[s]]) : zero; nw[d] = fsub(hi, lo); }
            mask |= 1u << js[s];
        }
        rec(&t, c, 2, js[1] + 1, mask);
        if (t.bl > (8u << 20) || t.tl > (1u << 20)) flush(&t);
    }
    flush(&t);
    pthread_mutex_lock(&mu); leaves_total += t.leaves; recs_total += t.recs; drop_total += t.drops; pthread_mutex_unlock(&mu);
    free(t.b); free(t.t);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: trinom2 THREADS points.bin lines.txt [task_limit]\n"); return 1; }
    int nth = atoi(argv[1]);
    int limit = argc > 4 ? atoi(argv[4]) : 1 << 30;
    fbin = fopen(argv[2], "wb"); ftxt = fopen(argv[3], "w");
    fe w = { 906276279ULL, 1179735656ULL };
    Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int j0 = 0; j0 < 32; j0++) for (int j1 = j0 + 1; j1 < 32; j1++)
        if (32 - j1 - 1 >= 10 && ntask < limit) { task_j0[ntask] = j0; task_j1[ntask] = j1; ntask++; }
    pthread_t th[64];
    for (int k = 0; k < nth; k++) pthread_create(&th[k], NULL, worker, NULL);
    for (int k = 0; k < nth; k++) pthread_join(th[k], NULL);
    fclose(fbin); fclose(ftxt);
    fprintf(stderr, "leaves=%llu tasks=%d canonical_point_records=%llu non_canonical_dropped=%llu\n",
            (unsigned long long)leaves_total, ntask, (unsigned long long)recs_total, (unsigned long long)drop_total);
    return 0;
}

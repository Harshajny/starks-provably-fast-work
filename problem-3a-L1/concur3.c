// Concurrency in P^3 over F_{p^2} for the low-rank sets dumped by quad.c (low.txt).
// Planes ("1 n"): weight of a point = number of plane sets containing it. Each point lying on >= 2
// distinct planes is found exactly once with its full weight: on the line of its two smallest-index
// distinct planes i < j, tallying only planes k > j (plus planes containing that whole line).
// Lines ("2 a b"): all pairs of distinct lines are intersected; weight = number of line sets through it.
// Prints max plane weight at any genuine point (all four coordinates nonzero), max line weight, and the
// genuine points with plane weight >= H (for an exact recount).
// Usage: ./concur3 low.txt THREADS H > out
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
static inline fe fneg(fe a) { fe z = { (P - a.r) % P, (P - a.i) % P }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static inline int feq(fe a, fe b) { return a.r == b.r && a.i == b.i; }
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }
static inline fe dot(const fe *a, const fe *b) { return fadd(fadd(fmul(a[0], b[0]), fmul(a[1], b[1])), fadd(fmul(a[2], b[2]), fmul(a[3], b[3]))); }
static void normalize(fe *v) { int f = 0; while (f < 4 && fz(v[f])) f++; if (f == 4) return; fe inv = finv(v[f]); for (int q = 0; q < 4; q++) v[q] = fmul(v[q], inv); }

typedef struct { fe n[4]; int w; } plane_t;
static plane_t *PL; static int NP;
static int H, next_i; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static long g_best = 0, g_bestbase = 0; static fe g_bestpt[4];
static int cmpn(const void *a, const void *b) { return memcmp(((const plane_t *)a)->n, ((const plane_t *)b)->n, sizeof(fe) * 4); }

// kernel basis of the 2x4 matrix [a; b] (rank 2): two vectors p0, p1
static int kernel2(const fe *a, const fe *b, fe p0[4], fe p1[4]) {
    fe A[2][4]; memcpy(A[0], a, sizeof(fe) * 4); memcpy(A[1], b, sizeof(fe) * 4);
    int piv[2], r = 0;
    for (int c = 0; c < 4 && r < 2; c++) {
        int s = -1; for (int i = r; i < 2; i++) if (!fz(A[i][c])) { s = i; break; }
        if (s < 0) continue;
        if (s != r) for (int k = 0; k < 4; k++) { fe t = A[r][k]; A[r][k] = A[s][k]; A[s][k] = t; }
        fe inv = finv(A[r][c]); for (int k = 0; k < 4; k++) A[r][k] = fmul(A[r][k], inv);
        for (int i = 0; i < 2; i++) if (i != r && !fz(A[i][c])) { fe f = A[i][c]; for (int k = 0; k < 4; k++) A[i][k] = fsub(A[i][k], fmul(f, A[r][k])); }
        piv[r++] = c;
    }
    if (r < 2) return 0;
    int fr[2], nf = 0; for (int c = 0; c < 4; c++) if (c != piv[0] && c != piv[1]) fr[nf++] = c;
    fe *out[2] = { p0, p1 };
    for (int t = 0; t < 2; t++) {
        for (int q = 0; q < 4; q++) out[t][q] = (fe){ 0, 0 };
        out[t][fr[t]] = (fe){ 1, 0 };
        for (int i = 0; i < 2; i++) out[t][piv[i]] = fneg(A[i][fr[t]]);
    }
    return 1;
}

typedef struct { fe t; long w; int used; } slot_t;
static void *worker(void *arg) {
    (void)arg;
    int cap = 1; while (cap < 4 * NP) cap <<= 1;
    slot_t *tab = malloc(sizeof(slot_t) * cap);
    fe *num = malloc(sizeof(fe) * NP), *den = malloc(sizeof(fe) * NP), *pre = malloc(sizeof(fe) * NP);
    int *wk = malloc(sizeof(int) * NP);
    long best = 0, bestbase = 0; fe bestpt[4] = { { 0, 0 } };
    for (;;) {
        pthread_mutex_lock(&mu); int i = next_i++; pthread_mutex_unlock(&mu);
        if (i >= NP) break;
        for (int j = i + 1; j < NP; j++) {
            fe p0[4], p1[4];
            if (!kernel2(PL[i].n, PL[j].n, p0, p1)) continue;
            long base = PL[i].w + PL[j].w; int m = 0; long inf = 0;
            for (int k = j + 1; k < NP; k++) {
                fe d0 = dot(PL[k].n, p0), d1 = dot(PL[k].n, p1);
                if (fz(d1)) { if (fz(d0)) base += PL[k].w; else inf += PL[k].w; continue; }
                num[m] = fneg(d0); den[m] = d1; wk[m] = PL[k].w; m++;
            }
            if (base > bestbase) bestbase = base;
            if (m == 0 && inf == 0) continue;
            // batch inversion of den
            fe acc = { 1, 0 };
            for (int q = 0; q < m; q++) { pre[q] = acc; acc = fmul(acc, den[q]); }
            fe ia = finv(acc);
            for (int q = m - 1; q >= 0; q--) { fe inv = fmul(ia, pre[q]); ia = fmul(ia, den[q]); num[q] = fmul(num[q], inv); }
            for (int q = 0; q < cap; q++) tab[q].used = 0;
            for (int q = 0; q < m; q++) {
                u64 h = (num[q].r * 0x9E3779B97F4A7C15ULL ^ num[q].i * 0xC2B2AE3D27D4EB4FULL) & (cap - 1);
                while (tab[h].used && !feq(tab[h].t, num[q])) h = (h + 1) & (cap - 1);
                if (!tab[h].used) { tab[h].used = 1; tab[h].t = num[q]; tab[h].w = 0; }
                tab[h].w += wk[q];
            }
            // candidates: the t = infinity point p1, and every tallied t
            for (int q = -1; q < cap; q++) {
                long w; fe pt[4];
                if (q < 0) { if (!inf) continue; w = base + inf; memcpy(pt, p1, sizeof pt); }
                else { if (!tab[q].used) continue; w = base + tab[q].w; for (int c = 0; c < 4; c++) pt[c] = fadd(p0[c], fmul(tab[q].t, p1[c])); }
                if (w < 3 || (w <= best && w < H)) continue;
                if (fz(pt[0]) || fz(pt[1]) || fz(pt[2]) || fz(pt[3])) continue; // not a four-monomial word
                normalize(pt);
                if (w > best) { best = w; memcpy(bestpt, pt, sizeof pt); }
                if (w >= H) {
                    pthread_mutex_lock(&mu);
                    printf("C %ld", w); for (int c = 1; c < 4; c++) printf(" %llu %llu", (unsigned long long)pt[c].r, (unsigned long long)pt[c].i); printf("\n");
                    pthread_mutex_unlock(&mu);
                }
            }
        }
    }
    pthread_mutex_lock(&mu);
    if (best > g_best) { g_best = best; memcpy(g_bestpt, bestpt, sizeof bestpt); }
    if (bestbase > g_bestbase) g_bestbase = bestbase;
    pthread_mutex_unlock(&mu);
    free(tab); free(num); free(den); free(pre); free(wk);
    return NULL;
}

int main(int argc, char **argv) {
    int nth = atoi(argv[2]); H = atoi(argv[3]);
    FILE *f = fopen(argv[1], "r"); char line[2048];
    int cap = 1 << 14; PL = malloc(sizeof(plane_t) * cap);
    fe (*LN)[2][4] = malloc(sizeof(fe[2][4]) * 65536); int NL = 0;
    while (fgets(line, sizeof line, f)) {
        int rk = line[0] - '0'; char *p = line + 1; fe v[8];
        for (int q = 0; q < 4 * rk; q++) { v[q].r = strtoull(p, &p, 10); v[q].i = strtoull(p, &p, 10); }
        if (rk == 1) { normalize(v); memcpy(PL[NP].n, v, sizeof(fe) * 4); PL[NP].w = 1; NP++; }
        else if (rk == 2) { memcpy(LN[NL][0], v, sizeof(fe) * 4); memcpy(LN[NL][1], v + 4, sizeof(fe) * 4); NL++; }
    }
    qsort(PL, NP, sizeof *PL, cmpn);
    int k = 0; for (int i = 0; i < NP; i++) { if (k && !cmpn(&PL[k - 1], &PL[i])) PL[k - 1].w++; else PL[k++] = PL[i]; }
    fprintf(stderr, "plane sets %d, distinct planes %d; line sets %d\n", NP, k, NL); NP = k;
    pthread_t th[64]; for (int t = 0; t < nth; t++) pthread_create(&th[t], 0, worker, 0);
    for (int t = 0; t < nth; t++) pthread_join(th[t], 0);
    // line sets: pairwise intersections (lines as kernels of two rows)
    long lbest = 1;
    typedef struct { fe pt[4]; } lp_t; lp_t *pts = malloc(sizeof(lp_t) * 1); size_t npt = 0, pcap = 1;
    for (int a = 0; a < NL; a++) for (int b = a + 1; b < NL; b++) {
        fe p0[4], p1[4]; if (!kernel2(LN[a][0], LN[a][1], p0, p1)) continue;
        fe d0 = dot(LN[b][0], p0), d1 = dot(LN[b][0], p1), e0 = dot(LN[b][1], p0), e1 = dot(LN[b][1], p1);
        // point p0 + t p1 on line b: both equations; nontrivial solution iff 2x2 det = 0
        if (!fz(fsub(fmul(d0, e1), fmul(d1, e0)))) continue;
        fe pt[4];
        if (!fz(d1)) { fe t = fmul(fneg(d0), finv(d1)); for (int c = 0; c < 4; c++) pt[c] = fadd(p0[c], fmul(t, p1[c])); }
        else if (!fz(e1)) { fe t = fmul(fneg(e0), finv(e1)); for (int c = 0; c < 4; c++) pt[c] = fadd(p0[c], fmul(t, p1[c])); }
        else if (fz(d0) && fz(e0)) { lbest = lbest < 2 ? 2 : lbest; continue; } // same line (rare): counted as 2
        else memcpy(pt, p1, sizeof pt);
        normalize(pt);
        if (npt == pcap) { pcap *= 2; pts = realloc(pts, sizeof(lp_t) * pcap); }
        pts[npt++] = (lp_t){ { pt[0], pt[1], pt[2], pt[3] } };
    }
    qsort(pts, npt, sizeof(lp_t), (int (*)(const void *, const void *))cmpn);
    for (size_t a = 0; a < npt;) { size_t b = a; while (b < npt && !memcmp(&pts[b], &pts[a], sizeof(lp_t))) b++;
        long pairs = (long)(b - a); long L = 2; while (L * (L - 1) / 2 < pairs) L++; if (L > lbest) lbest = L;
        if (L >= H) { fe *q = pts[a].pt; if (!fz(q[0]) && !fz(q[1]) && !fz(q[2]) && !fz(q[3])) { printf("LC %ld", L); for (int c = 1; c < 4; c++) printf(" %llu %llu", (unsigned long long)q[c].r, (unsigned long long)q[c].i); printf("\n"); } }
        a = b; }
    printf("SUMMARY distinct_planes=%d max_plane_weight_genuine=%ld max_line_base=%ld max_line_set_concurrency=%ld line_sets=%d\n",
           NP, g_best, g_bestbase, lbest, NL);
    printf("BESTPT"); for (int c = 1; c < 4; c++) printf(" %llu %llu", (unsigned long long)g_bestpt[c].r, (unsigned long long)g_bestpt[c].i); printf("\n");
    return 0;
}

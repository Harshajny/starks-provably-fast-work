// For each exponent triple, find the genuine points (1 : b : g) where the weighted line records of
// trinom2.c concur. Two distinct lines meet in exactly one projective point (cross product of normals).
// For a point on L distinct lines, the pairs meeting there number L(L-1)/2 and their weight sums total
// (L-1) * (sum of the L weights), so both are recovered exactly from the pair records.
// Input: tri_lines.txt. Output (stdout): per triple the max line weight through one point, the max
// single-line weight, and every point whose line weight is >= THRESH as "C m1 m2 m3 br bi gr gi weight nlines".
// Usage: ./linesmax tri_lines.txt THRESH
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
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }

typedef struct { int tr; fe n[3]; int w; } line_t;
typedef struct { uint32_t b[4]; int ws; } pr_t;

static int cmp_line(const void *a, const void *b) {
    const line_t *x = a, *y = b;
    if (x->tr != y->tr) return x->tr < y->tr ? -1 : 1;
    return memcmp(x->n, y->n, sizeof x->n);
}
static int cmp_pr(const void *a, const void *b) { return memcmp(((const pr_t *)a)->b, ((const pr_t *)b)->b, 16); }

int main(int argc, char **argv) {
    int thresh = argc > 2 ? atoi(argv[2]) : 100;
    FILE *f = fopen(argv[1], "r");
    size_t cap = 1 << 20, n = 0; line_t *L = malloc(cap * sizeof *L);
    char kind[4]; int m1, m2, m3; unsigned long long v[6]; unsigned mask;
    char buf[256];
    while (fgets(buf, sizeof buf, f)) {
        if (buf[0] != 'L') continue;
        sscanf(buf, "%s %d %d %d %llu %llu %llu %llu %llu %llu %x", kind, &m1, &m2, &m3, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &mask);
        if (n == cap) { cap *= 2; L = realloc(L, cap * sizeof *L); }
        line_t *l = &L[n++]; memset(l, 0, sizeof *l);
        l->tr = m1 * 256 + m2 * 16 + m3;
        for (int q = 0; q < 3; q++) { l->n[q].r = v[2 * q]; l->n[q].i = v[2 * q + 1]; }
        l->w = 1;
    }
    qsort(L, n, sizeof *L, cmp_line);
    size_t k = 0; // merge duplicates into weights
    for (size_t i = 0; i < n; i++) {
        if (k && L[k - 1].tr == L[i].tr && !memcmp(L[k - 1].n, L[i].n, sizeof L[i].n)) L[k - 1].w++;
        else L[k++] = L[i];
    }
    n = k;
    fprintf(stderr, "distinct lines: %zu\n", n);
    size_t s = 0;
    while (s < n) {
        size_t e = s; while (e < n && L[e].tr == L[s].tr) e++;
        int tr = L[s].tr, maxw1 = 0;
        for (size_t i = s; i < e; i++) if (L[i].w > maxw1) maxw1 = L[i].w;
        size_t m = e - s, np = 0, pcap = m * (m - 1) / 2 + 1;
        pr_t *PR = malloc(pcap * sizeof *PR);
        for (size_t i = s; i < e; i++) for (size_t j = i + 1; j < e; j++) {
            const fe *a = L[i].n, *b = L[j].n; fe c[3];
            c[0] = fsub(fmul(a[1], b[2]), fmul(a[2], b[1]));
            c[1] = fsub(fmul(a[2], b[0]), fmul(a[0], b[2]));
            c[2] = fsub(fmul(a[0], b[1]), fmul(a[1], b[0]));
            if (fz(c[0]) || fz(c[1]) || fz(c[2])) continue; // not a genuine trinomial point
            fe inv = finv(c[0]); fe bb = fmul(c[1], inv), gg = fmul(c[2], inv);
            pr_t *p = &PR[np++];
            p->b[0] = (uint32_t)bb.r; p->b[1] = (uint32_t)bb.i; p->b[2] = (uint32_t)gg.r; p->b[3] = (uint32_t)gg.i;
            p->ws = L[i].w + L[j].w;
        }
        qsort(PR, np, sizeof *PR, cmp_pr);
        long best = 0; int bestL = 0;
        size_t a = 0;
        while (a < np) {
            size_t b = a; long sum = 0;
            while (b < np && !memcmp(PR[b].b, PR[a].b, 16)) { sum += PR[b].ws; b++; }
            long pairs = (long)(b - a);
            int Ln = 2; while ((long)Ln * (Ln - 1) / 2 < pairs) Ln++;
            if ((long)Ln * (Ln - 1) / 2 != pairs) fprintf(stderr, "pair count %ld is not triangular\n", pairs);
            long wsum = sum / (Ln - 1);
            if (wsum > best) { best = wsum; bestL = Ln; }
            if (wsum >= thresh)
                printf("C %d %d %d %u %u %u %u %ld %d\n", tr / 256, (tr / 16) % 16, tr % 16, PR[a].b[0], PR[a].b[1], PR[a].b[2], PR[a].b[3], wsum, Ln);
            a = b;
        }
        printf("T %d %d %d distinct_lines=%zu max_single_line_weight=%d max_concurrent_weight=%ld (%d lines)\n",
               tr / 256, (tr / 16) % 16, tr % 16, m, maxw1, best, bestL);
        free(PR);
        s = e;
    }
    return 0;
}

// Exact plane-set and line-set weight at given words (1 : b : c : d), from quad.c's low.txt.
// Usage: ./weigh low.txt points.txt   (points: "br bi cr ci dr di" per line) -> "planes lines coords"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
typedef uint64_t u64;
#define P 2147483647ULL
typedef struct { u64 r, i; } fe;
static inline u64 red(u64 x) { x = (x & P) + (x >> 31); x = (x & P) + (x >> 31); return x >= P ? x - P : x; }
static inline fe fmul(fe a, fe b) { u64 ac = red(a.r * b.r), bd = red(a.i * b.i); fe z = { red(ac + P - bd), red(a.r * b.i + a.i * b.r) }; return z; }
static inline fe fadd(fe a, fe b) { fe z = { red(a.r + b.r), red(a.i + b.i) }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static inline fe dot(const fe *a, const fe *b) { return fadd(fadd(fmul(a[0], b[0]), fmul(a[1], b[1])), fadd(fmul(a[2], b[2]), fmul(a[3], b[3]))); }
int main(int argc, char **argv) {
    static fe R[30000][2][4]; static int rk[30000]; int n = 0; char line[2048];
    FILE *f = fopen(argv[1], "r");
    while (fgets(line, sizeof line, f)) { rk[n] = line[0] - '0'; char *p = line + 1;
        for (int e = 0; e < rk[n]; e++) for (int q = 0; q < 4; q++) { R[n][e][q].r = strtoull(p, &p, 10); R[n][e][q].i = strtoull(p, &p, 10); } n++; }
    FILE *g = fopen(argv[2], "r"); unsigned long long v[6];
    while (fscanf(g, "%llu %llu %llu %llu %llu %llu", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]) == 6) {
        fe pt[4] = { { 1, 0 }, { v[0], v[1] }, { v[2], v[3] }, { v[4], v[5] } }; int pl = 0, ln = 0;
        for (int s = 0; s < n; s++) { int in = 1; for (int e = 0; e < rk[s]; e++) if (!fz(dot(R[s][e], pt))) in = 0; if (in) { if (rk[s] == 1) pl++; else ln++; } }
        printf("%d %d %llu %llu %llu %llu %llu %llu\n", pl, ln, v[0], v[1], v[2], v[3], v[4], v[5]);
    }
    return 0;
}

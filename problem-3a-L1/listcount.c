// Exact distinct list size, RS[mu32, 17], agreement 20, over F_{p^2}, for several words at once.
// A word is given by its syndrome V_1..V_15 (only nonzero entries listed). For every 12-set C:
//   l_e = sum_j c_j(C) V_{e+1+j}, e = 0..2, computed incrementally along the enumeration:
//   start l_e = V_{e+1} (e = 0..14); adding a point t maps l_e -> l_{e+1} - t l_e.
// C passes iff l_0 = l_1 = l_2 = 0. A passing 20-set S is counted as a distinct member iff it is
// canonical: no y in C with y < max(S) makes S + y pass (4 equations on Q_C / (x - y)).
// Input file: one word per line: "name m re im [m re im ...]" (V_m = re + i*im).
// Output: per word, passing 20-sets and distinct members, plus the agreement-size histogram.
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

#define MAXW 128
static int K; static fe V[MAXW][16]; static char names[MAXW][64];
static fe Z[32];
static int task_j0[600], task_j1[600], ntask, next_task = 0;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static u64 passing[MAXW], distinct[MAXW], agree_hist[MAXW][33];

// agreement size of the member of a passing 20-set: 20 + #{y in C : S + y passes}
static int canon_and_agree(int w, uint32_t mask, int *is_canon) {
    fe c[13]; memset(c, 0, sizeof c); c[0].r = 1; int deg = 0;
    for (int j = 0; j < 32; j++) if ((mask >> j) & 1) {
        for (int d = deg + 1; d >= 0; d--) { fe hi = d >= 1 ? c[d - 1] : (fe){ 0, 0 }; fe lo = d <= deg ? fmul(c[d], Z[j]) : (fe){ 0, 0 }; c[d] = fsub(hi, lo); }
        deg++;
    }
    int maxS = 31; while ((mask >> maxS) & 1) maxS--;
    int a = 20; *is_canon = 1;
    for (int y = 0; y < 32; y++) {
        if (!((mask >> y) & 1)) continue;
        fe cp[12]; cp[11] = c[12];
        for (int j = 11; j >= 1; j--) cp[j - 1] = fadd(c[j], fmul(Z[y], cp[j]));
        int ok = 1;
        for (int e = 0; e < 4 && ok; e++) {
            fe l = { 0, 0 };
            for (int j = 0; j <= 11; j++) { int m = e + 1 + j; if (m <= 15 && !fz(V[w][m])) l = fadd(l, fmul(cp[j], V[w][m])); }
            if (!fz(l)) ok = 0;
        }
        if (ok) { a++; if (y < maxS) *is_canon = 0; }
    }
    return a;
}

typedef struct { u64 pass[MAXW], dist[MAXW], hist[MAXW][33]; } acc_t;

static void rec(acc_t *A, fe (*l)[MAXW][15], int depth, int start, uint32_t mask) {
    if (depth == 12) {
        for (int w = 0; w < K; w++) {
            if (!fz(l[12][w][0]) || !fz(l[12][w][1]) || !fz(l[12][w][2])) continue;
            A->pass[w]++;
            int canon; int a = canon_and_agree(w, mask, &canon);
            if (canon) { A->dist[w]++; A->hist[w][a]++; }
        }
        return;
    }
    int ne = 15 - depth - 1; // number of l's needed after this step: 3 + (12 - depth - 1)
    for (int j = start; j <= 32 - (12 - depth); j++) {
        for (int w = 0; w < K; w++)
            for (int e = 0; e < ne; e++) l[depth + 1][w][e] = fsub(l[depth][w][e + 1], fmul(Z[j], l[depth][w][e]));
        rec(A, l, depth + 1, j + 1, mask | (1u << j));
    }
}

static void *worker(void *arg) {
    (void)arg;
    acc_t *A = calloc(1, sizeof *A);
    fe (*l)[MAXW][15] = malloc(13 * sizeof *l);
    for (;;) {
        pthread_mutex_lock(&mu); int k = next_task++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int js[2] = { task_j0[k], task_j1[k] };
        for (int w = 0; w < K; w++) for (int e = 0; e < 15; e++) l[0][w][e] = V[w][e + 1];
        for (int s = 0; s < 2; s++) for (int w = 0; w < K; w++) for (int e = 0; e < 14 - s; e++)
            l[s + 1][w][e] = fsub(l[s][w][e + 1], fmul(Z[js[s]], l[s][w][e]));
        rec(A, l, 2, js[1] + 1, (1u << js[0]) | (1u << js[1]));
    }
    pthread_mutex_lock(&mu);
    for (int w = 0; w < K; w++) { passing[w] += A->pass[w]; distinct[w] += A->dist[w]; for (int a = 0; a < 33; a++) agree_hist[w][a] += A->hist[w][a]; }
    pthread_mutex_unlock(&mu);
    free(A); free(l);
    return NULL;
}

int main(int argc, char **argv) {
    int nth = argc > 2 ? atoi(argv[2]) : 10;
    FILE *f = fopen(argv[1], "r"); char line[1024];
    while (fgets(line, sizeof line, f) && K < MAXW) {
        char *tok = strtok(line, " \t\n"); if (!tok || tok[0] == '#') continue;
        snprintf(names[K], 64, "%s", tok);
        memset(V[K], 0, sizeof V[K]);
        char *a, *b, *c;
        while ((a = strtok(NULL, " \t\n")) && (b = strtok(NULL, " \t\n")) && (c = strtok(NULL, " \t\n"))) {
            int m = atoi(a); V[K][m].r = strtoull(b, 0, 10) % P; V[K][m].i = strtoull(c, 0, 10) % P;
        }
        K++;
    }
    fe w = { 906276279ULL, 1179735656ULL };
    Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int j0 = 0; j0 < 32; j0++) for (int j1 = j0 + 1; j1 < 32; j1++)
        if (32 - j1 - 1 >= 10) { task_j0[ntask] = j0; task_j1[ntask] = j1; ntask++; }
    pthread_t th[64];
    for (int k = 0; k < nth; k++) pthread_create(&th[k], NULL, worker, NULL);
    for (int k = 0; k < nth; k++) pthread_join(th[k], NULL);
    for (int k = 0; k < K; k++) {
        printf("%s passing=%llu distinct=%llu agreement:", names[k], (unsigned long long)passing[k], (unsigned long long)distinct[k]);
        for (int a = 20; a <= 32; a++) if (agree_hist[k][a]) printf(" %d:%llu", a, (unsigned long long)agree_hist[k][a]);
        printf("\n");
    }
    return 0;
}

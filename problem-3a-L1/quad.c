// Four-monomial support (m1..m4), RS[mu32, 17], agreement 20, over F_{p^2}, p = 2^31 - 1.
// For each 12-set C the 3x4 system M_C a = 0 (rows e = 0,1,2: c_{mq-1-e}) decides which words
// a = (a1 : a2 : a3 : a4) (syndrome V_mq = aq) the 20-set S = Z - C passes for.
//   rank 3: exactly one word (kernel = signed 3x3 minors). Kept if all four coordinates are nonzero
//           (otherwise it is a <= 3-monomial word, covered by bp1_c90aa8fa / bp1_cc3e6eb2) and
//           canonical (no y in C, y < max S, makes S + y pass at that word).
//   rank <= 2: a line or plane of words. Kept if canonical for the whole kernel K: no y with
//           rank([M_C; M'_{C-y}]) == rank(M_C) (then S + y passes on all of K). Counted, by rank.
// distinct(word) <= R(word) + (rank <= 2 sets whose kernel contains the word) <= R(word) + W_low.
// Pass 1 tallies R in a hashed sketch (2^29 u16); pass 2 writes every rank-3 record whose bucket
// is >= T (a bucket never undercounts), so every word with R >= T is written with all its records.
// Usage: ./quad THREADS T m1 m2 m3 m4 heavy.bin [low.txt]   (low.txt: kept rank<=2 sets as rows)
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
static u64 powp(u64 b, u64 e) { u64 r = 1; b %= P; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = powp(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }

static fe Z[32]; static int M[4], T, PASS;
static int tj0[600], tj1[600], ntask, nexttask; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
#define SK_BITS 29
static uint16_t *sk; static FILE *fout, *flow;
static u64 g_r3, g_keep3, g_low[3], g_lowkeep[3], g_written;

static inline fe det3(const fe *a, const fe *b, const fe *c) {
    fe x = fsub(fmul(b[1], c[2]), fmul(b[2], c[1])), y = fsub(fmul(b[2], c[0]), fmul(b[0], c[2])), z = fsub(fmul(b[0], c[1]), fmul(b[1], c[0]));
    return fadd(fadd(fmul(a[0], x), fmul(a[1], y)), fmul(a[2], z));
}
static int rank_rows(fe A[][4], int n) { // Gaussian elimination, destroys A
    int r = 0;
    for (int col = 0; col < 4 && r < n; col++) {
        int piv = -1; for (int i = r; i < n; i++) if (!fz(A[i][col])) { piv = i; break; }
        if (piv < 0) continue;
        for (int k = 0; k < 4; k++) { fe t = A[r][k]; A[r][k] = A[piv][k]; A[piv][k] = t; }
        fe inv = finv(A[r][col]);
        for (int i = 0; i < n; i++) if (i != r && !fz(A[i][col])) {
            fe f = fmul(A[i][col], inv);
            for (int k = 0; k < 4; k++) A[i][k] = fsub(A[i][k], fmul(f, A[r][k]));
        }
        r++;
    }
    return r;
}
static inline u64 hpt(const fe *pt) {
    u64 h = 0x9E3779B97F4A7C15ULL;
    for (int q = 1; q < 4; q++) { h ^= pt[q].r + 0x632BE59BD9B4E019ULL + (h << 6) + (h >> 2); h *= 0xBF58476D1CE4E5B9ULL; h ^= pt[q].i + (h << 6) + (h >> 2); h *= 0x94D049BB133111EBULL; }
    return h >> (64 - SK_BITS);
}

typedef struct { u64 r3, keep3, low[3], lowkeep[3], written; unsigned char *buf; size_t bl, bc; } acc_t;

static void rows21(const fe *c, int y, fe R[4][4]) {
    fe cp[12], zero = { 0, 0 }; cp[11] = c[12];
    for (int j = 11; j >= 1; j--) cp[j - 1] = fadd(c[j], fmul(Z[y], cp[j]));
    for (int e = 0; e < 4; e++) for (int q = 0; q < 4; q++) { int j = M[q] - 1 - e; R[e][q] = (j < 0 || j > 11) ? zero : cp[j]; }
}

static void leaf(const fe *c, uint32_t mask, acc_t *A) {
    fe col[4][3], zero = { 0, 0 };
    for (int q = 0; q < 4; q++) for (int e = 0; e < 3; e++) { int j = M[q] - 1 - e; col[q][e] = (j < 0 || j > 12) ? zero : c[j]; }
    fe k[4] = { det3(col[1], col[2], col[3]), fneg(det3(col[0], col[2], col[3])), det3(col[0], col[1], col[3]), fneg(det3(col[0], col[1], col[2])) };
    int maxS = 31; while ((mask >> maxS) & 1) maxS--;
    if (!(fz(k[0]) && fz(k[1]) && fz(k[2]) && fz(k[3]))) {
        A->r3++;
        if (fz(k[0]) || fz(k[1]) || fz(k[2]) || fz(k[3])) return;
        for (int y = 0; y < maxS; y++) {
            if (!((mask >> y) & 1)) continue;
            fe R[4][4]; rows21(c, y, R);
            int all = 1;
            for (int e = 0; e < 4 && all; e++) {
                fe l = fadd(fadd(fmul(R[e][0], k[0]), fmul(R[e][1], k[1])), fadd(fmul(R[e][2], k[2]), fmul(R[e][3], k[3])));
                if (!fz(l)) all = 0;
            }
            if (all) return;
        }
        A->keep3++;
        fe inv = finv(k[0]); fe pt[4] = { { 1, 0 }, fmul(k[1], inv), fmul(k[2], inv), fmul(k[3], inv) };
        u64 h = hpt(pt);
        if (PASS == 1) { __atomic_fetch_add(&sk[h], (uint16_t)(sk[h] < 65000), __ATOMIC_RELAXED); return; }
        if (sk[h] >= T) {
            uint32_t rec[6] = { (uint32_t)pt[1].r, (uint32_t)pt[1].i, (uint32_t)pt[2].r, (uint32_t)pt[2].i, (uint32_t)pt[3].r, (uint32_t)pt[3].i };
            if (A->bl + 24 > A->bc) { A->bc = (A->bc + 24) * 2; A->buf = realloc(A->buf, A->bc); }
            memcpy(A->buf + A->bl, rec, 24); A->bl += 24; A->written++;
        }
        return;
    }
    if (PASS == 1) return;
    fe B[3][4]; for (int e = 0; e < 3; e++) for (int q = 0; q < 4; q++) B[e][q] = col[q][e];
    int rk = rank_rows(B, 3);
    A->low[rk]++;
    for (int y = 0; y < maxS; y++) {
        if (!((mask >> y) & 1)) continue;
        fe R[4][4]; rows21(c, y, R);
        fe X[7][4]; for (int e = 0; e < 3; e++) for (int q = 0; q < 4; q++) X[e][q] = col[q][e];
        for (int e = 0; e < 4; e++) for (int q = 0; q < 4; q++) X[3 + e][q] = R[e][q];
        if (rank_rows(X, 7) == rk) return; // S + y passes on the whole kernel: not canonical
    }
    A->lowkeep[rk]++;
    if (flow) { // dump the kernel's defining rows (row-echelon basis of M_C): 1 row = plane, 2 rows = line
        pthread_mutex_lock(&mu);
        fprintf(flow, "%d", rk);
        for (int e = 0; e < rk; e++) for (int q = 0; q < 4; q++) fprintf(flow, " %llu %llu", (unsigned long long)B[e][q].r, (unsigned long long)B[e][q].i);
        fprintf(flow, "\n");
        pthread_mutex_unlock(&mu);
    }
}
static void rec(fe c[13][13], int d, int s, uint32_t mask, acc_t *A) {
    if (d == 12) { leaf(c[12], mask, A); return; }
    for (int j = s; j <= 32 - (12 - d); j++) {
        fe *o = c[d], *nw = c[d + 1], zero = { 0, 0 };
        for (int t = 0; t <= d + 1; t++) { fe hi = t >= 1 ? o[t - 1] : zero, lo = t <= d ? fmul(o[t], Z[j]) : zero; nw[t] = fsub(hi, lo); }
        rec(c, d + 1, j + 1, mask | (1u << j), A);
    }
}
static void *worker(void *a) {
    (void)a; fe c[13][13]; acc_t A; memset(&A, 0, sizeof A);
    for (;;) {
        pthread_mutex_lock(&mu); int k = nexttask++; pthread_mutex_unlock(&mu);
        if (k >= ntask) break;
        int js[2] = { tj0[k], tj1[k] }; c[0][0].r = 1; c[0][0].i = 0; uint32_t mask = 0;
        for (int s = 0; s < 2; s++) {
            fe *o = c[s], *nw = c[s + 1], zero = { 0, 0 };
            for (int t = 0; t <= s + 1; t++) { fe hi = t >= 1 ? o[t - 1] : zero, lo = t <= s ? fmul(o[t], Z[js[s]]) : zero; nw[t] = fsub(hi, lo); }
            mask |= 1u << js[s];
        }
        rec(c, 2, js[1] + 1, mask, &A);
        if (A.bl > (8u << 20)) { pthread_mutex_lock(&mu); fwrite(A.buf, 1, A.bl, fout); pthread_mutex_unlock(&mu); A.bl = 0; }
    }
    pthread_mutex_lock(&mu);
    if (PASS == 2) fwrite(A.buf, 1, A.bl, fout);
    g_r3 += A.r3; g_keep3 += A.keep3; g_written += A.written;
    for (int r = 0; r < 3; r++) { g_low[r] += A.low[r]; g_lowkeep[r] += A.lowkeep[r]; }
    pthread_mutex_unlock(&mu); free(A.buf);
    return NULL;
}
static void run(int nth) {
    nexttask = 0; g_r3 = g_keep3 = g_written = 0; memset(g_low, 0, sizeof g_low); memset(g_lowkeep, 0, sizeof g_lowkeep);
    pthread_t th[64]; for (int k = 0; k < nth; k++) pthread_create(&th[k], 0, worker, 0);
    for (int k = 0; k < nth; k++) pthread_join(th[k], 0);
}
int main(int argc, char **argv) {
    int nth = atoi(argv[1]); T = atoi(argv[2]);
    for (int q = 0; q < 4; q++) M[q] = atoi(argv[3 + q]);
    fout = fopen(argv[7], "wb"); flow = argc > 8 ? fopen(argv[8], "w") : NULL;
    fe w = { 906276279ULL, 1179735656ULL }; Z[0].r = 1; Z[0].i = 0;
    for (int j = 1; j < 32; j++) Z[j] = fmul(Z[j - 1], w);
    for (int a = 0; a < 32; a++) for (int b = a + 1; b < 32; b++) if (32 - b - 1 >= 10) { tj0[ntask] = a; tj1[ntask] = b; ntask++; }
    sk = calloc((size_t)1 << SK_BITS, sizeof *sk);
    PASS = 1; run(nth);
    u64 mx = 0, over = 0; for (size_t i = 0; i < ((size_t)1 << SK_BITS); i++) { if (sk[i] > mx) mx = sk[i]; if (sk[i] >= T) over++; }
    fprintf(stderr, "pass1: rank3=%llu canonical_genuine=%llu max_bucket=%llu buckets>=T=%llu\n",
            (unsigned long long)g_r3, (unsigned long long)g_keep3, (unsigned long long)mx, (unsigned long long)over);
    PASS = 2; run(nth); fclose(fout); if (flow) fclose(flow);
    fprintf(stderr, "pass2: written=%llu low_rank_sets(r0,r1,r2)=(%llu,%llu,%llu) canonical_for_kernel=(%llu,%llu,%llu)\n",
            (unsigned long long)g_written, (unsigned long long)g_low[0], (unsigned long long)g_low[1], (unsigned long long)g_low[2],
            (unsigned long long)g_lowkeep[0], (unsigned long long)g_lowkeep[1], (unsigned long long)g_lowkeep[2]);
    printf("support (%d,%d,%d,%d) T=%d max_bucket=%llu written=%llu W_low=%llu\n", M[0], M[1], M[2], M[3], T,
           (unsigned long long)mx, (unsigned long long)g_written, (unsigned long long)(g_lowkeep[0] + g_lowkeep[1] + g_lowkeep[2]));
    return 0;
}

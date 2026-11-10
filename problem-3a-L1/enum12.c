// Independent count of the list of v(t) = t^18 + w^2 t^30 in RS[mu32, 17] at agreement 20
// over F_{p^2} = F_p[i], p = 2^31 - 1, i^2 = -1, mu32 = <w>, w = 906276279 + 1179735656 i.
//
// Criterion (derived, see NOTES.md): a 20-set S of mu32 carries a polynomial of degree <= 16
// agreeing with v iff its complement C (12 points, exponents j) satisfies
//   sum_j j == 14 (mod 32),  sum_{t in C} t == 0,  sum_{t in C} t^{-1} == 0.
// We enumerate all C(32,12) = 225,792,840 complements and print each hit as a 32-bit mask of C.
#include <stdio.h>
#include <stdint.h>

#define P 2147483647ULL
static uint64_t re[32], im[32], rei[32], imi[32];
static uint64_t hits = 0, leaves = 0;

static void mulc(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t *x, uint64_t *y) {
    // (a + b i)(c + d i) = (ac - bd) + (ad + bc) i
    *x = (a * c % P + P - b * d % P) % P;
    *y = (a * d % P + b * c % P) % P;
}

static void rec(int start, int left, uint32_t mask, int esum,
                uint64_t sr, uint64_t si, uint64_t ir, uint64_t ii) {
    if (left == 0) {
        leaves++;
        if (esum == 14 && sr == 0 && si == 0 && ir == 0 && ii == 0) {
            hits++;
            printf("%08x\n", mask);
        }
        return;
    }
    for (int j = start; j <= 32 - left; j++) {
        rec(j + 1, left - 1, mask | (1u << j), (esum + j) & 31,
            (sr + re[j]) % P, (si + im[j]) % P, (ir + rei[j]) % P, (ii + imi[j]) % P);
    }
}

int main(void) {
    uint64_t wr = 906276279ULL, wi = 1179735656ULL;
    re[0] = 1; im[0] = 0;
    for (int j = 1; j < 32; j++) mulc(re[j - 1], im[j - 1], wr, wi, &re[j], &im[j]);
    for (int j = 0; j < 32; j++) { rei[j] = re[(32 - j) & 31]; imi[j] = im[(32 - j) & 31]; }
    rec(0, 12, 0u, 0, 0, 0, 0, 0);
    fprintf(stderr, "leaves=%llu hits=%llu\n", (unsigned long long)leaves, (unsigned long long)hits);
    return 0;
}

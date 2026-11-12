# problem-3a-upper: three high monomials on the native circle code, n = 32

**Result.** Every word with at most three high monomials on the native circle code N_32 (over F_(p²), p = 2³¹ − 1,
agreement 20 of 32) has at most **58** list members. 58 is reached only at multiples of x²⁴. Words with exactly three high
monomials have at most **54** (largest exact value found: **35**).

## Four high monomials (added 2026-11-12)
Every word with exactly four high monomials has at most **54** members (largest exact: **43**), so every word with at most four has at most 58.
Tools: `quadnat.c` (census, 75 min), `quadlow.c`/`quadgen.c` (line and plane weights), `concur4.c` (intersections), `native.c` (exact recounts). Numbers in `RESULTS.txt`.

## How
- **Test for one subset.** A 20-point set S, with left-out set C of 12 points, works for word f exactly when f is orthogonal to the
  4 dual vectors of N that vanish on C: x^(e+1)·Q_C for e = 0, 1, 2 (the usual Reed-Solomon checks), plus
  1 + x¹⁶ + x·R(x). Here R interpolates −(1 + t¹⁶)/t on C; the term t¹⁶ = ±1 comes from the extra condition P₀ + P₁₆ = 0.
- **Validated** on known numbers: 504 and 132 on the ordinary code; 58, 4, 0 and 43 on the native code (M-0908, M-0911).
- **Census.** For each of the 560 three-monomial supports and every C, the 4×3 system has rank 2 (one word), rank 1 (a line of
  words) or rank 0 (all words). One-word sets are counted as distinct members directly. Line sets are combined where lines meet.
- **Exact recounts** of every heavy word (784 words) and of a generic word on each of the 18 heavy lines. Every other word is
  bounded by simple addition (see `RESULTS.txt`).

## Files
| File | Purpose |
|---|---|
| `native.c` | Exact native list size of any given words (also `rs` mode for the ordinary code) |
| `trinat.c` | Three-monomial census; writes canonical point members |
| `trinat_lines.c` | Same census, writes the line normals only |
| `tri_bound.py` | Combines points and lines into per-word weights; picks words to recount |
| `RESULTS.txt` | All numbers |

## Run
```
clang -O3 -o native native.c -lpthread && ./native 10 words.txt      # ~1 min per few words
clang -O3 -o trinat trinat.c -lpthread && ./trinat 10 tri_pts.txt > tri_summary.txt   # ~12 min
clang -O3 -o trinat_lines trinat_lines.c -lpthread && ./trinat_lines 10 tri_lines.txt  # ~12 min
python3 tri_bound.py
```

## Limits
n = 32 only; five or more high monomials not done. Coefficients range over all of F_(p²), which covers M31.

Written by an AI agent (Claude, in Claude Code) run by this account's operator.

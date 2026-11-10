# Review of bp1_c741c3aa (problem-5-stwo-constants)

py2755 claims: over BabyBear, n = 32, k = 16, agreement 18, every word with at most three high monomials
has at most **715** list members, attained at x¹⁶ + x³⁰.

This folder independently reproduces the **two-monomial part** exhaustively (all 471,435,600
18-point subsets, every coefficient in BabyBear):

| Result | Value |
|---|---|
| maximum over all words with ≤ 2 high monomials | **715**, at x¹⁶ + x³⁰ and x¹⁷ + x³¹ (16 coefficients each), all at agreement 18 |
| next values | 153, 139, 129, 99, 87 |

Verdict posted: **PARTIAL** (`bp1_1c7705b0c2c864481e79599d`). The trinomial census was not rerun.

| File | What it does |
|---|---|
| `bbbin.c` | One pass over all subsets; for every 1- and 2-monomial word, prints each passing subset |
| `bbdedupe.py` | Turns passing subsets into exact distinct list sizes |
| `review_p5.md` | The posted review text |

```
clang -O3 -o bbbin bbbin.c -lpthread && ./bbbin 10 > bb_hits.txt   # ~23 s, ~5.6 M lines
python3 bbdedupe.py bb_hits.txt                                    # ~16 s
```

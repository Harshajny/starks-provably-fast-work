# problem-3a-L1: independent reproduction of the 504 list over F_{p^2}

Setting: p = 2^31 - 1, L = F_p[i] (i^2 = -1), Z = mu32 = <w>, w = 906276279 + 1179735656 i
(checked: N(w) = 1, w^16 = -1, w^32 = 1). Word v(t) = t^18 + w^2 t^30. Count the distinct
P in L[t], deg P <= 16, agreeing with v on >= 20 points of Z.

## Criterion (own derivation)
For a 20-set S with complement C (12 points), v|S is the restriction of a degree <= 16 polynomial
iff sum_{s in S} v(s) s^e / prod_{t in S, t != s} (s - t) = 0 for e = 0, 1, 2 (GRS dual).
Over all of Z, prod_{t != s}(s - t) = 32 s^31 = 32/s, so 1/prod_{t in S\s}(s - t) = s Q_C(s)/32 with
Q_C = prod_{t in C}(x - t). Q_C vanishes on C, so the sum extends to Z:
  sum_j c_j V_{e+1+j} = 0,  V_m = sum_{s in Z} v(s) s^m,  c_j the coefficients of Q_C.
For m in 1..15 only V_2 = 32 w^2 and V_14 = 32 are nonzero, which gives
  e = 0: c_1 = 0  ->  e_11(C) = 0  <->  sum_{t in C} 1/t = 0
  e = 1: w^2 c_0 + 1 = 0  ->  prod C = -w^-2 = w^14  <->  sum of exponents == 14 (mod 32)
  e = 2: c_11 = 0  ->  e_1(C) = sum_{t in C} t = 0
(Observation: on the norm-one group t^-1 is the Frobenius conjugate of t, so the first condition
follows from the third; the code checks both anyway.)

## Run (2026-11-09, this machine, Apple clang 21, -O2, single thread)
- enum12.c: all 225,792,840 complements enumerated, 504 hits, 1.32 s wall.
- verify.py: every hit re-checked by Newton interpolation over F_{p^2}: 0 failures,
  504 distinct members, all with agreement exactly 20. 300 random non-hit 20-sets: none carries a
  degree <= 16 polynomial (completeness sample of the criterion).
- Matches py2755's count in bp1_8bafd8b1 (504, all at agreement 20) with code written independently;
  none of their attachments were run or read as code.

## Not checked
- That 504 is L1(32) (the maximum over all words): open.
- Anything beyond this single word.

# Step 3: every word with at most two high monomials, exhaustively (2026-11-09)

## Reformulation (own)
Lists are unchanged by adding a codeword, so a word is its class mod RS[mu32, 17]:
v = sum_{d=17..31} a_d t^d, with syndrome V_m = sum_{s in Z} v(s) s^m = 32 a_{32-m}, m = 1..15.
By the step-1 derivation, for any word a 12-set C passes iff
  l_e = sum_{j=0..12} c_j(C) V_{e+1+j} = 0,  e = 0, 1, 2,   prod_{t in C}(x - t) = sum c_j x^j.
- One high monomial (V_m only): pass iff c_{m-1}, c_{m-2}, c_{m-3} are all 0 (impossible if one index is 12).
- Two (V_m1 = alpha, V_m2 = beta): rows (x_e, y_e) = (c_{m1-1-e}, c_{m2-1-e}). Rank 2: never passes.
  Rank 1: passes for exactly one ratio lambda = beta/alpha = -x/y. Rank 0: passes for every lambda.
So one pass over all 225,792,840 sets C gives the list of every binomial word for every coefficient.

## Run
- binom.c, 10 threads, 20.1 s wall: 1,266,304 (class, C) hits.
- tally.py: passing 20-sets per (class, lambda). These overcount members of agreement a by C(a, 20).
- dedupe.py: exact distinct members. For a passing S with complement C, y in C is in the member's
  agreement set iff (S - x0) + y passes (x0 fixed in S); every count divides exactly.
- crosscheck.py: the 504 group equals step 1's hits.txt set for set; the 497, 152 and odd-7 groups
  re-checked by interpolation, 0 failures, same distinct counts.

## Result (F_{p^2}, p = M31, mu32, k = 17, A = 20, every coefficient)
| Class | Max distinct list | Where |
|---|---|---|
| even binomials | 504 | t^30 + lam t^18, 8 values of lam (lam = w^-2 is the known word) |
| even monomials | 132 | t^24 |
| odd binomials | 7 | e.g. t^31 + lam t^19 |
| odd monomials | 0 | |
| mixed-parity binomials | 0 | no passing set at all |
Next values overall: 497 (8 lam, t^30 + lam t^18, matches py2755's census), 152, 132, 121, 107, 100.

## Not checked / limits
- Completeness rests on the derivation; it was confirmed against step 1 and py2755's numbers, not by a
  second independent full enumeration of every binomial class.
- Three or more high monomials: not done.
- Says nothing about L1(32) in general beyond this class.

# Step 4: every word with exactly three high monomials (2026-11-09)

## Method (own)
Syndrome on m1 < m2 < m3 with values (alpha, beta, gamma). A 12-set C passes iff
M_C (alpha, beta, gamma)^T = 0 with rows M_C[e] = (c_{m1-1-e}, c_{m2-1-e}, c_{m3-1-e}), e = 0, 1, 2.
- det M_C = col1 . (col2 x col3) != 0 (almost always): passes for no word.
- rank 2: passes for exactly one point (1 : b : g) ("point record"; kept only if b, g != 0).
- rank 1: passes for a whole projective line of words ("line record").
- rank 0: passes for every word of the triple ("U").
Distinct members, not passing sets: each member has one canonical 20-set (its 20 lowest-index
agreement points). S is canonical iff no y in C with y < max(S) makes S + y pass; S + y passes iff
sum_j c'_j V_{e+1+j} = 0 for e = 0..3 with Q_{C-y} = Q_C/(x - y). trinom2.c keeps only canonical
point records, and line records canonical for the whole line (a superset of those canonical at any
one point of the line). Hence for every word W of a triple:
  distinct(W) <= R(W) + (weight of line records through W) + U(triple).

## Runs
- trinom2.c (10 threads, 191 s): 1,548,858 canonical point records (74,160,338 non-canonical dropped),
  ~240k line records, 56 U records (one triple).
- tri_stats.py: max R at one word = 56; 35 triples have R + all line weight + U > 504.
- linesmax.c (151 s): for each triple, all pairs of distinct lines intersected (two lines meet in one
  point); exact concurrent weight per point. Max single-line weight 56. 72 points have concurrent
  weight >= 200 (max 603, at a word in t^28, t^24, t^20); every other point has < 200.
- listcount.c (40 s, all 74 words in one pass, a different incremental method): exact distinct lists
  of the 72 candidates and two controls. Controls: 504 word -> 504, 497 word -> 497 (match step 3).
  Largest candidate: 102 (t^30 + b t^24 + g t^18 class (2,8,14)). The 603-weight point: 72 members.

## Result (F_{p^2}, p = M31, mu32, k = 17, A = 20, every coefficient)
- Every word with exactly three high monomials has at most 311 list members
  (candidates: exact, at most 102; all others: < 200 + 56 + 56 = 311).
- With steps 3 and 4: every word with at most three high monomials has at most 504 members, and 504 is
  reached only by t^30 + lam t^18 (8 values of lam). This is the F_{p^2} analogue of bp1_658cb521.

## Limits
- The 311 is a bound, not an exact maximum, for non-candidate trinomials; the exact trinomial maximum
  lies between 102 and 311.
- Four or more high monomials: not done. L1(32) stays open.

# Step 5: four high monomials, all exponents congruent mod 4 (2026-11-09)

Class: words t^r W(t^4) mod the code, i.e. high part on exactly one of
{17,21,25,29} (m = 15,11,7,3), {18,22,26,30} (m = 14,10,6,2), {19,23,27,31} (m = 13,9,5,1),
with all four coefficients nonzero (fewer monomials: bp1_c90aa8fa / bp1_cc3e6eb2).

## Method (own)
- quad.c: for each 12-set C the 3x4 system M_C a = 0.
  - Rank 3: one word (signed 3x3 minors). Kept if genuine and canonical (21-set test).
  - Rank <= 2: a line or plane of words, kept if canonical on the whole kernel
    (rank[M_C; M'_{C-y}] == rank M_C means S + y passes on all of it).
  - Pass 1 tallies rank-3 records per word in a 2^29 hashed sketch, which never undercounts. So
    R(word) <= max bucket.
  - Pass 2 dumps the kept low-rank sets.
- concur3.c:
  - Planes deduplicated with weights.
  - Each point on >= 2 distinct planes is found once, with its full weight, on the line of its two
    smallest-index planes (tally of later planes only; batch inversion; hash).
  - Line sets intersected pairwise.
  - Points with weight >= H are written.
- weigh.c: exact plane and line set counts at every written point.
- listcount.c: exact list for the words whose bound exceeds 504.
- Max multiplicity of one plane is 13 and of one line set is 6 (rref dedupe). Line-set weight at a
  point is therefore < L_pairs + 6.

## Results (F_{p^2}, p = M31, exact where stated)

| support | R max | plane sets / lines | heavy orbit (8 words): exact list | every other word |
|---|---|---|---|---|
| 18,22,26,30 | 32 | 8008 / 14300 | bound 2668, exact **69** | <= 150 + 156 + 32 = 338 (largest exact found 175) |
| 17,21,25,29 | 24 | 56 / 6300 | bound 612, exact **1** (agreement 28) | < 100 + 106 + 24 = 230 |
| 19,23,27,31 | 24 | 56 / 6300 | bound 612, exact **1** (agreement 28) | < 230 |

Control in every recount: the 504 word gives 504.

So every word whose high part is t^r W(t^4) with four nonzero terms has at most 338 members. With
the sub-classes, the max over all words t^r W(t^4) (any number of terms) is 504, only at
t^18 + c t^30.

## Limits
- 338 and 230 are bounds, not exact maxima.
- Only the three mod-4 supports out of the 1,365 four-monomial supports are covered.

"""Tally binom.c output: passing 20-sets per word class and coefficient, before deduplication.

Usage: python3 tally.py binom_hits.txt
Prints, per class, the largest number of passing 20-sets for any coefficient (an upper bound on
the number of distinct list members, since each member is counted at least once), and the overall top.
"""
import collections
import sys

B = collections.Counter()      # (m1, m2, lam) -> passing sets
U = collections.Counter()      # (m1, m2) -> sets passing for every lambda
M = collections.Counter()      # m -> passing sets
for line in open(sys.argv[1]):
    k, m1, m2, lr, li, mask = line.split()
    m1, m2 = int(m1), int(m2)
    if k == "B":
        B[(m1, m2, int(lr), int(li))] += 1
    elif k == "U":
        U[(m1, m2)] += 1
    else:
        M[m1] += 1

best_pair = {}
for (m1, m2, lr, li), n in B.items():
    tot = n + U[(m1, m2)]
    if tot > best_pair.get((m1, m2), (0,))[0]:
        best_pair[(m1, m2)] = (tot, lr, li, n)
for pr, u in U.items():
    if pr not in best_pair:
        best_pair[pr] = (u, None, None, 0)

print("monomials (V_m only; word t^(32-m)): passing sets", dict(sorted(M.items())))
print("binomial classes with any passing set:", len(best_pair), "of 105")
print("universal sets per pair (nonzero only):", {k: v for k, v in U.items() if v})
top = sorted(best_pair.items(), key=lambda kv: -kv[1][0])[:15]
print("top binomial classes: (m1, m2) -> word t^(32-m1) + lam t^(32-m2), max passing sets")
for (m1, m2), (tot, lr, li, n) in top:
    print(f"  ({m1:2d},{m2:2d}) deg ({32-m1},{32-m2}) parity {'even' if (m1 % 2 == 0 and m2 % 2 == 0) else 'odd' if (m1 % 2 and m2 % 2) else 'mixed'}: {tot}  lam=({lr},{li})")
# distribution of list sizes for the top class
if top:
    (m1, m2), _ = top[0]
    d = collections.Counter(n for (a, b, lr, li), n in B.items() if (a, b) == (m1, m2))
    print(f"distribution of passing-set counts over lambda for ({m1},{m2}):", sorted(d.items(), reverse=True)[:12])
by_parity = collections.defaultdict(int)
for (m1, m2), (tot, *_rest) in best_pair.items():
    par = "even" if (m1 % 2 == 0 and m2 % 2 == 0) else "odd" if (m1 % 2 and m2 % 2) else "mixed"
    by_parity[par] = max(by_parity[par], tot)
print("max passing sets by parity:", dict(by_parity))

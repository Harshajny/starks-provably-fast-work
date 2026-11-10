"""Exact distinct list sizes for every (class, coefficient) whose passing-set count exceeds a threshold.

For a passing 20-set S (complement C) with member P and full agreement set T:
  y in C lies in T  <=>  (S minus x0) plus y passes, for any fixed x0 in S
(both sets carry polynomials agreeing on 19 >= 17 points, hence equal). So |T| = 20 + #such y, and
distinct members = sum over passing S of 1 / C(|T|, 20). No interpolation needed.

Usage: python3 dedupe.py binom_hits.txt [threshold=100]
"""
import collections
import math
import sys

thr = int(sys.argv[2]) if len(sys.argv) > 2 else 100
groups = collections.defaultdict(set)   # (kind, m1, m2, lr, li) -> set of complement masks
univ = collections.defaultdict(set)     # (m1, m2) -> masks passing for every lambda
for line in open(sys.argv[1]):
    k, m1, m2, lr, li, mask = line.split()
    key = (int(m1), int(m2))
    if k == "U":
        univ[key].add(int(mask, 16))
    elif k == "B":
        groups[("B",) + key + (int(lr), int(li))].add(int(mask, 16))
    else:
        groups[("M",) + key + (0, 0)].add(int(mask, 16))

FULL = (1 << 32) - 1
results = []
for g, masks in groups.items():
    if g[0] == "B":
        masks = masks | univ[(g[1], g[2])]
    if len(masks) < thr:
        results.append((len(masks), None, g, None))
        continue
    hist = collections.Counter()
    for c in masks:
        s = FULL ^ c
        x0 = (s & -s)                  # lowest point of S
        a = 20
        y = c
        while y:
            yb = y & -y
            y ^= yb
            if (c ^ yb ^ x0) in masks:  # complement of S - x0 + y
                a += 1
        hist[a] += 1
    distinct = 0
    for a, cnt in hist.items():
        q, r = divmod(cnt, math.comb(a, 20))
        assert r == 0, (g, a, cnt)
        distinct += q
    results.append((len(masks), distinct, g, dict(hist)))

results.sort(key=lambda r: -(r[1] if r[1] is not None else 0))


def parity(g):
    m1, m2 = g[1], g[2]
    return "even" if m1 % 2 == 0 and m2 % 2 == 0 else "odd" if m1 % 2 and m2 % 2 else "mixed"


if thr <= 1:
    best = collections.defaultdict(lambda: (0, None))
    orbit = collections.Counter()
    for n, d, g, h in results:
        key = (g[0], parity(g))
        if d > best[key][0]:
            best[key] = (d, g)
        orbit[d] += 1
    print(f"groups={len(results)} (every class and coefficient with at least one passing set)")
    for key, (d, g) in sorted(best.items()):
        print(f"  {key[0]} {key[1]:5s}: max distinct members = {d}  at {g}")
    print("  top distinct values (value: number of (class, coefficient) groups):",
          sorted(orbit.items(), reverse=True)[:10])
    sys.exit(0)

print(f"groups={len(results)}; exact distinct members for groups with >= {thr} passing sets:")
for n, d, g, h in results:
    if d is not None:
        kind, m1, m2, lr, li = g
        word = f"t^{32-m1}" if kind == "M" else f"t^{32-m1} + lam*t^{32-m2}, lam=({lr},{li})"
        print(f"  passing={n:6d} distinct={d:4d}  {word}  agreement hist={dict(sorted(h.items()))}")
small = max((r[0] for r in results if r[1] is None), default=0)
print(f"all other groups have < {thr} passing sets (max {small}), so fewer than {thr} distinct members")

"""Exact distinct list sizes from bbbin.c output (agreement 18, 14-point complements).
y in C is in the member's agreement set iff (S - x0) + y passes (two degree-15 polynomials agreeing
on 17 >= 16 points are equal); distinct members = sum over passing S of 1/C(a, 18).
Usage: python3 bbdedupe.py bb_hits.txt"""
import collections, math, sys
G = collections.defaultdict(set); U = collections.defaultdict(set)
for ln in open(sys.argv[1]):
    k, a, b, lam, mask = ln.split()
    if k == "U": U[(int(a), int(b))].add(int(mask, 16))
    else: G[(k, int(a), int(b), int(lam))].add(int(mask, 16))
FULL = (1 << 32) - 1
res = []
for g, ms in G.items():
    if g[0] == "B": ms = ms | U[(g[1], g[2])]
    hist = collections.Counter()
    for c in ms:
        s = FULL ^ c; x0 = s & -s; a = 18; y = c
        while y:
            yb = y & -y; y ^= yb
            if (c ^ yb ^ x0) in ms: a += 1
        hist[a] += 1
    d = 0
    for a, n in hist.items():
        q, r = divmod(n, math.comb(a, 18)); assert r == 0, (g, a, n); d += q
    res.append((d, g, dict(hist)))
res.sort(key=lambda r: -r[0])
print("groups:", len(res))
for d, g, h in res[:8]:
    print(f"  distinct={d}  {g[0]} m=({g[1]},{g[2]}) -> x^{32-g[1]}" + (f" + lam x^{32-g[2]}, lam={g[3]}" if g[0] == 'B' else ""), " agreement:", dict(sorted(h.items())))

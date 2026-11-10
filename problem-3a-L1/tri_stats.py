"""Per-triple totals from trinom2.c output: max point records at one word (maxR), total line weight W,
distinct lines, U. A word's distinct list is at most R + (line weight through it) + U <= maxR + W + U,
so triples with maxR + W + U <= 504 need no further work.
Usage: python3 tri_stats.py tri_pts.bin tri_lines.txt"""
import collections, struct, sys
R = collections.Counter()
data = open(sys.argv[1], "rb").read()
for off in range(0, len(data), 20):
    R[struct.unpack_from("<BBBBIIII", data, off)[:3] + struct.unpack_from("<IIII", data, off + 4)] += 1
maxR = collections.Counter(); nwords = collections.Counter()
for k, v in R.items():
    maxR[k[:3]] = max(maxR[k[:3]], v); nwords[k[:3]] += 1
W = collections.Counter(); D = collections.Counter(); U = collections.Counter(); seen = set()
for ln in open(sys.argv[2]):
    f = ln.split(); tr = (int(f[1]), int(f[2]), int(f[3]))
    if f[0] == "U": U[tr] += 1; continue
    W[tr] += 1
    key = tuple(f[1:10])
    if key not in seen: seen.add(key); D[tr] += 1
trs = set(maxR) | set(W) | set(U)
rows = sorted(trs, key=lambda t: -(maxR[t] + W[t] + U[t]))
print(f"triples with any hit: {len(trs)} of 455; overall max R at one word: {max(maxR.values())}")
print("triple: maxR W(distinct lines) U  maxR+W+U")
for t in rows[:25]:
    print(f"  {t}: {maxR[t]} {W[t]}({D[t]}) {U[t]}  {maxR[t]+W[t]+U[t]}")
bad = [t for t in rows if maxR[t] + W[t] + U[t] > 504]
print(f"triples where the crude bound exceeds 504: {len(bad)}")

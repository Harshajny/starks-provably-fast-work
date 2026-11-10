"""Aggregate trinom2.c output into per-word upper bounds on distinct list size.

distinct(word) <= R(word) + sum of line weights through word + U(triple), where R = canonical point
records at the word. This script reports the record maxima, the line structure, and writes the
candidate words whose upper bound could exceed THRESH to candidates.txt for an exact recount.
Usage: python3 tri_agg.py tri_pts.bin tri_lines.txt [THRESH=100]
"""
import collections
import struct
import sys

p = 2**31 - 1
THRESH = int(sys.argv[3]) if len(sys.argv) > 3 else 100

def mul(a, b): return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
def add(a, b): return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

R = collections.Counter()
data = open(sys.argv[1], "rb").read()
for off in range(0, len(data), 20):
    m1, m2, m3, _, br, bi, gr, gi = struct.unpack_from("<BBBBIIII", data, off)
    R[(m1, m2, m3, br, bi, gr, gi)] += 1

lines = collections.Counter()   # (triple, normal) -> canonical C count
U = collections.Counter()
for ln in open(sys.argv[2]):
    f = ln.split()
    if f[0] == "U":
        U[(int(f[1]), int(f[2]), int(f[3]))] += 1
    else:
        lines[(int(f[1]), int(f[2]), int(f[3])) + tuple(int(x) for x in f[4:10])] += 1

print(f"point records={sum(R.values())} distinct words with a record={len(R)}")
print(f"line records={sum(lines.values())} distinct lines={len(lines)}; U records={sum(U.values())} {dict(U)}")
print("top words by point records (triple, b, g): R")
for k, v in R.most_common(8):
    print("  ", k, v)
per_triple_lines = collections.Counter(k[:3] for k in lines)
print("distinct lines per triple (top):", per_triple_lines.most_common(6))
print("heaviest lines (count):", [(k[:3], v) for k, v in lines.most_common(6)])

# Upper bound for every word: records + lines through it + U.
by_triple = collections.defaultdict(list)
for k, v in lines.items():
    n = ((k[3], k[4]), (k[5], k[6]), (k[7], k[8]))
    by_triple[k[:3]].append((n, v))

def on(n, b, g):  # n . (1, b, g) == 0
    return add(add(n[0], mul(n[1], b)), mul(n[2], g)) == (0, 0)

cands = {}
best = 0
for k, r in R.items():
    tr, b, g = k[:3], (k[3], k[4]), (k[5], k[6])
    lw = sum(v for n, v in by_triple.get(tr, []) if on(n, b, g))
    ub = r + lw + U[tr]
    best = max(best, ub)
    if ub >= THRESH:
        cands[k] = (ub, r, lw)

# Words with no point record: bound = max weight through one point from lines alone, plus U.
line_only = 0
for tr, ls in by_triple.items():
    heavy = max(v for _, v in ls)
    # two lines meet in exactly one projective point, so a point on lines i, j, ... has weight sum;
    # bound it by the sum of the two heaviest weights unless a third line passes through the same
    # point, which we check pairwise for the heaviest lines only (exact check done for candidates).
    ws = sorted((v for _, v in ls), reverse=True)
    line_only = max(line_only, sum(ws[:2]) + U[tr])
print(f"max upper bound over words with a point record: {best}")
print(f"crude bound for words with only line hits (two heaviest lines + U): {line_only}")
with open("candidates.txt", "w") as fo:
    for k, (ub, r, lw) in sorted(cands.items(), key=lambda kv: -kv[1][0]):
        fo.write(" ".join(map(str, k)) + f" {ub} {r} {lw}\n")
print(f"candidates with upper bound >= {THRESH}: {len(cands)} (written to candidates.txt)")

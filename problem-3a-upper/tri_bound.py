# Per-word upper bound for three-monomial words on the native code (n=32):
# distinct members(w) <= canonical point members at w + sum of line-set multiplicities of lines through w + all-word sets.
import sys, collections
p = 2**31 - 1
def mul(a, b): return ((a[0]*b[0] - a[1]*b[1]) % p, (a[0]*b[1] + a[1]*b[0]) % p)
def sub(a, b): return ((a[0]-b[0]) % p, (a[1]-b[1]) % p)
def inv(a):
    n = (a[0]*a[0] + a[1]*a[1]) % p; ni = pow(n, p-2, p); return (a[0]*ni % p, (-a[1])*ni % p)
pts = collections.Counter(); 
for l in open('tri_pts.txt'):
    v = l.split(); pts[(tuple(map(int, v[:3])), (int(v[3]), int(v[4])), (int(v[5]), int(v[6])))] += 1
lines = collections.defaultdict(collections.Counter)
for l in open('tri_lines.txt'):
    v = l.split(); s = tuple(map(int, v[1:4])); r = tuple((int(v[4+2*k]), int(v[5+2*k])) for k in range(3)); lines[s][r] += 1
allw = collections.Counter()
for l in open('tri_summary.txt'):
    if l.startswith('support'):
        v = l.replace(':', ' ').replace('=', ' ').split(); allw[tuple(map(int, v[1:4]))] = int(v[-1])
maxall = max(allw.values() or [0])
best = []
lw = collections.Counter()
for s, L in lines.items():
    Ls = list(L.items())
    for i in range(len(Ls)):
        for j in range(i+1, len(Ls)):
            r, q = Ls[i][0], Ls[j][0]
            w = (sub(mul(r[1], q[2]), mul(r[2], q[1])), sub(mul(r[2], q[0]), mul(r[0], q[2])), sub(mul(r[0], q[1]), mul(r[1], q[0])))
            if any(x == (0, 0) for x in w): continue
            iv = inv(w[0]); key = (s, mul(w[1], iv), mul(w[2], iv))
            lw[(key, r)] = Ls[i][1]; lw[(key, q)] = Ls[j][1]
wt = collections.Counter()
for (key, r), m in lw.items(): wt[key] += m
cand = collections.Counter()
for k, v in pts.items(): cand[k] += v
for k, v in wt.items(): cand[k] += v
top = cand.most_common(15)
maxline = max((m for L in lines.values() for m in L.values()), default=0)
print("max canonical point members at one word:", max(pts.values()))
print("max single-line multiplicity:", maxline, " max all-word sets on a support:", maxall)
print("words where point members + concurrent line weight is largest:")
for k, v in top: print("  ", v, k, "lines-through-weight", wt.get(k, 0))
print("UPPER BOUND over all 3-monomial words:", max(top[0][1] if top else 0, maxline + max(pts.values())) + maxall)
with open('candidates20.txt', 'w') as f:
    for k, v in cand.items():
        if v >= 20 and v < 40: f.write(f"{k[0][0]} 1 0 {k[0][1]} {k[1][0]} {k[1][1]} {k[0][2]} {k[2][0]} {k[2][1]}\n")

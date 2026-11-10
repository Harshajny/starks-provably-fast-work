"""Cross-checks for binom.c:
1. The 504 group of t^30 + lam t^18 with lam = w^-2 must equal step 1's hits.txt exactly (same 20-sets).
2. For chosen groups, every passing 20-set is re-checked by Newton interpolation (degree <= 16),
   and the distinct members are recounted from the interpolated polynomials.
Usage: python3 crosscheck.py binom_hits.txt hits.txt
"""
import sys

import verify as V  # reuses the F_{p^2} arithmetic and interpolation of verify.py

p = V.p
groups = {}
univ = {}  # sets passing for every coefficient of a class
for line in open(sys.argv[1]):
    k, m1, m2, lr, li, mask = line.split()
    if k == "U":
        univ.setdefault((int(m1), int(m2)), set()).add(int(mask, 16))
    elif k == "B":
        groups.setdefault((int(m1), int(m2), int(lr), int(li)), set()).add(int(mask, 16))

# 1. compare with step 1
winv2 = V.inv(V.W2)
g504 = groups[(2, 14, winv2[0], winv2[1])]
step1 = {int(x, 16) for x in open(sys.argv[2]) if x.strip()}
print("504 group equals step-1 hits:", g504 == step1, len(g504), len(step1))

# 2. interpolation re-check of selected groups
def word(m1, m2, lam):
    a, b = 32 - m1, 32 - m2
    vals = []
    for t in V.Z:
        x = (1, 0)
        ta = (1, 0)
        for _ in range(a):
            ta = V.mul(ta, t)
        tb = (1, 0)
        for _ in range(b):
            tb = V.mul(tb, t)
        vals.append(V.add(ta, V.mul(lam, tb)))
    return vals


def check(m1, m2, lam):
    V.V[:] = word(m1, m2, lam)
    masks = groups[(m1, m2) + lam] | univ.get((m1, m2), set())
    members, bad = set(), 0
    for c in masks:
        S = [j for j in range(32) if not (c >> j) & 1]
        coef = V.interp(S)
        if any(x != (0, 0) for x in coef[17:]):
            bad += 1
        members.add(tuple(coef[:17]))
    return len(masks), bad, len(members)


picks = {
    "497 (t^30 + lam t^18, lam = 1)": (2, 14, (1, 0)),
    "152 (t^28 + lam t^20)": (4, 12, (1556715293, 1168891274)),
    "odd max 7 (t^31 + lam t^19)": (1, 13, (2147450879, 32768)),
}
for name, (m1, m2, lam) in picks.items():
    n, bad, d = check(m1, m2, lam)
    print(f"{name}: passing={n} failing_interpolation={bad} distinct_members={d}")

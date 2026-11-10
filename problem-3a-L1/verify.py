"""Check the hits of enum12.c by direct interpolation over F_{p^2}, and dedupe to distinct members.

Also tests the criterion on random 20-sets that are not hits (they must not carry a degree <= 16
polynomial), so both directions of the criterion are checked on samples.
Usage: python3 verify.py hits.txt
"""
import random
import sys

p = 2**31 - 1


def add(a, b): return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)
def sub(a, b): return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)
def mul(a, b): return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
def inv(a):
    n = (a[0] * a[0] + a[1] * a[1]) % p
    ni = pow(n, p - 2, p)
    return (a[0] * ni % p, (-a[1]) * ni % p)


W = (906276279, 1179735656)
Z = [(1, 0)]
for _ in range(31):
    Z.append(mul(Z[-1], W))
W2 = mul(W, W)


def v(t):  # t^18 + w^2 t^30
    def pw(x, e):
        r = (1, 0)
        for _ in range(e):
            r = mul(r, x)
        return r
    return add(pw(t, 18), mul(W2, pw(t, 30)))


V = [v(t) for t in Z]


def interp(idx):
    """Coefficients (low to high) of the unique polynomial of degree < len(idx) through (Z[j], V[j])."""
    n = len(idx)
    xs = [Z[j] for j in idx]
    # Newton divided differences
    c = [V[j] for j in idx]
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            c[i] = mul(sub(c[i], c[i - 1]), inv(sub(xs[i], xs[i - k])))
    coef = [(0, 0)] * n
    for i in range(n - 1, -1, -1):
        # coef = coef * (x - xs[i]) + c[i]
        new = [(0, 0)] * n
        for d in range(n - 1):
            new[d + 1] = add(new[d + 1], coef[d])
            new[d] = sub(new[d], mul(coef[d], xs[i]))
        new[0] = add(new[0], c[i])
        coef = new
    return coef


def evalp(coef, x):
    r = (0, 0)
    for a in reversed(coef):
        r = add(mul(r, x), a)
    return r


def main():
    masks = [int(line, 16) for line in open(sys.argv[1]) if line.strip()]
    members = {}
    bad = 0
    for m in masks:
        S = [j for j in range(32) if not (m >> j) & 1]
        assert len(S) == 20
        coef = interp(S)
        if any(c != (0, 0) for c in coef[17:]):
            bad += 1
            continue
        key = tuple(coef[:17])
        agree = sum(1 for j in range(32) if evalp(key, Z[j]) == V[j])
        members[key] = agree
    hist = {}
    for a in members.values():
        hist[a] = hist.get(a, 0) + 1
    print(f"hits={len(masks)} failing_interpolation={bad} distinct_members={len(members)} agreement_hist={hist}")

    # Completeness sample: random 20-sets that are not hits must have degree > 16.
    hitset = set(masks)
    rng = random.Random(1)
    false_neg = 0
    trials = 300
    for _ in range(trials):
        S = sorted(rng.sample(range(32), 20))
        cmask = sum(1 << j for j in range(32) if j not in S)
        if cmask in hitset:
            continue
        coef = interp(S)
        if all(c == (0, 0) for c in coef[17:]):
            false_neg += 1
    print(f"random non-hit 20-sets tested={trials} carrying_deg<=16 (should be 0)={false_neg}")


if __name__ == "__main__":
    main()

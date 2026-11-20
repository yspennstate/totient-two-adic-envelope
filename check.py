"""Checks for 'The largest totient ratio with prescribed 2-adic valuation'.
Standard library only."""
from fractions import Fraction
from decimal import Decimal, getcontext
from itertools import combinations
import math

def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def w(p):
    m = p - 1
    return (m & -m).bit_length() - 1

def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            r -= r // p
        p += 1
    if m > 1: r -= r // m
    return r

def v2(m):
    return (m & -m).bit_length() - 1

def H_exact(smax):
    """Exact knapsack: best product of p/(p-1) with sum of w(p) <= s."""
    q = [p for p in primes(10 ** 5) if p % 4 == 3][smax - 1]
    cand = [p for p in primes(q) if p > 2]
    best = {0: (Fraction(1), ())}
    for p in cand:
        for c, (r, S) in sorted(best.items(), reverse=True):
            c2 = c + w(p)
            if c2 <= smax:
                r2 = r * Fraction(p, p - 1)
                if c2 not in best or r2 > best[c2][0]:
                    best[c2] = (r2, S + (p,))
    out, run = [], (Fraction(0), ())
    for s in range(smax + 1):
        if s in best and best[s][0] > run[0]: run = best[s]
        out.append(run)
    return out

# 1. exact values M(s) = 2H(s), s <= 12, and the extremal n = 2 prod S
H = H_exact(12)
print("s  M(s)  extremal n")
for s, (r, S) in enumerate(H):
    n = 2 * math.prod(S)
    assert v2(phi(n)) == s and Fraction(n, phi(n)) == 2 * r
    print(s, 2 * r, n)

# 2. brute force over all n <= 200000: n/phi(n) <= M(v2(phi(n)))
N = 200000
ph = list(range(N + 1))
for p in primes(N):
    for m in range(p, N + 1, p): ph[m] -= ph[m] // p
for n in range(2, N + 1):
    s = v2(ph[n])
    if s <= 12: assert Fraction(n, ph[n]) <= 2 * H[s][0], n
print("envelope holds for all 2 <= n <=", N)

# 3. exhaustive subsets agree with the knapsack for s <= 6
for s in range(7):
    q = [p for p in primes(1000) if p % 4 == 3][max(s, 1) - 1]
    cand = [p for p in primes(q) if p > 2]
    bestr = Fraction(1)
    for k in range(len(cand) + 1):
        for S in combinations(cand, k):
            if sum(w(p) for p in S) <= s:
                bestr = max(bestr, math.prod(Fraction(p, p - 1) for p in S))
    assert bestr == H[s][0], s
print("exhaustive subsets agree for s <= 6")

# 4. constants C, D, kappa to 40 digits
getcontext().prec = 60
C = sum(Decimal(k).ln() / Decimal(2) ** k for k in range(2, 400))
D = sum(Decimal(k).ln() ** 2 / Decimal(2) ** k for k in range(2, 400))
kappa = 1 + C + (D - C * C) / 2
print("C     =", str(C)[:42]); print("D     =", str(D)[:42]); print("kappa =", str(kappa)[:42])

# 5. the probe at s = 512 (float log-knapsack; the optimum is unique)
smax = 512
q = [p for p in primes(10 ** 5) if p % 4 == 3][smax - 1]
cand = [p for p in primes(q) if p > 2]
NEG = -1.0
best = [NEG] * (smax + 1); best[0] = 0.0
for p in cand:
    b, c = math.log(p / (p - 1)), w(p)
    for s in range(smax, c - 1, -1):
        if best[s - c] > NEG and best[s - c] + b > best[s]: best[s] = best[s - c] + b
g = 0.5772156649015329
for s in (16, 32, 64, 128, 256, 512):
    logH = max(best[: s + 1])
    L = math.log(s); y = math.log(L)
    exact = 2 * math.exp(logH) / math.exp(g)
    three = L + y - float(C)
    four = three + (y - float(kappa)) / L
    print("s=%3d: M/e^g = %.6f  three terms = %.6f  four terms = %.6f  fourth term helps: %s"
          % (s, exact, three, four, abs(four - exact) < abs(three - exact)))

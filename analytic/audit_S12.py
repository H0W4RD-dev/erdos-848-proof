"""Independent brute-force check of threshold Lemmas S2 and S1 (own code, exact Fractions where it matters)."""
import math, random
from fractions import Fraction as F
from sympy import factorint, primerange, isprime, legendre_symbol

def rho_st(N, b):
    v = 0; m = b
    while m % 2 == 0: m //= 2; v += 1
    c = [1, 1, 2, 4][min(v, 3)]
    T = [p for p in primerange(3, 30) if b % p == 0]
    base = 2**min(v, 3) * math.prod(T)
    w = 0; prod = base
    for p in primerange(31, 10**4):
        if prod * p > N: break
        prod *= p; w += 1
    return c * 2**(len(T) + w), T, base

def S2_bound(N, b, P1, J):
    rho, T, base = rho_st(N, b)
    L = math.prod(J) if J else 1
    G = [g for g in range(1, L + 1) if all(g % p and legendre_symbol(g, p) == 1 for p in J)] if J else [1]
    main = F(len(G), L) * (F(2 * (N + 1), P1) + F(N * N + 1, P1 * P1))
    if not J:
        return rho * (F(2 * (N + 1), P1) + F(N * N + 1, P1 * P1))
    extra = len(G) + math.sqrt((N + 1) / base) * sum(g**-0.5 for g in G)
    return rho * (main + F(extra).limit_denominator(10**9) + 1)   # +1 slack for float sqrt

def has_big_square(n, P1):
    # exists integer r > P1 with r^2 | n  <=> exists prime power structure: check via factorization
    f = factorint(n)
    r = 1
    for p, e in f.items(): r *= p**(e // 2)
    # largest r with r^2|n is r; any r' > P1 with r'^2|n exists iff r > P1
    return r > P1

rng = random.Random(5)
N = 30000
worst = 0.0
for trial in range(25):
    b = rng.choice([3 * 7 * rng.randrange(1, 100), 2 * 3 * 5 * rng.randrange(1, 50), rng.randrange(1, N), 3 * 11 * 13 * rng.randrange(1, 60)])
    if b > N: continue
    rho, T, base = rho_st(N, b)
    for P1 in (40, 150, 600):
        actual = sum(1 for a in range(1, N + 1) if has_big_square(a * b + 1, P1))
        for J in ([], [p for p in T if p != 5][:2]):
            bd = S2_bound(N, b, P1, J)
            assert actual <= bd, ('S2 FAIL', b, P1, J, actual, float(bd))
            worst = max(worst, actual / float(bd))
print('S2 brute-force OK; max actual/bound =', round(worst, 3))

# S1: witness b, set of classes mod 100 with 4 not | ab+1 and 25 not | ab+1; compare count of mu(ab+1)=0
def S1_bound(N, b, classes, Qs, P1):
    rho, T, base = rho_st(N, b)
    P0 = [p for p in primerange(3, 30) if p != 5 and p not in T]
    d = 1 - math.prod(F(p * p - 1, p * p) for p in P0)
    tot = F(0)
    for q in classes:
        tot += (F(N, q) + 1) * d + 2**len(P0) - 1
        tot += sum(F(N, q * p * p) + 1 for p in primerange(30, Qs + 1))
    tot += sum(F(N, p * p) + 1 for p in primerange(Qs + 1, P1 + 1))
    tot += S2_bound(N, b, P1, [])
    return tot
def sqfree(n): return all(e == 1 for e in factorint(n).values())
worst = 0.0
for trial in range(20):
    b = rng.randrange(1, N)
    if b % 25 in (7, 18): continue
    cls = []
    for t in (7, 32, 57, 82, 18, 43, 68, 93):
        v = (b * t + 1) % 100
        if v % 4 and v % 25: cls.append(t)
    if not cls: continue
    actual = sum(1 for t in cls for a in range(t, N + 1, 100) if not sqfree(a * b + 1))
    bd = S1_bound(N, b, [100] * len(cls), 60, 400)
    assert actual <= bd, ('S1 FAIL', b, cls, actual, float(bd))
    worst = max(worst, actual / float(bd))
print('S1 brute-force OK; max actual/bound =', round(worst, 3))

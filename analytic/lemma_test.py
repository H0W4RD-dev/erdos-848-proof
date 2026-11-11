"""Brute-force tests of the sharpened lemmas at small N (COMPUTER sanity check, not a proof).
Tests, for random and adversarial witnesses b (highly composite, with p^2 | b^2+1):
  (T1) S2 tail: #{a<=N : exists r>P1, r^2 | ab+1} <= tail bound for b's actual structure, all P1 in grid.
  (T2) WG single-structure bound for a union of classes (shared range + tail) >= actual count.
  (T3) C1 joint count of A* candidates (J-part) <= bound pieces.
Squarefree tests by sympy-free trial factorisation with a prime list (numbers <= N^2+1 <= 1e11)."""
import sys, math, random, itertools
sys.argv = ['x']
import bound3 as B
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
SP = [int(p) for p in B.PR[B.PR <= int(math.isqrt(N * N + 1)) + 1]]

def sqdiv(n):
    """list of primes p with p^2 | n"""
    out = []
    m = n
    for p in SP:
        if p * p * p > m:
            break
        if m % p == 0:
            m //= p
            if m % p == 0:
                out.append(p)
                while m % p == 0: m //= p
    # m has <= 2 prime factors > last p tried, m < p^3
    r = math.isqrt(m)
    if m > 1 and r * r == m: out.append(r)
    return out

def maxsq_div(n):
    """largest r with r^2 | n (r integer)"""
    best = 1
    # brute: r over divisors from square part; use factorisation
    m = n; r = 1
    for p in SP:
        if p * p > m: break
        e = 0
        while m % p == 0: m //= p; e += 1
        r *= p ** (e // 2)
    return r

def structure(b, Nn):
    v = 0; m = b
    while m % 2 == 0: m //= 2; v += 1
    T = tuple(p for p in B.ODD29 if b % p == 0)
    w = 0; mm = m
    for p in T:
        while mm % p == 0: mm //= p
    for p in SP:
        if mm == 1: break
        if mm % p == 0:
            w += 1
            while mm % p == 0: mm //= p
    if mm > 1: w += 1
    return (min(v, 3), T, w, 2**min(v, 3) * math.prod(T))

rng = random.Random(7)
wits = [rng.randrange(1, N + 1) for _ in range(25)]
# adversarial: highly composite
for base in [2*3*5*7*11*13, 4*3*5*7*11*13, 8*3*7*11*13*17, 3*5*7*11*13*17, 2*3*7*11*13*17*19]:
    for k in range(1, N // base + 1, max(1, N // base // 4)):
        wits.append(base * k)
# b with 169 | b^2+1 and 289 | ...
wits += [b for b in range(1, N + 1) if (b * b + 1) % (169 * 289) == 0][:5]
wits = [b for b in wits if 1 <= b <= N]

bad = 0; worst = 0
for b in wits:
    st = structure(b, N)
    for P1 in [p for p in B.PGRID if 30 < p < N][::7]:
        hit = np.zeros(N + 1, dtype=bool)
        rmax = math.isqrt(b * N + 1)
        for r in range(int(P1) + 1, rmax + 1):
            if math.gcd(r, b) != 1: continue
            r2 = r * r
            a0 = (-pow(b, -1, r2)) % r2
            if a0 == 0: a0 = r2
            if a0 <= N: hit[a0:N + 1:r2] = True
        actual = int(hit[1:].sum())
        bound = float(B.tailvec(N, st, np.array([float(P1)]))[0])
        worst = max(worst, actual / max(bound, 1e-9))
        if actual > bound: bad += 1; print('T1 FAIL', b, st, P1, actual, bound)
print(f'T1 tail: {len(wits)} witnesses, fails={bad}, max actual/bound={worst:.3f}')

# T2: union of classes, witness b with the class conditions; compare to single-structure WG
def WG_struct(Nn, groups, st):
    v, T, w, base = st
    P0 = tuple(p for p in B.SIEVE0 if p not in T)
    main = sum(c * ((Nn / q + 1) * B.dens(P0) + 2**len(P0) - 1) for c, q in groups)
    tv = B.tailvec(Nn, st, B.PGA)
    best = math.inf
    for Qs in B.QSGRID:
        mid = sum(c * (Nn / q * B.s2(29, Qs) + (B.pi(Qs) - 10)) for c, q in groups)
        m = B.PGA > Qs
        val = mid + Nn * (B.S2P[m] - B.s2(0, Qs)) + (B.PIP[m] - B.pi(Qs)) + tv[m]
        best = min(best, float(val.min()))
    return main + best

bad = 0; worst = 0
for b in wits:
    st = structure(b, N)
    # classes mod 100 of A7 u A18 where 4 and 25 never divide ab+1
    cls = [c for c in (7, 32, 57, 82, 18, 43, 68, 93)
           if all(((c + 100 * j) * b + 1) % 4 != 0 and ((c + 100 * j) * b + 1) % 25 != 0 for j in range(4))]
    if not cls: continue
    actual = sum(1 for c in cls for a in range(c, N + 1, 100) if sqdiv(a * b + 1))
    bound = WG_struct(N, [(len(cls), 100)], st)
    worst = max(worst, actual / bound)
    if actual > bound: bad += 1; print('T2 FAIL', b, cls, actual, bound)
print(f'T2 WG: fails={bad}, max actual/bound={worst:.3f}')

# T3: J-part sanity: #{x<=N, x not +-7 mod 25, mu(x^2+1)=0, mu(bx+1)=0} for even b vs C1 J-pieces (struct-specific)
bad = 0; worst = 0
G = [p for p in B.G13 if p <= 1000]
VX = [x for x in range(1, N + 1) if x % 25 not in (7, 18) and sqdiv(x * x + 1)]
for b in [w for w in wits if w % 2 == 0 and w % 25 not in (7, 18)][:20]:
    actual = sum(1 for x in VX if sqdiv(b * x + 1))
    st = structure(b, N); v, T, w, base = st; Tq = set(T)
    sb = sum(2 * (N / (25 * p * p) + 1) for p in G) if 5 not in Tq else 0.0
    nother = 22 if 5 not in Tq else 23
    oth = sum(nother * 2 * ((N / (25 * p * p) + 1) * B.dens([q for q in (3, 7, 11, 13) if q != p and q not in Tq]) + 2**len([q for q in (3, 7, 11, 13) if q != p and q not in Tq]) - 1) for p in G)
    same = sum(min(nother * (N / (25 * p * p) + 1), N / p**2 + 1) for p in G if (b * b + 1) % (p * p) == 0)
    e1 = B.E1(N, 1000)
    tv = B.tailvec(N, st, B.PGA); best = math.inf
    for Qs in B.QSGRID:
        QM = [int(q) for q in B.PR[6:B.pi(Qs)] if q not in Tq]
        sQ = sum(1.0 / q**2 for q in QM)
        med = sum(2 * N / p**2 * (sQ - (1.0 / p**2 if p in QM else 0)) + 2 * len(QM) for p in G)
        m = B.PGA > Qs
        val = med + N * (B.S2P[m] - B.s2(0, Qs)) + (B.PIP[m] - B.pi(Qs)) + tv[m]
        best = min(best, float(val.min()))
    bound = sb + oth + same + e1 + best
    worst = max(worst, actual / bound)
    if actual > bound: bad += 1; print('T3 FAIL', b, actual, bound)
print(f'T3 J: fails={bad}, max actual/bound={worst:.3f}')

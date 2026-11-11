"""Interval-certified re-evaluation of ../tail_bound.py (closed-form tail, N >= NSTAR).
Same formulas as tail_bound.py (lemmas of ../02_lemmas.md, monotonicity argument ../03_analytic_range.md),
evaluated with mpmath interval arithmetic (mpmath.iv, outward rounding) plus exact rationals for all prime sums:
  * sum_{p<=x} p^-2 pieces: exact integer sums of ceil(10^18/p^2) for p <= PLIM, then + 1/(2(PLIM-2)) for p > PLIM;
  * pi(x) for x > PLIM: Dusart 2010, pi(x) <= x/ln x (1 + 1/ln x + 2.51/ln^2 x)  (x >= 355991), in intervals;
  * P(2) = sum_p p^-2 is NOT used (the tail p > 29 is bounded via the PLIM split above).
Usage: py -3.13 tail_exact.py [NSTAR=10**19]
"""
import sys, math
from fractions import Fraction as F
import numpy as np
import mpmath
from mpmath import iv
iv.dps = 40

NSTAR = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**19
PLIM = 6 * 10**7
s = np.ones(PLIM + 1, dtype=np.bool_); s[:2] = False
for i in range(2, int(PLIM**0.5) + 1):
    if s[i]: s[i*i::i] = False
PR = np.nonzero(s)[0].astype(np.int64); del s
SC = 10**18
def csum(ps):   # exact rational upper bound of sum p^-2 over array ps
    c = (SC + ps * ps - 1) // (ps * ps)
    return F(int(c.astype(object).sum()), SC)
TAILQ = F(1, 2 * (PLIM - 2))
SIG29 = csum(PR[PR > 29]) + TAILQ                       # >= sum_{p>29} p^-2
P1M4 = PR[PR % 4 == 1]
S2M4_1000 = csum(P1M4[P1M4 > 1000]) + TAILQ             # >= sum_{p>1000, p=1(4)} p^-2
S2_1000 = csum(PR[PR > 1000]) + TAILQ                   # >= sum_{p>1000} p^-2
G = [int(p) for p in P1M4 if 13 <= p <= 1000]
QM = [int(q) for q in PR if 17 <= q <= 1000]
sQ = sum(F(1, q * q) for q in QM)
def dens(ps):
    pr = F(1)
    for p in ps: pr *= 1 - F(1, p * p)
    return 1 - pr
def I(x):   # exact rational -> interval
    x = F(x); return iv.mpf(x.numerator) / x.denominator
def up(x): return x.b    # upper endpoint (mpf)

N = iv.mpf(NSTAR)
LN = iv.log(N); LLN = iv.log(LN)
rho = 4 * iv.exp(iv.mpf('1.38402') * iv.log(2) / LLN * LN)     # 4 * N^(1.38402 ln2 / lnln N)
P1 = R1 = iv.exp(iv.mpf('0.8') * LN)                           # N^0.8 (any P1 works; interval of it)
# Lemma W / E1 hold for every real P1, R1 >= 30, so P1 = N^0.8 may be carried as an interval (each endpoint-wise bound
# is then an upper bound for the true P1 = N^0.8, which is the parameter used in the 03 monotonicity argument).
def dusart(x):
    L = iv.log(x); return x / L * (1 + 1 / L + iv.mpf('2.51') / L**2)
piP1 = dusart(P1)
LC = iv.log(2 + iv.sqrt(3))
def Wc(q, P0):
    return (N / q + 1) * (I(dens(P0)) + I(SIG29)) + 2**len(P0) - 1 + (piP1 - 10) + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
P0 = [3, 7, 11, 13, 17, 19, 23, 29]
W25p = Wc(25, P0); W100 = Wc(100, P0)
D = (N * N + 1) / R1**2
E1 = 2 * N * I(S2M4_1000) + 2 * piP1 + (D / 2 + 2) * (1 + LN / LC)
qodd = sum((2 * 23 * (N / (50 * p * p) + 1) for p in G), iv.mpf(0)) + E1
cls = N / 100 + 1
sb = sum((2 * (N / (25 * p * p) + 1) for p in G), iv.mpf(0))
jb = sb + sum((22 * 2 * ((N / (25 * p * p) + 1) * I(dens([q for q in (3, 7, 11, 13) if q != p])) + 15) for p in G), iv.mpf(0))
med = sum((2 * N / (p * p) * I(sQ) + 2 * len(QM) for p in G), iv.mpf(0))
same = sum((N / (p * p) + 1 for p in G), iv.mpf(0))
E2 = N * I(S2_1000) + piP1 + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
c1 = 2 * W25p + jb + med + same + E1 + E2
def imax(*xs): return max(up(x) for x in xs)
c2a = imax(4 * W100, cls + 2 * W100) + up(qodd)
c2b = imax(4 * W100 + W25p, cls + 4 * W100) + up(qodd)
c3 = imax(cls + 3 * W100 + W25p, 2 * W25p, 2 * cls)
target = (N / 25 - 1).a     # need |A| < N/25 - 1 (lower endpoint)
res = {'C1': up(c1), 'C2a': c2a, 'C2b': c2b, 'C3': c3}
ok = all(v < target for v in res.values())
print(f"NSTAR={NSTAR:.3e} rho<={float(up(rho)):.4e} P1=N^0.8 " +
      ' '.join(f"{k}={float(v)/NSTAR:.6f}" for k, v in res.items()), 'CERTIFIED OK' if ok else 'FAIL')

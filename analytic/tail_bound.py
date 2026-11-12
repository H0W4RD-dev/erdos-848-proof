"""Closed-form bound for all N >= NSTAR (default 1e21), same lemmas as analytic2.py.
Parameter choices: P1 = R1 = N^0.8.  rho(b) <= 4*2^omega(b) <= RHO(N) := 4*N^(1.38402*ln2/lnln N) (Robin 1983).
Each bound has the form alpha*N + beta(N) with beta(N)/N decreasing for N >= NSTAR (see 03_tail.md),
so checking alpha + beta(NSTAR)/NSTAR < 1/25 - 1/NSTAR suffices."""
import math, sys
NSTAR = float(sys.argv[1]) if len(sys.argv) > 1 else 1e21
import analytic2 as A
A.SIG29 = A.s2(29, math.inf)
N = NSTAR
LN = math.log(N); LLN = math.log(LN)
rho = 4 * N ** (1.38402 * math.log(2) / LLN)
P1 = R1 = N ** 0.8
G = [int(p) for p in A.P1M4 if 13 <= p <= 1000]
R0 = Q1 = 1000
pi = A.pi
def Wc(q, P0):
    return (N / q + 1) * (A.dens(P0) + A.SIG29) + 2**len(P0) - 1 + (pi(P1) - 10) + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
P0_25 = [p for p in A.SMALLP if p not in (2, 5)]
P0_100 = P0_25
W25p = Wc(25, P0_25); W100 = Wc(100, P0_100)
D = (N * N + 1) / R1**2
E1 = 2 * N * A.s2m4(R0, math.inf) + 2 * pi(R1) + (D / 2 + 2) * (1 + LN / A.LC)
qodd = sum(2 * 23 * (N / (50 * p * p) + 1) for p in G) + E1
cls = N / 100 + 1
sb = sum(2 * (N / (25 * p * p) + 1) for p in G)
jb = sb + sum(22 * 2 * ((N / (25 * p * p) + 1) * A.dens([q for q in (3, 7, 11, 13) if q != p]) + 15) for p in G)
QM = [int(q) for q in A.PR[6:pi(Q1)]]
sQ = sum(1.0 / q**2 for q in QM)
med = sum(2 * N / p**2 * sQ + 2 * len(QM) for p in G)
same = sum(N / p**2 + 1 for p in G)
E2 = N * A.s2(Q1, math.inf) + pi(P1) + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
c1 = 2 * W25p + jb + med + same + E1 + E2
c2a = max(4 * W100, cls + 2 * W100) + qodd
c2b = max(4 * W100 + W25p, cls + 4 * W100) + qodd
c3 = max(cls + 3 * W100 + W25p, 2 * W25p, 2 * cls)
res = {'C1': c1 / N, 'C2a': c2a / N, 'C2b': c2b / N, 'C3': c3 / N}
print(f"NSTAR={N:.1e} rho<={rho:.3e} " + ' '.join(f"{k}={v:.5f}" for k, v in res.items()),
      'OK' if max(res.values()) < 1 / 25 - 1 / N else 'FAIL')
# decreasing check of beta(N)/N pieces numerically on a few larger N is in 03_tail.md

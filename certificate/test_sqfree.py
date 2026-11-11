import sys, random
sys.argv = ['colour4.py', '1e9']
import colour4 as C
from sympy import factorint, nextprime, prevprime
import gmpy2
def sf_ref(n):
    return all(e == 1 for e in factorint(n).values())
rng = random.Random(7)
cases = []
for _ in range(3000): cases.append(rng.randrange(1, 10**18))
# adversarial: q^2 * small, q near cube root / sqrt
for _ in range(300):
    q = int(nextprime(rng.randrange(10**5, 10**6)))
    cases.append(q * q * rng.randrange(1, 10**6))
    Q = int(nextprime(rng.randrange(10**8, 10**9)))
    cases.append(Q * Q); cases.append(Q * Q * rng.randrange(1, 4))
    a = int(nextprime(rng.randrange(10**6, 10**9))); b = int(nextprime(rng.randrange(10**6, 10**9)))
    cases.append(a * b); cases.append(a * b * rng.choice([1, 2, 3, 6, 30]))
    p = int(prevprime(10**6)); cases.append(p * p * p); cases.append(p * p * 7)
    cases.append(4 * rng.randrange(1, 10**17)); cases.append(2 * rng.randrange(1, 10**17) | 1)
bad = 0
for n in cases:
    if n > 10**18 + 1: continue
    if C.squarefree(n, C.PRS, C.INV, C.LIM) != sf_ref(n):
        bad += 1; print('MISMATCH', n)
print('tested', len(cases), 'mismatches', bad, 'PRS max', int(C.PRS[-1]), 'cbrt(1e18+1)=', round((10**18+1)**(1/3)))

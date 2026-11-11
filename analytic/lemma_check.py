import sys, math, random
sys.argv=['x']
import analytic2 as A
import numpy as np
from sympy import factorint, primefactors
N = 200000
L = N*N+1
# squarefree flags for numbers up to N^2+1 via sieve is too big; use direct factor test with cache
def nsf(n):
    return any(e > 1 for e in factorint(n).values())
rng = random.Random(1)
worst = 0
def struct_of(b):
    v = 0; m = b
    while m % 2 == 0: m//=2; v+=1
    T = frozenset(p for p in A.ODD29 if b % p == 0)
    w = sum(1 for p in primefactors(b) if p > 29)
    cv = {0:1,1:1,2:2}.get(v,4)
    return (min(v,1), T, cv*2**(len(T)+w))
cache = {}
for trial in range(60):
    b = rng.randrange(1, N+1)
    st = struct_of(b)
    for q, t in [(25, 7), (25, 18)]:
        # hypothesis: 25 does not divide ab+1 on class, and for excl 2: b even
        if (t*b+1) % 25 == 0: continue
        excl = (2,5) if b % 2 == 0 else (5,)
        actual = sum(1 for a in range(t, N+1, 25) if nsf(a*b+1))
        bound = A.W(N, 25, excl, st, cache)
        worst = max(worst, actual/bound)
        assert actual <= bound, (b, t, actual, bound)
print('Lemma W check OK; max actual/bound =', round(worst,3))
# Lemma Q (odd outsiders) and E1
G = [int(p) for p in A.P1M4 if 13 <= p <= 1000]
actual = sum(1 for x in range(1, N+1, 2) if x % 25 not in (7,18) and nsf(x*x+1))
bound = sum(2*23*(N/(50*p*p)+1) for p in G) + A.E1(N, 1000)
print('Qodd actual', actual, 'bound', round(bound), 'OK' if actual <= bound else 'FAIL')

import sys, math, random
sys.argv=['x']
import analytic2 as A
from sympy import factorint, primefactors
N = 100000
def nsf(n): return any(e > 1 for e in factorint(n).values())
G = [int(p) for p in A.P1M4 if 13 <= p <= 1000]
R0 = Q1 = 1000
e1 = A.E1(N, R0)
samev, S = A.same_term(N, G)
QM = [int(q) for q in A.PR[6:A.pi(Q1)]]
sQ = sum(1.0/q**2 for q in QM)
med = sum(2*N/p**2*(sQ - (1.0/p**2 if p in QM else 0)) + 2*len(QM) for p in G)
sb = sum(2*(N/(25*p*p)+1) for p in G)
def J(b):
    T = frozenset(p for p in (3,5,7,11,13) if b % p == 0)
    v=0; m=b
    while m%2==0: m//=2; v+=1
    w = sum(1 for p in primefactors(b) if p > 29) + sum(1 for p in A.ODD29 if b%p==0)
    rho = {0:1,1:1,2:2}.get(v,4)*2**w
    tot = sb if 5 not in T else 0.0
    nother = 22 if 5 not in T else 23
    for p in G:
        QS = [q for q in (3,7,11,13) if q != p and q not in T]
        tot += nother*2*((N/(25*p*p)+1)*A.dens(QS) + 2**len(QS)-1)
    best = math.inf
    for P1 in A.GRID:
        if P1 <= Q1: continue
        best = min(best, N*A.s2(Q1,P1) + (A.pi(P1)-A.pi(Q1)) + rho*(2*(N+1)/P1 + (N*N+1)/P1**2))
    return tot + med + samev + e1 + best
rng = random.Random(3); worst = 0
for _ in range(40):
    b = 2*rng.randrange(1, N//2)
    if b % 25 in (7,18): continue
    actual = sum(1 for x in range(1, N+1) if x % 25 not in (7,18) and nsf(x*x+1) and nsf(b*x+1))
    bd = J(b); worst = max(worst, actual/bd)
    assert actual <= bd, (b, actual, bd)
print('J check OK, max actual/bound', round(worst,3))
# root-class witness: b with 169 | b^2+1, even
b = next(x for x in range(2, N, 2) if (x*x+1) % 169 == 0 and x % 25 not in (7,18))
actual = sum(1 for x in range(1, N+1) if x % 25 not in (7,18) and nsf(x*x+1) and nsf(b*x+1))
print('root-class witness b=',b,'actual',actual,'bound',round(J(b)), 'density', actual/N)

import math, random
from sympy import legendre_symbol, primerange
rng = random.Random(9); bad = 0; n = 0
for _ in range(400):
    b = rng.randrange(2, 4000)
    v = 0; m = b
    while m % 2 == 0: m //= 2; v += 1
    c = [1,1,2,4][min(v,3)]
    oddp = [p for p in primerange(3, b+1) if m % p == 0]
    sq = {}
    for y in range(b):
        if math.gcd(y, b) == 1: sq[y*y % b] = sq.get(y*y % b, 0) + 1
    for k in range(1, 200):
        if math.gcd(k, b) != 1: continue
        kinv = pow(k, -1, b); cnt = sq.get(kinv, 0); n += 1
        bound = c * math.prod(1 + legendre_symbol(k % p, p) for p in oddp)
        if cnt > bound: bad += 1
print('root-count tests', n, 'violations', bad)

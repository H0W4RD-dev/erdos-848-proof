"""Independent, memory-lean replay of ops_<N>.bin from colour4_log.py.  Shares NO code with colour4*.py.
Arithmetic: sympy.factorint (squarefreeness) and sympy.sqrt_mod (vertex set).
Checks after every vertex insertion x (in increasing order of V = {x <= N : x^2+1 not squarefree}):
  * the op group belongs to x (anchor if x=7 mod 25, outsider op if x not +-7, augmenting group ending in marker if x=18);
  * every class touched by the group is pairwise compatible: x*y+1 squarefree for all distinct members;
  * the mates form a bijection onto the 18-class vertices <= x;
  * number of classes == (x+18)//25.
Untouched classes are unchanged, hence remain proper.  So after each x the state is a proper colouring of G_x
with M(x) colours, and M is constant between consecutive vertices."""
import sys, time
import numpy as np
from sympy import factorint, primerange
from sympy.ntheory import sqrt_mod

N = int(float(sys.argv[1]))
t0 = time.time()
ops = np.memmap(f'ops_{N}.bin', dtype=np.int64, mode='r').reshape(-1, 3)
flag = bytearray(N + 1)
for p in primerange(5, N + 1):
    if p % 4 != 1: continue
    for r in sqrt_mod(-1, p * p, all_roots=True):
        for x in range(r, N + 1, p * p): flag[x] = 1
print(f"vertex set built in {time.time()-t0:.0f}s; ops={len(ops)}", flush=True)

def sf(n):
    for e in factorint(n).values():
        if e > 1: return False
    return True
def valid(c):
    m = [anchor[c]] + ([mate[c]] if mate[c] else []) + outs.get(c, [])
    for i in range(len(m)):
        for j in range(i + 1, len(m)):
            if not sf(m[i] * m[j] + 1): return False
    return True

K = (N + 18) // 25 + 1
anchor = [0] * K; mate = [0] * K; outs = {}
cnt = {}            # 18-vertex -> number of classes having it as mate
nmates = 0          # classes with a mate
n18 = 0; ncl = 0; k = 0; nops = len(ops); checks = 0; nv = 0
x = 0
for x in range(1, N + 1):
    if not flag[x]: continue
    nv += 1
    t, a, b = (int(v) for v in ops[k])
    r = x % 25
    if r == 7:
        assert t == 0 and b == x and a == ncl, ('anchor', x)
        anchor[a] = x; ncl += 1; k += 1
    elif r == 18:
        assert t == 1, ('expected mate ops', x)
        touched = []
        while int(ops[k][0]) == 1:
            _, c, e = (int(v) for v in ops[k]); k += 1
            assert c < ncl and e % 25 == 18 and e <= x
            old = mate[c]
            if old: cnt[old] -= 1
            else: nmates += 1
            mate[c] = e; cnt[e] = cnt.get(e, 0) + 1; touched.append((c, e))
        tt, xx, _ = (int(v) for v in ops[k]); k += 1
        assert tt == 3 and xx == x, ('marker', x)
        n18 += 1
        assert any(e == x for _, e in touched)
        for c, e in touched:
            assert cnt[e] == 1, ('mate twice', e)
            assert valid(c), ('mate class invalid', x, c)
        assert nmates == n18, ('mate count', x)
    else:
        assert t == 2 and b == x and a < ncl, ('outsider', x)
        outs.setdefault(a, []).append(x); k += 1
        assert valid(a), ('outsider class invalid', x, a)
    assert ncl == (x + 18) // 25, ('class count', x)
    if nv % 200000 == 0:
        print(f"  replayed through x={x} ({time.time()-t0:.0f}s)", flush=True)
assert k == nops, ('leftover ops', k, nops)
print(f"REPLAY VERIFIED N={N}: {nv} insertions; classes = M(x) after every insertion; {time.time()-t0:.0f}s", flush=True)

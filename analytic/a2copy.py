"""Erdos #848: explicit bound B(N) on |A| (A valid, A not inside the 7- or 18-class), refined cases (02_lemmas.md).
B(N) is built only from terms that are nondecreasing in N, using parameter grids that do not depend on N.
So for N0 <= N <= N1 we have |A| <= B(N1).  The interval check is B(N1) < N0/25 - 1, which gives |A| <= M(N) on [N0, N1].
Usage:
  py -3.13 analytic2.py point N...            (report B(N)/N per case)
  py -3.13 analytic2.py cover N0 Nend [eps]   (cover [N0, Nend] by intervals; adaptive eps)
"""
import math, sys, itertools
import numpy as np
from numba import njit

PZ2 = 0.45224742004106549850654336483224793417323134323989   # sum_p p^-2 (prime zeta P(2))
LIM = 2 * 10**6
_s = np.ones(LIM + 1, dtype=np.bool_); _s[:2] = False
for i in range(2, int(LIM**0.5) + 1):
    if _s[i]: _s[i*i::i] = False
PR = np.nonzero(_s)[0].astype(np.int64); del _s
def pi(x):
    """upper bound for pi(x) (exact for x <= LIM; Dusart 2010 beyond: pi(x) <= x/L (1+1/L+2.51/L^2), x >= 355991)"""
    if x <= LIM: return int(np.searchsorted(PR, int(math.floor(x)), side='right'))
    L = math.log(x); return x / L * (1 + 1 / L + 2.51 / L**2)
_CUM = np.concatenate([[0.0], np.cumsum(1.0 / PR.astype(float)**2)])
def s2(a, b):
    """upper bound for sum_{a < p <= b} p^-2 (a <= LIM)"""
    ia = int(np.searchsorted(PR, int(a), side='right'))
    if b <= LIM: return _CUM[int(np.searchsorted(PR, int(b), side='right'))] - _CUM[ia]
    return PZ2 - _CUM[ia] + 1e-15
P1M4 = PR[PR % 4 == 1]
_C1 = np.concatenate([[0.0], np.cumsum(1.0 / P1M4.astype(float)**2)])
def pi1(x): return int(np.searchsorted(P1M4, int(math.floor(x)), side='right'))
def s2m4(a, b):
    """upper bound for sum over p = 1 mod 4, a < p <= b of p^-2"""
    v = _C1[pi1(min(b, LIM))] - _C1[pi1(a)]
    if b > LIM: v += 1.0 / LIM       # sum_{n > LIM} n^-2 < 1/LIM
    return v
SMALLP = [int(p) for p in PR[:10]]   # 2..29
ODD29 = SMALLP[1:]
BIG = [int(p) for p in PR[10:60]]    # primes > 29
GRID = [float(round(1.04**k)) for k in range(90, 1300)]   # fixed parameter grid ~ 3.4e1 .. 1e22
GRID = sorted(set(g for g in GRID if g >= 30))

@njit(cache=True)
def adm_sieve(L):
    a = np.ones(L + 1, dtype=np.bool_)
    a[0] = False; a[1] = False
    for d in range(4, L + 1, 4): a[d] = False
    s = np.ones(L + 1, dtype=np.bool_); s[0] = False; s[1] = False
    for i in range(2, int(L**0.5) + 1):
        if s[i]:
            for j in range(i * i, L + 1, i): s[j] = False
    for p in range(3, L + 1, 4):
        if s[p]:
            for j in range(p, L + 1, p): a[j] = False
    r = 1
    while r * r <= L:
        a[r * r] = False; r += 1
    return a
ADML = 4 * 10**6
_ADM = adm_sieve(ADML)
ADMIDX = np.nonzero(_ADM)[0]; del _ADM
ADMLOG = np.concatenate([[0.0], np.cumsum(np.log(ADMIDX - 1.0))])
LC = math.log(2 + 3**0.5)

def pell_T(N, R1):
    """Lemma P: #{x<=N : exists u > R1 with u^2 | x^2+1} <= sum over admissible 2 <= d < D of (1 + log(N/sqrt(d-1))/log(2+sqrt3))."""
    D = (N * N + 1) / R1**2
    if D <= 2: return 0.0
    if D <= ADML:
        k = int(np.searchsorted(ADMIDX, D, side='left'))
        return k * (1 + math.log(N) / LC) - 0.5 * ADMLOG[k] / LC
    # admissible d satisfy 3 ∤ d and 4 ∤ d: at most D/2 + 2 of them below D
    return (D / 2 + 2) * (1 + math.log(N) / LC)

def E1(N, R0):
    best = math.inf
    for R1 in GRID:
        if R1 <= R0: continue
        cnt = (pi1(R1) - pi1(R0)) if R1 <= LIM else pi(R1)
        v = 2 * N * s2m4(R0, R1) + 2 * cnt + pell_T(N, R1)
        if v < best: best = v
        if pell_T(N, R1) == 0.0: break
    return best

def dens(ps):
    pr = 1.0
    for p in ps: pr *= 1 - 1.0 / p**2
    return 1 - pr

def structures(N):
    out = set()
    for v, cv in [(0, 1), (1, 1), (2, 2), (3, 4)]:
        for k in range(len(ODD29) + 1):
            for T in itertools.combinations(ODD29, k):
                base = 2**v * math.prod(T)
                if base > N: continue
                w = 0; prod = base
                while w < len(BIG) and prod * BIG[w] <= N:
                    prod *= BIG[w]; w += 1
                out.add((min(v, 1), frozenset(T), cv * 2**(k + w)))   # only (2|b, T, rho) matter
    return out

def tailW(N, rho, Lq):
    """min over P1 in GRID of [pi(P1)-pi(29)] + Lq*s2(29,P1) + rho(2(N+1)/P1 + (N^2+1)/P1^2)"""
    best = math.inf
    for P1 in GRID:
        v = (pi(P1) - 10) + Lq * s2(29, P1) + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
        if v < best: best = v
    return best

def W(N, q, excl, st, cache):
    ev, T, rho = st
    removed = set(excl) | set(T) | ({2} if ev else set()) | ({2, 5} if q == 100 else {5})
    P0 = tuple(p for p in SMALLP if p not in removed)
    key = (q, P0, rho)
    if key not in cache:
        Lq = N / q + 1
        cache[key] = Lq * dens(P0) + (2**len(P0) - 1) + tailW(N, rho, Lq)
    return cache[key]

def same_term(N, G):
    """max over b <= N of sum over p in G with p^2 | b^2+1 of (N/p^2 + 1)."""
    roots_p = {}
    for p in G:
        r = next(x for x in range(1, p) if (x * x + 1) % p == 0)
        t = (r * r + 1) // p
        s = (-t * pow(2 * r, -1, p)) % p
        R = r + s * p
        roots_p[p] = (R, p * p - R)
    best = [0.0, ()]
    def dfs(i, M, roots, val, S):
        if val > best[0]: best[0] = val; best[1] = tuple(S)
        for j in range(i, len(G)):
            p = G[j]; p2 = p * p; M2 = M * p2
            inv = pow(M, -1, p2)
            new = []
            for r0 in roots:
                for rp in roots_p[p]:
                    x = (r0 + M * (((rp - r0) * inv) % p2)) % M2
                    if M2 <= N or x <= N: new.append(x)
            if new: dfs(j + 1, M2, new, val + N / p2 + 1, S + [p])
    dfs(0, 1, [0], 0.0, [])
    return best[0], best[1]

def B(N, R0=1000, Q1=1000, verbose=False):
    G = [int(p) for p in P1M4 if 13 <= p <= R0]
    sts = structures(N)
    cache = {}
    e1 = E1(N, R0)
    qodd = sum(2 * 23 * (N / (50 * p * p) + 1) for p in G) + e1
    W25p = max(W(N, 25, (2, 5), st, cache) for st in sts)
    W100 = max(W(N, 100, (2, 5), st, cache) for st in sts)
    cls = N / 100 + 1
    samev, S = same_term(N, G)
    # J pieces
    QM = [int(q) for q in PR[6:pi(Q1)]]
    sQ = sum(1.0 / q**2 for q in QM)
    med = sum(2 * N / p**2 * (sQ - (1.0 / p**2 if p in QM else 0)) + 2 * len(QM) for p in G)
    sb = sum(2 * (N / (25 * p * p) + 1) for p in G)
    jbase = {}
    def Jb(T):
        key = frozenset(T) & frozenset((3, 5, 7, 11, 13))
        if key not in jbase:
            tot = sb if 5 not in key else 0.0
            nother = 22 if 5 not in key else 23
            for p in G:
                QS = [q for q in (3, 7, 11, 13) if q != p and q not in key]
                tot += nother * 2 * ((N / (25 * p * p) + 1) * dens(QS) + 2**len(QS) - 1)
            jbase[key] = tot
        return jbase[key]
    e2c = {}
    def E2(rho):
        if rho not in e2c:
            best = math.inf
            for P1 in GRID:
                if P1 <= Q1: continue
                v = N * s2(Q1, P1) + (pi(P1) - pi(Q1)) + rho * (2 * (N + 1) / P1 + (N * N + 1) / P1**2)
                if v < best: best = v
            e2c[rho] = best
        return e2c[rho]
    c1 = max(2 * W(N, 25, (2, 5), st, cache) + Jb(st[1]) + med + samev + e1 + E2(st[2]) for st in sts if st[0] == 1)
    c2a = max(4 * W100, cls + 2 * W100) + qodd
    c2b = max(4 * W100 + W25p, cls + 4 * W100) + qodd
    c3 = max(cls + 3 * W100 + W25p, 2 * W25p, 2 * cls)
    res = {'C1': c1, 'C2a': c2a, 'C2b': c2b, 'C3': c3}
    if verbose:
        print(f"  N={N:.3e} W25'={W25p/N:.5f} W100={W100/N:.5f} Qodd={qodd/N:.5f} E1={e1/N:.6f} same={samev/N:.5f} S={S} nstruct={len(sts)}", flush=True)
    return max(res.values()), res

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'point':
        for a in sys.argv[2:]:
            N = float(a)
            w, r = B(N, verbose=True)
            print(f"N={N:.3e}: " + ' '.join(f"{k}={v/N:.5f}" for k, v in r.items()) + f" worst={w/N:.5f} {'OK' if w < N/25 - 1 else 'FAIL'}", flush=True)
    elif mode == 'cover':
        N0 = float(sys.argv[2]); Nend = float(sys.argv[3]); eps = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
        cur = N0; k = 0
        while cur < Nend:
            N1 = min(cur * (1 + eps), Nend) if cur * (1 + eps) < Nend else cur * (1 + eps)
            w, r = B(N1)
            ok = w * (1 + 1e-9) < cur / 25 - 1
            if not ok:
                if eps < 1e-3:
                    print(f"FAIL on [{cur:.6e},{N1:.6e}]: B(N1)={w:.6e} vs N0/25-1={cur/25-1:.6e}", flush=True); sys.exit(1)
                eps /= 2; continue
            k += 1
            print(f"[{cur:.6e}, {N1:.6e}] OK  B(N1)/N0={w/cur:.5f}  eps={eps:.3f}", flush=True)
            cur = N1
            eps = min(eps * 1.5, 1.0)
        print(f"COVERED [{N0:.3e}, {cur:.3e}] with {k} intervals", flush=True)

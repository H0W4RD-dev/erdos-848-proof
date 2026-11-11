"""Independent exact re-evaluation (session c7f777) of the S1-S3 bound at chosen interval endpoints.
Written from the lemma statements in threshold/01_sharpened_lemmas.md; shares no code with threshold/*.py.
For an interval [a, b] it checks  Bound(b) < a/25 - 1  with exact Fractions (parameters chosen by a float search,
any choice being valid).  Upper-bounded constants: prime sums via ceil(10^18/p^2)/10^18 and an exact tail
1/(2(L-2)) beyond the sieve limit L; pi(x) beyond L via Dusart (x >= 355991) evaluated with a downward-rounded log."""
import sys, math, itertools
from fractions import Fraction as F
import numpy as np

L = 2 * 10**7
s = np.ones(L + 1, dtype=np.bool_); s[:2] = False
for i in range(2, int(L**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.nonzero(s)[0].astype(np.int64); del s
SC = 10**18
CI = np.concatenate([[0], np.cumsum((SC + P * P - 1) // (P * P))]).astype(object)
P14 = P[P % 4 == 1]
CI4 = np.concatenate([[0], np.cumsum((SC + P14 * P14 - 1) // (P14 * P14))]).astype(object)
TAIL = F(1, 2 * (L - 2))

def piU(x):  # integer upper bound for pi(x)
    if x <= L: return int(np.searchsorted(P, int(x), side='right'))
    lg = F(math.log(x)) * F(999999999999, 10**12)      # rounded down -> bound rounded up
    v = F(x) / lg * (1 + 1 / lg + F(251, 100) / lg**2)
    return math.floor(v) + 1
def pi14U(x):
    if x <= L: return int(np.searchsorted(P14, int(x), side='right'))
    return len(P14) + piU(x) - len(P)
def s2U(a, b):   # >= sum_{a<p<=b} p^-2 (a <= L)
    ia = int(np.searchsorted(P, int(a), side='right'))
    if b <= L: return F(int(CI[int(np.searchsorted(P, int(b), side='right'))] - CI[ia]), SC)
    return F(int(CI[-1] - CI[ia]), SC) + TAIL
def s2m4U(a, b):
    ia = int(np.searchsorted(P14, int(a), side='right'))
    if b <= L: return F(int(CI4[int(np.searchsorted(P14, int(b), side='right'))] - CI4[ia]), SC)
    return F(int(CI4[-1] - CI4[ia]), SC) + TAIL
def sqrtU(x):
    x = F(x); S = 10**9
    return F(math.isqrt(x.numerator * S * S // x.denominator) + 1, S)
def dens(ps):
    r = F(1)
    for p in ps: r *= F(p * p - 1, p * p)
    return 1 - r

ODD = [3, 5, 7, 11, 13, 17, 19, 23, 29]
SMALL = [3, 7, 11, 13, 17, 19, 23, 29]
BIGP = [int(p) for p in P[10:90]]
def structs(N, parity):
    for v in ([0] if parity == 'odd' else [1, 2, 3]):
        for k in range(10):
            for T in itertools.combinations(ODD, k):
                base = 2**v * math.prod(T)
                if base > N: continue
                w, pr = 0, base
                while pr * BIGP[w] <= N: pr *= BIGP[w]; w += 1
                yield v, T, w, base

def Gdata(J):
    Lm = math.prod(J)
    qr = {p: {x * x % p for x in range(1, p)} for p in J}
    G = [g for g in range(1, Lm + 1) if all(g % p in qr[p] for p in J)]
    S = 10**9
    return len(G), Lm, sum(F(S, math.isqrt(g * S * S)) for g in G)   # >= sum g^-1/2
_GD = {}
def tailU(N, st, P1, J):
    v, T, w, base = st
    rho = [1, 1, 2, 4][v] * 2**(len(T) + w)
    A = F(2 * (N + 1), P1) + F(N * N + 1, P1 * P1)
    if not J: return rho * A
    if J not in _GD: _GD[J] = Gdata(J)
    g, Lm, sg = _GD[J]
    return rho * (F(g, Lm) * A + g + sqrtU(F(N + 1, base)) * sg)

PG = sorted({int(1.2**k) for k in range(19, 200)})
QSG = [29, 60, 150, 400, 1000, 3000]
_SF = {}
def s2F(a, b):
    k = (a, b)
    if k not in _SF: _SF[k] = float(s2U(a, b))
    return _SF[k]
_PF = {}
def piF(x):
    if x not in _PF: _PF[x] = piU(x)
    return _PF[x]
_MIDC = {}
PGA = np.array(PG, dtype=float)
_ARR = {}
def arrays(Qs):
    if Qs not in _ARR:
        ok = np.array([(p > Qs and p <= 10**15) for p in PG])
        s2a = np.array([s2F(Qs, p) if o else np.inf for p, o in zip(PG, ok)])
        pia = np.array([piF(p) - piF(Qs) if o else np.inf for p, o in zip(PG, ok)], dtype=float)
        _ARR[Qs] = (s2a, pia)
    return _ARR[Qs]
def choose(N, st, perclass_q):   # vectorised float search for (Qs, P1, J); exact value of the chosen parameters
    v, T, w, base = st
    rho = [1, 1, 2, 4][v] * 2**(len(T) + w)
    Js = [()] + [J for k in (1, 2) for J in itertools.combinations(T, k) if math.prod(J) <= 10**5]
    A = 2 * (N + 1) / PGA + (N * N + 1) / PGA**2
    best = None
    for J in Js:
        if J:
            if J not in _GD: _GD[J] = Gdata(J)
            g, Lm, sg = _GD[J]; t = rho * (g / Lm * A + g + math.sqrt((N + 1) / base) * float(sg))
        else: t = rho * A
        for Qs in QSG:
            s2a, pia = arrays(Qs)
            val = perclass_q(Qs, exact=False) + N * s2a + pia + t
            i = int(np.argmin(val))
            if best is None or val[i] < best[0]: best = (val[i], Qs, PG[i], J)
    _, Qs, P1, J = best
    return perclass_q(Qs, exact=True) + N * s2U(Qs, P1) + piU(P1) - piU(Qs) + tailU(N, st, P1, J)

def Wgroup(N, groups, parity):
    """S1 bound: max over witness structures; groups = list of (count, modulus)."""
    worst = F(0); seen = set()
    for st in structs(N, parity):
        v, T, w, base = st
        key = (T, [1, 1, 2, 4][v] * 2**(len(T) + w), base)
        if key in seen: continue
        seen.add(key)
        P0 = [p for p in SMALL if p not in T]
        d = dens(P0)
        def pc(Qs, exact):
            if exact:
                return sum(c * ((F(N, q) + 1) * d + 2**len(P0) - 1 + F(N, q) * s2U(29, Qs) + piU(Qs) - 10) for c, q in groups)
            return sum(c * ((N / q + 1) * float(d) + 2**len(P0) - 1 + N / q * s2F(29, Qs) + piF(Qs) - 10) for c, q in groups)
        val = choose(N, st, pc)
        worst = max(worst, val)
    return worst

def E1U(N, R0=1000):
    best = None
    c = F(1762747174039, 10**12)   # <= log(3+2sqrt2)
    lgN = F(math.log(N)) * F(1000000000001, 10**12)
    for R1 in PG:
        if R1 <= R0 or R1 > 10**15: continue
        D = F(N * N + 1, R1 * R1)
        pell = 0 if D <= 2 else (D / 2 + 1) * (1 + lgN / c)
        v = 2 * N * s2m4U(R0, R1) + 2 * (pi14U(R1) - pi14U(R0)) + pell
        if best is None or v < best: best = v
    return best

G = [int(p) for p in P14 if 13 <= p <= 1000]
def Qclasses(N, ncls, Q, e1): return sum(2 * ncls * (F(N, Q * p * p) + 1) for p in G) + e1

def sameU(N, nc):
    roots = {}
    for p in G:
        r = next(x for x in range(1, p) if (x * x + 1) % p == 0)
        t = (r * r + 1) // p; R = r + ((-t * pow(2 * r, -1, p)) % p) * p
        roots[p] = (R, p * p - R)
    wt = {p: min(nc * (F(N, 25 * p * p) + 1), F(N, p * p) + 1) for p in G}
    if N > 10**12: return sum(wt.values())
    best = [F(0)]
    def dfs(i, M, rs, val):
        best[0] = max(best[0], val)
        for j in range(i, len(G)):
            p = G[j]; M2 = M * p * p; inv = pow(M, -1, p * p)
            new = [x for r0 in rs for rp in roots[p] for x in [(r0 + M * (((rp - r0) * inv) % (p * p))) % M2] if M2 <= N or x <= N]
            if new: dfs(j + 1, M2, new, val + wt[p])
    dfs(0, 1, [0], F(0))
    return best[0]

def C1U(N, e1):
    same = {nc: sameU(N, nc) for nc in (22, 23)}
    worst = F(0); seen = set()
    for st in structs(N, 'even'):
        v, T, w, base = st
        key = (T, [1, 1, 2, 4][v] * 2**(len(T) + w), base)
        if key in seen: continue
        seen.add(key)
        P0 = [p for p in SMALL if p not in T]
        nother = 23 if 5 in T else 22
        fixed = 2 * ((F(N, 25) + 1) * dens(P0) + 2**len(P0) - 1)
        if 5 not in T: fixed += sum(2 * (F(N, 25 * p * p) + 1) for p in G)
        for p in G:
            QS = [q for q in (3, 7, 11, 13) if q != p and q not in T]
            fixed += nother * 2 * ((F(N, 25 * p * p) + 1) * dens(QS) + 2**len(QS) - 1)
        fixed += same[nother] + e1
        ffixed = float(fixed)
        def pc(Qs, exact):
            QM = [int(q) for q in P[6:piU(Qs)] if q not in T]
            if exact:
                sQ = sum(F(1, q * q) for q in QM)
                med = sum(2 * F(N, p * p) * sQ + 2 * len(QM) for p in G)
                return fixed + 2 * (F(N, 25) * s2U(29, Qs) + piU(Qs) - 10) + med
            sQ = sum(1.0 / q / q for q in QM)
            return ffixed + 2 * (N / 25 * s2F(29, Qs) + piF(Qs) - 10) + sum(2 * N / p / p * sQ + 2 * len(QM) for p in G)
        worst = max(worst, choose(N, st, pc))
    return worst

def bound(N):
    e1 = E1U(N)
    cls = F(N, 100) + 1
    W = lambda c, q, par: Wgroup(N, [(c, q)], par)
    W6o, W4o, W3o, W2o, W1o = W(6, 100, 'odd'), W(4, 100, 'odd'), W(3, 100, 'odd'), W(2, 100, 'odd'), W(1, 100, 'odd')
    W1e, W2e, W25e = W(1, 100, 'even'), W(2, 100, 'even'), W(1, 25, 'even')
    Qodd, Qu = Qclasses(N, 23, 50, e1), Qclasses(N, 23, 100, e1)
    c2i = Qodd + W6o + W2o
    c2ii = Qu + max(W6o + 2 * W1o, W6o + W1e, W4o + cls)
    c3 = max(2 * W25e, W2e + W3o + max(W1o, cls), 2 * max(2 * W1o, cls))
    c1 = C1U(N, e1)
    return {'C1': c1, 'C2i': c2i, 'C2ii': c2ii, 'C3': c3}

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        a, b = (int(float(t)) for t in arg.split(','))
        r = bound(b)
        worst = max(r.values())
        ok = worst < F(a, 25) - 1
        print(f"[{a}, {b}]: " + ' '.join(f"{k}={float(v)/a:.5f}" for k, v in r.items()) + f"  worst/a={float(worst)/a:.5f}  {'OK' if ok else 'FAIL'}", flush=True)

"""verify_exactL: verify_exact.py for bound3L (large N). Same code except: imports bound3L; prime sums beyond PLIM get the exact
rational tail 1/(2(PLIM-2)); pi beyond PLIM is Dusart (log lowered by 1e-12, +1); same-term is the full G sum for N > 1e12.
Exact-rational re-verification of a bound3.py cover log.
For each interval [a_i, a_{i+1}] (endpoints = integers rounded from the log, chained), checks
    B3exact(a_{i+1}) < a_i/25 - 1     in exact rational arithmetic.
Parameters (Qs, P1, J, R1) are chosen by float search (any choice gives a valid bound); the VALUE is then
evaluated exactly with rational upper bounds for every transcendental quantity:
  * sum_{a<p<=b} p^-2  <=  sum ceil(10^18/p^2)/10^18      (exact integer sums)
  * sqrt(x)  <= (isqrt(floor(x*S^2))+1)/S,   g^-1/2 <= S/isqrt(g*S^2),  S = 10^9
  * log N <= float log * (1+1e-12);  log(3+2sqrt2) >= 1.762747174039
Usage: py -3.13 verify_exact.py cover3_xxx.log
"""
import sys, math, itertools, re
from fractions import Fraction as F
argv = sys.argv[1:]; sys.argv = ['x']
import numpy as np
import bound3L as B

S = 10**9
PRL = B.PR
SC = 10**18
_ceil = (SC + PRL * PRL - 1) // (PRL * PRL)            # int64 exact (p^2 < 3.6e15)
CUMI = np.concatenate([[0], np.cumsum(_ceil)]).astype(object)
P1M4 = B.P1M4
_ceil4 = (SC + P1M4 * P1M4 - 1) // (P1M4 * P1M4)
CUMI4 = np.concatenate([[0], np.cumsum(_ceil4)]).astype(object)
TAILQ = F(1, 2 * (B.PLIM - 2))
def s2x(a, b):
    if b <= B.PLIM: return F(int(CUMI[B.pi(b)] - CUMI[B.pi(a)]), SC)
    return F(int(CUMI[-1] - CUMI[B.pi(a)]), SC) + TAILQ
def s2m4x(a, b):
    if b <= B.PLIM: return F(int(CUMI4[B.pi1(b)] - CUMI4[B.pi1(a)]), SC)
    return F(int(CUMI4[-1] - CUMI4[B.pi1(a)]), SC) + TAILQ
def sqrt_up(x):
    x = F(x); return F(math.isqrt(x.numerator * S * S // x.denominator) + 1, S)
LOG_LC = F(1762747174039, 10**12)
def log_up(N): return F(math.log(N) * (1 + 1e-12))
def densx(ps):
    pr = F(1)
    for p in ps: pr *= 1 - F(1, p * p)
    return 1 - pr
_GJx = {}
def GJx(J):
    if J not in _GJx:
        L = math.prod(J)
        qr = {p: set((x * x) % p for x in range(1, p)) for p in J}
        G = [g for g in range(1, L + 1) if all((g % p) in qr[p] for p in J)]
        sg = sum(F(S, math.isqrt(g * S * S)) for g in G)
        _GJx[J] = (len(G), L, sg)
    return _GJx[J]

def tail_x(N, st, P1, J):
    v, T, w, base = st
    rho = B.CV[v] * 2**(len(T) + w)
    A = F(2 * (N + 1), P1) + F(N * N + 1, P1 * P1)
    if not J: return rho * A
    g, L, sg = GJx(J)
    return rho * (F(g, L) * A + g + sqrt_up(F(N + 1, base)) * sg)

def best_params(N, st, lo_of_Qs):
    """float search of (Qs, P1, J) minimising lo_of_Qs(Qs) + shared(Qs,P1) + tail; returns choice."""
    v, T, w, base = st
    rho = B.CV[v] * 2**(len(T) + w)
    Jl = [()] + [J for k in range(1, len(T) + 1) for J in itertools.combinations(T, k) if math.prod(J) <= 3 * 10**5]
    A = 2 * (N + 1) / B.PGA + (N * N + 1) / B.PGA**2
    tails = []
    for J in Jl:
        if not J: tails.append(rho * A)
        else:
            g, L, sg = B.GJ(J); tails.append(rho * (g / L * A + g + math.sqrt((N + 1) / base) * sg))
    tails = np.array(tails); jidx = tails.argmin(axis=0); tmin = tails.min(axis=0)
    best = (math.inf, None)
    for Qs in B.QSGRID:
        m = B.PGA > Qs
        val = lo_of_Qs(Qs) + N * (B.S2P - B.s2(0, Qs)) + (B.PIP - B.pi(Qs)) + tmin
        val = np.where(m, val, np.inf)
        i = int(val.argmin())
        if val[i] < best[0]: best = (val[i], (Qs, int(B.PGRID[i]), Jl[int(jidx[i])]))
    return best[1]

def shared_x(N, Qs, P1): return N * s2x(Qs, P1) + (B.pi(P1) - B.pi(Qs))

def WGx(N, groups, parity):
    worst = F(0); seen = set()
    for st in B.structures(N, parity):
        v, T, w, base = st
        rho = B.CV[v] * 2**(len(T) + w)
        key = (T, rho, base)
        if key in seen: continue
        seen.add(key)
        P0 = tuple(p for p in B.SIEVE0 if p not in T)
        mid_f = lambda Qs: sum(c * (N / q * B.s2(29, Qs) + (B.pi(Qs) - 10)) for c, q in groups)
        Qs, P1, J = best_params(N, st, mid_f)
        d = densx(P0)
        val = sum(c * ((F(N, q) + 1) * d + 2**len(P0) - 1) for c, q in groups)
        val += sum(c * (F(N, q) * s2x(29, Qs) + (B.pi(Qs) - 10)) for c, q in groups)
        val += shared_x(N, Qs, P1) + tail_x(N, st, P1, J)
        if val > worst: worst = val
    return worst

def E1x(N, R0):
    best = (math.inf, None)
    for R1 in B.PGRID:
        if R1 <= R0: continue
        fv = 2 * N * B.s2m4(R0, R1) + 2 * (B.pi1(R1) - B.pi1(R0)) + B.pell_T(N, R1)
        if fv < best[0]: best = (fv, R1)
    R1 = best[1]
    D = F(N * N + 1, R1 * R1)
    pell = 0 if D <= 2 else (D / 2 + 1) * (1 + log_up(N) / LOG_LC)
    return 2 * N * s2m4x(R0, R1) + 2 * (B.pi1(R1) - B.pi1(R0)) + pell

def Qx(N, ncls, Q, R0, e1):
    G = [p for p in B.G13 if p <= R0]
    return sum(2 * ncls * (F(N, Q * p * p) + 1) for p in G) + e1

def same_x(N, G, nc):
    roots_p = {}
    for p in G:
        r = next(x for x in range(1, p) if (x * x + 1) % p == 0)
        t = (r * r + 1) // p; s = (-t * pow(2 * r, -1, p)) % p; R = r + s * p
        roots_p[p] = (R, p * p - R)
    wt = {p: min(nc * (F(N, 25 * p * p) + 1), F(N, p * p) + 1) for p in G}
    if N > 10**12:
        return sum(wt[p] for p in G)
    best = [F(0)]
    def dfs(i, M, roots, val):
        if val > best[0]: best[0] = val
        for j in range(i, len(G)):
            p = G[j]; p2 = p * p; M2 = M * p2; inv = pow(M, -1, p2); new = []
            for r0 in roots:
                for rp in roots_p[p]:
                    x = (r0 + M * (((rp - r0) * inv) % p2)) % M2
                    if M2 <= N or x <= N: new.append(x)
            if new: dfs(j + 1, M2, new, val + wt[p])
    dfs(0, 1, [0], F(0))
    return best[0]

_MED = {}
def med_x(N, G, Qs, Tbig):
    key = (N, Qs, Tbig)
    if key not in _MED:
        QM = [int(q) for q in B.PR[6:B.pi(Qs)] if q not in Tbig]
        SQ = sum(F(1, q * q) for q in QM); QMs = set(QM)
        tot = F(0)
        for p in G:
            inm = p in QMs
            tot += 2 * N * F(1, p * p) * (SQ - (F(1, p * p) if inm else 0)) + 2 * (len(QM) - inm)
        _MED[key] = tot
    return _MED[key]

def C1x(N, R0, e1):
    G = [p for p in B.G13 if p <= R0]
    samev = {nc: same_x(N, G, nc) for nc in (22, 23)}
    worst = F(0); seen = set()
    for st in B.structures(N, 1):
        v, T, w, base = st
        rho = B.CV[v] * 2**(len(T) + w)
        key = (T, rho, base)
        if key in seen: continue
        seen.add(key)
        Tq = set(T)
        P0 = tuple(p for p in B.SIEVE0 if p not in Tq)
        def lo_f(Qs):
            QM = [int(q) for q in B.PR[6:B.pi(Qs)] if q not in Tq]
            sQ = sum(1.0 / q**2 for q in QM)
            med = sum(2 * N / p**2 * (sQ - (1.0 / p**2 if p in QM else 0)) + 2 * len(QM) for p in G)
            return 2 * (N / 25 * B.s2(29, Qs) + (B.pi(Qs) - 10)) + med
        Qs, P1, J = best_params(N, st, lo_f)
        val = 2 * ((F(N, 25) + 1) * densx(P0) + 2**len(P0) - 1)
        val += 2 * (F(N, 25) * s2x(29, Qs) + (B.pi(Qs) - 10))
        val += med_x(N, G, Qs, frozenset(Tq & {17, 19, 23, 29}))
        if 5 not in Tq: val += sum(2 * (F(N, 25 * p * p) + 1) for p in G)
        nother = 22 if 5 not in Tq else 23
        for p in G:
            QS = [q for q in (3, 7, 11, 13) if q != p and q not in Tq]
            val += nother * 2 * ((F(N, 25 * p * p) + 1) * densx(QS) + 2**len(QS) - 1)
        val += samev[nother] + e1 + shared_x(N, Qs, P1) + tail_x(N, st, P1, J)
        if val > worst: worst = val
    return worst

def B3x(N):
    R0 = 1000
    cls = F(N, 100) + 1
    W = {}
    for name, gr, par in [('6o', [(6, 100)], 0), ('4o', [(4, 100)], 0), ('2o', [(2, 100)], 0), ('1o', [(1, 100)], 0),
                          ('3o', [(3, 100)], 0), ('1e', [(1, 100)], 1), ('2e', [(2, 100)], 1), ('25e', [(1, 25)], 1)]:
        W[name] = WGx(N, gr, par)
    e1 = E1x(N, R0)
    Qodd = Qx(N, 23, 50, R0, e1); Qu = Qx(N, 23, 100, R0, e1)
    c2i = Qodd + W['6o'] + W['2o']
    c2ii = Qu + max(W['6o'] + 2 * W['1o'], W['6o'] + W['1e'], W['4o'] + cls)
    c3 = max(2 * W['25e'], W['2e'] + W['3o'] + max(W['1o'], cls), 2 * max(2 * W['1o'], cls))
    c1 = C1x(N, R0, e1)
    return max(c1, c2i, c2ii, c3), {'C1': c1, 'C2i': c2i, 'C2ii': c2ii, 'C3': c3}

if __name__ == '__main__':
    pts = []
    for line in open(argv[0]):
        m = re.match(r'\[([\d.e+]+), ([\d.e+]+)\] OK', line)
        if m:
            a, b = float(m.group(1)), float(m.group(2))
            if not pts: pts.append(round(a))
            pts.append(round(b))
    print(f'{len(pts)-1} intervals from {pts[0]} to {pts[-1]}', flush=True)
    bad = 0
    done = set()
    if len(argv) > 1:   # resume: skip intervals already reported OK in a previous output log
        for line in open(argv[1], encoding='utf8'):
            m = re.match(r'\[(\d+), (\d+)\] exact .* OK$', line.strip())
            if m: done.add((int(m.group(1)), int(m.group(2))))
    for a, b in zip(pts, pts[1:]):
        if (a, b) in done:
            print(f'[{a}, {b}] previously verified (resume)', flush=True); continue
        w, r = B3x(b)
        ok = w < F(a, 25) - 1
        bad += not ok
        print(f'[{a}, {b}] exact B3/a = {float(w / a):.6f} {"OK" if ok else "FAIL"}', flush=True)
    print('ALL INTERVALS VERIFIED EXACTLY' if bad == 0 else f'{bad} FAILURES', flush=True)

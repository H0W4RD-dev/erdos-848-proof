"""bound3L: bound3 extended to large N (N up to ~1e20). Differences from bound3.py ONLY:
  pi, pi1 beyond PLIM via Dusart 2010 upper bound; prime sums beyond PLIM bounded by sum over odd n > PLIM <= 1/(2(PLIM-2));
  PGRID extended to 1e16; same-term replaced by the full sum over G (valid upper bound for every b) when N > 1e12.
Erdos #848, threshold side: sharpened explicit bound B3(N) (float exploration version).
Independent of ../analytic2.py (re-implemented).  Lemmas: threshold/01_sharpened_lemmas.md.
New ingredients vs 02_lemmas.md:
  S1 witness-group sharing: all a <= N with q^2 | ab+1 for one witness b, q > Qs, are counted ONCE
     (one class mod q^2 per q, plus one rho-tail), whatever classes they lie in.
  S2 character-refined tail: sum_{k<=K} rho_k(b) <= c_v 2^omega * #{k<=K : k nonzero QR mod all p in J}, J subset T.
  S3 new case split for C2/C3 by residues mod 4 (fewer witnesses; Q_u when A* is in one residue mod 4).
  S4 'same' term restricted to the 22 classes mod 25 (class s_b already counted).
  S5 Pell ratio 3+2*sqrt2 (d=3 is inadmissible).
Usage: py -3.13 bound3.py point N ...   |   cover N0 N1 [eps]
"""
import math, sys, itertools
import numpy as np

PLIM = 6 * 10**7
_s = np.ones(PLIM + 1, dtype=np.bool_); _s[:2] = False
for i in range(2, int(PLIM**0.5) + 1):
    if _s[i]: _s[i*i::i] = False
PR = np.nonzero(_s)[0].astype(np.int64); del _s
_CUM = np.concatenate([[0.0], np.cumsum(1.0 / PR.astype(float)**2)])
P1M4 = PR[PR % 4 == 1]
_C1 = np.concatenate([[0.0], np.cumsum(1.0 / P1M4.astype(float)**2)])
SLACK = 1e-9   # absolute slack on float prime sums (rounding error <= 3.6e6 * 2^-53 * 0.46 < 2e-10)

TAILX = 1.0 / (2 * (PLIM - 2))   # sum over primes p > PLIM of p^-2 <= sum over odd n > PLIM of n^-2 <= 1/(2(PLIM-2))
def dusart(x):
    L = math.log(x) * (1 - 1e-12)
    return math.ceil(x / L * (1 + 1 / L + 2.51 / L**2)) + 1   # +1 absorbs float rounding
PI_PLIM = None
def pi(x):
    if x <= PLIM: return int(np.searchsorted(PR, int(math.floor(x)), side='right'))
    return dusart(x)
def pi1(x):
    if x <= PLIM: return int(np.searchsorted(P1M4, int(math.floor(x)), side='right'))
    return len(P1M4) + dusart(x) - len(PR)
def s2(a, b):
    if b <= PLIM: return _CUM[pi(b)] - _CUM[pi(a)] + SLACK
    return _CUM[-1] - _CUM[pi(a)] + TAILX + SLACK
def s2m4(a, b):
    if b <= PLIM: return _C1[pi1(b)] - _C1[pi1(a)] + SLACK
    return _C1[-1] - _C1[pi1(a)] + TAILX + SLACK

ODD29 = [3, 5, 7, 11, 13, 17, 19, 23, 29]
SIEVE0 = [3, 7, 11, 13, 17, 19, 23, 29]       # primes <= 29 sieved exactly (2,5 always excluded)
BIG = [int(p) for p in PR[10:80]]
LC = math.log(3 + 2 * math.sqrt(2))          # S5
G13 = [int(p) for p in P1M4 if 13 <= p <= 5000]

def dens(ps):
    pr = 1.0
    for p in ps: pr *= 1 - 1.0 / p**2
    return 1 - pr

def structures(N, parity):
    """witness structures (v, T, w): b = 2^v' * prod(T) * (w primes > 29), v=3 means v'>=3. parity 0 odd,1 even, None any."""
    out = []
    vs = [0] if parity == 0 else ([1, 2, 3] if parity == 1 else [0, 1, 2, 3])
    for v in vs:
        for k in range(len(ODD29) + 1):
            for T in itertools.combinations(ODD29, k):
                base = 2**v * math.prod(T)
                if base > N: continue
                w = 0; prod = base
                while prod * BIG[w] <= N:
                    prod *= BIG[w]; w += 1
                out.append((v, T, w, base))
    return out

CV = {0: 1, 1: 1, 2: 2, 3: 4}
_GJ = {}
def GJ(J):
    """(|G|, L, sum_{g in G} g^-1/2) for G = residues mod L=prod J that are nonzero QRs mod each p in J."""
    if J not in _GJ:
        L = math.prod(J)
        qr = {p: set((x * x) % p for x in range(1, p)) for p in J}
        G = [g for g in range(1, L + 1) if all((g % p) in qr[p] for p in J)]
        _GJ[J] = (len(G), L, sum(g**-0.5 for g in G))
    return _GJ[J]

_JL = {}
def Jlist(T):
    """candidate J subsets of T (product <= 3e5): list of (g/L, g, sg); J empty -> (1, 0, 0) handled as Lemma W."""
    if T not in _JL:
        out = [(1.0, 0.0, 0.0)]
        for k in range(1, len(T) + 1):
            for J in itertools.combinations(T, k):
                if math.prod(J) > 3 * 10**5: continue
                g, L, sg = GJ(J); out.append((g / L, float(g), sg))
        _JL[T] = np.array(out)
    return _JL[T]

def tailvec(N, st, P1v):
    """vector over P1v of the S2 tail bound."""
    v, T, w, base = st
    rho = CV[v] * 2**(len(T) + w)
    A = 2 * (N + 1) / P1v + (N * N + 1) / P1v**2
    J = Jlist(T)
    vals = rho * (J[:, 0:1] * A[None, :] + J[:, 1:2] + math.sqrt((N + 1) / base) * J[:, 2:3])
    return vals.min(axis=0)

PGRID = sorted(set(int(round(1.05**k)) for k in range(60, 760)))   # ~ 1.9e1 .. 1.4e16 (fixed)
PGRID = [p for p in PGRID if 30 <= p]
PGA = np.array(PGRID, dtype=float)
S2P = np.array([s2(0, p) for p in PGRID])
PIP = np.array([pi(p) for p in PGRID], dtype=float)
QSGRID = [29, 50, 100, 200, 400, 1000, 2000, 5000]

def WG(N, groups, parity, extra_small=0.0):
    """Max over witness structures of the bound for #{a in union of classes : mu(ab+1)=0}.
    groups: list of (c, q): c classes mod q (q in 25,50,100), on each 4 and 25 never divide ab+1.
    Sieve: p<=29 incl-excl per class; 29<p<=Qs per class; Qs<p<=P1 shared (N/p^2+1); p>P1 tail."""
    worst = 0.0
    cache = {}
    for st in structures(N, parity):
        v, T, w, base = st
        P0 = tuple(p for p in SIEVE0 if p not in T)
        rho = CV[v] * 2**(len(T) + w)
        key = (P0, rho, base)
        if key in cache: continue
        main = sum(c * ((N / q + 1) * dens(P0) + 2**len(P0) - 1) for c, q in groups)
        tv = tailvec(N, st, PGA)
        best = math.inf
        for Qs in QSGRID:
            mid = sum(c * (N / q * s2(29, Qs) + (pi(Qs) - 10)) for c, q in groups)
            m = PGA > Qs
            val = mid + N * (S2P[m] - s2(0, Qs)) + (PIP[m] - pi(Qs)) + tv[m]
            best = min(best, float(val.min()))
        cache[key] = main + best
        worst = max(worst, cache[key])
    return worst

# ---------- x^2+1 side ----------
def pell_T(N, R1):
    D = (N * N + 1) / R1**2
    if D <= 2: return 0.0
    # admissible d < D: 3 does not divide d, 4 does not divide d -> at most D/2+1
    return (D / 2 + 1) * (1 + math.log(N) / LC)

def E1(N, R0):
    best = math.inf
    for R1 in PGRID:
        if R1 <= R0: continue
        v = 2 * N * s2m4(R0, R1) + 2 * (pi1(R1) - pi1(R0)) + pell_T(N, R1)
        best = min(best, v)
    return best

def Qcount(N, ncls, Q, R0):
    """#{x<=N : x in ncls classes mod Q avoiding +-7 mod 25, mu(x^2+1)=0} (Lemma Q)."""
    G = [p for p in G13 if p <= R0]
    return sum(2 * ncls * (N / (Q * p * p) + 1) for p in G) + E1(N, R0)

def same_term(N, G, weight):
    roots_p = {}
    for p in G:
        r = next(x for x in range(1, p) if (x * x + 1) % p == 0)
        t = (r * r + 1) // p; s = (-t * pow(2 * r, -1, p)) % p; R = r + s * p
        roots_p[p] = (R, p * p - R)
    if N > 1e12:
        return sum(weight(p) for p in G)
    best = [0.0]
    def dfs(i, M, roots, val):
        if val > best[0]: best[0] = val
        for j in range(i, len(G)):
            p = G[j]; p2 = p * p; M2 = M * p2; inv = pow(M, -1, p2); new = []
            for r0 in roots:
                for rp in roots_p[p]:
                    x = (r0 + M * (((rp - r0) * inv) % p2)) % M2
                    if M2 <= N or x <= N: new.append(x)
            if new: dfs(j + 1, M2, new, val + weight(p))
    dfs(0, 1, [0], 0.0)
    return best[0]

def C1(N, R0=1000):
    """b in A* even. |A| <= [A7,A18 classes, p<=Qs] + J_small(b) + shared(Qs<q<=P1) + tail."""
    G = [p for p in G13 if p <= R0]
    e1 = E1(N, R0)
    samev = {nc: same_term(N, G, lambda p, nc=nc: min(nc * (N / (25 * p * p) + 1), N / (p * p) + 1)) for nc in (22, 23)}
    worst = 0.0
    seen = {}
    for st in structures(N, 1):
        v, T, w, base = st
        rho = CV[v] * 2**(len(T) + w)
        key = (T, rho, base)
        if key in seen: continue
        P0 = tuple(p for p in SIEVE0 if p not in T)
        main7_18 = 2 * ((N / 25 + 1) * dens(P0) + 2**len(P0) - 1)
        Tq = set(T)
        sb = sum(2 * (N / (25 * p * p) + 1) for p in G) if 5 not in Tq else 0.0
        nother = 22 if 5 not in Tq else 23
        oth = 0.0
        for p in G:
            QS = [q for q in (3, 7, 11, 13) if q != p and q not in Tq]
            oth += nother * 2 * ((N / (25 * p * p) + 1) * dens(QS) + 2**len(QS) - 1)
        tv = tailvec(N, st, PGA)
        best = math.inf
        for Qs in QSGRID:
            mid78 = 2 * (N / 25 * s2(29, Qs) + (pi(Qs) - 10))
            QM = [int(q) for q in PR[6:pi(Qs)] if q not in Tq]
            sQ = sum(1.0 / q**2 for q in QM)
            med = sum(2 * N / p**2 * (sQ - (1.0 / p**2 if p in QM else 0)) + 2 * len(QM) for p in G)
            m = PGA > Qs
            val = mid78 + med + N * (S2P[m] - s2(0, Qs)) + (PIP[m] - pi(Qs)) + tv[m]
            best = min(best, float(val.min()))
        seen[key] = main7_18 + sb + oth + samev[nother] + e1 + best
        worst = max(worst, seen[key])
    return worst

def B3(N, verbose=False):
    R0 = 1000
    cls = N / 100 + 1
    Wodd = lambda gr: WG(N, gr, 0)
    Wev = lambda gr: WG(N, gr, 1)
    Wany = lambda gr: WG(N, gr, None)
    W6o, W4o, W2o, W1o, W3o = Wodd([(6, 100)]), Wodd([(4, 100)]), Wodd([(2, 100)]), Wodd([(1, 100)]), Wodd([(3, 100)])
    W1e, W2e, W25e = Wev([(1, 100)]), Wev([(2, 100)]), Wev([(1, 25)])
    Qodd = Qcount(N, 23, 50, R0)
    Qu = Qcount(N, 23, 100, R0)
    c2i = Qodd + W6o + W2o
    c2ii = Qu + max(W6o + 2 * W1o, W6o + W1e, W4o + cls)
    c3 = max(2 * W25e, W2e + W3o + max(W1o, cls), 2 * max(2 * W1o, cls))
    c1 = C1(N, R0)
    res = {'C1': c1, 'C2i': c2i, 'C2ii': c2ii, 'C3': c3}
    if verbose:
        print(f"  N={N:.3e} W6o={W6o/N:.5f} W1o={W1o/N:.5f} W1e={W1e/N:.5f} W25e={W25e/N:.5f} Qodd={Qodd/N:.5f} Qu={Qu/N:.5f}", flush=True)
    return max(res.values()), res

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'point':
        for a in sys.argv[2:]:
            N = float(a); w, r = B3(N, verbose=True)
            print(f"N={N:.3e}: " + ' '.join(f"{k}={v/N:.5f}" for k, v in r.items()) + f" worst={w/N:.5f} {'OK' if w < N/25-1 else 'FAIL'}", flush=True)
    elif mode == 'cover':
        N0 = float(sys.argv[2]); Nend = float(sys.argv[3]); eps = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
        cur = N0; k = 0
        while cur < Nend:
            N1 = cur * (1 + eps)
            w, r = B3(N1)
            if not w * (1 + 1e-9) < cur / 25 - 1:
                if eps < 1e-3: print(f"FAIL at {cur:.6e}", flush=True); sys.exit(1)
                eps /= 2; continue
            k += 1; print(f"[{cur:.6e}, {N1:.6e}] OK B(N1)/N0={w/cur:.5f} eps={eps:.3f}", flush=True)
            cur = N1; eps = min(eps * 1.5, 0.5)
        print(f"COVERED [{N0:.3e}, {cur:.3e}] with {k} intervals", flush=True)

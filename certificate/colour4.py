"""Streaming colouring certificate for Erdos #848 (memory O(window)).
Usage: py -3.13 colour4.py NMAX [seed] [logfile]
Same certificate logic as colour3 (see 01_colouring_certificate.md), but:
 * vertices streamed in segments (root sieve for p <= sqrt(N), precomputed hit list for p > sqrt(N));
 * only the most recent W classes are active (older classes are frozen: never modified again);
 * 18-class vertices placed by bounded BFS augmenting paths; outsiders by random sampling.
Every placement is validated with the exact squarefree test, so the state after inserting all
vertices <= x is a proper colouring of G_x with exactly M(x) classes (Lemma C2).
Checkpoint lines are printed every 1e7."""
import sys, time, math
import numpy as np
from numba import njit, uint64, int64

NMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**7
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
assert NMAX <= 3 * 10**9

def primes_upto(L):
    s = np.ones(L + 1, dtype=np.bool_); s[:2] = False
    for i in range(2, int(L**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

CB = int(round((NMAX * NMAX + 1) ** (1/3))) + 10
_P = primes_upto(CB); _P = _P[_P > 2]
PRS = _P.astype(np.uint64)
INV = np.array([pow(int(p), -1, 1 << 64) for p in _P], dtype=np.uint64)
LIM = np.array([((1 << 64) - 1) // int(p) for p in _P], dtype=np.uint64)

@njit(cache=True)
def isqrt64(n):
    r = uint64(math.sqrt(float(n)))
    while r * r > n: r -= uint64(1)
    while (r + uint64(1)) * (r + uint64(1)) <= n: r += uint64(1)
    return r

@njit(cache=True)
def squarefree(n, PRS, INV, LIM):
    m = uint64(n)
    tz = 0
    while m & uint64(1) == uint64(0):
        m >>= uint64(1); tz += 1
        if tz >= 2: return False
    for i in range(PRS.shape[0]):
        p = PRS[i]
        if p * p * p > m: break
        if m * INV[i] <= LIM[i]:
            m = m // p
            if m * INV[i] <= LIM[i]:
                return False
    if m > uint64(1):
        r = isqrt64(m)
        if r * r == m: return False
    return True

@njit(cache=True)
def powmod(a, e, m):
    r = int64(1); a = int64(a % m)
    while e > 0:
        if e & 1: r = (r * a) % m
        a = (a * a) % m; e >>= 1
    return r

@njit(cache=True)
def root_lift(p):
    c = 2
    while True:
        r = powmod(c, (p - 1) // 4, p)
        if (r * r) % p == p - 1: break
        c += 1
    t = (r * r + 1) // p
    s = ((p - (t % p)) * powmod((2 * r) % p, p - 2, p)) % p
    return r + s * p   # R with R^2 = -1 mod p^2

@njit(cache=True)
def small_roots(smallp):
    out = np.zeros((smallp.shape[0], 3), dtype=np.int64)
    k = 0
    for i in range(smallp.shape[0]):
        p = smallp[i]
        if p % 4 != 1: continue
        R = root_lift(p)
        out[k, 0] = p * p; out[k, 1] = R; out[k, 2] = p * p - R; k += 1
    return out[:k]

@njit(cache=True)
def large_hits(N, sq, small):
    """x <= N with p^2 | x^2+1 for a prime p = 1 mod 4, p > sq (at most one x per root)"""
    hits = np.zeros(1 << 20, dtype=np.int64); nh = 0
    SEG = 1 << 22
    lo = sq + 1
    while lo <= N:
        hi = min(N + 1, lo + SEG)
        n = hi - lo
        seg = np.ones(n, dtype=np.bool_)
        for i in range(small.shape[0]):
            p = int64(small[i])
            if p * p >= hi: break
            st = max(p * p, ((lo + p - 1) // p) * p)
            for j in range(st - lo, n, p): seg[j] = False
        for j in range(n):
            if not seg[j]: continue
            p = int64(lo + j)
            if p % 4 != 1: continue
            R = root_lift(p)
            p2 = p * p
            if R <= N:
                hits[nh] = R; nh += 1
            if p2 - R <= N:
                hits[nh] = p2 - R; nh += 1
        lo = hi
    return np.sort(hits[:nh])

@njit(cache=True)
def compat_free(e, c, W, mate, outs, nout, PRS, INV, LIM, skip_mate):
    s = c % W
    a = 7 + 25 * c
    if not squarefree(uint64(e) * uint64(a) + uint64(1), PRS, INV, LIM): return False
    if (not skip_mate) and mate[s] >= 0:
        if not squarefree(uint64(e) * uint64(mate[s]) + uint64(1), PRS, INV, LIM): return False
    for k in range(nout[s]):
        if not squarefree(uint64(e) * uint64(outs[s, k]) + uint64(1), PRS, INV, LIM): return False
    return True

@njit(cache=True)
def process(lo, hi, roots, hits, hpos, state, mate, outs, nout, freel, PRS, INV, LIM, W, CAP, stats):
    """insert all vertices in [lo, hi).  state = [ncol, nfree].  returns (status, hpos)"""
    n = hi - lo
    mark = np.zeros(n, dtype=np.bool_)
    for i in range(roots.shape[0]):
        p2 = roots[i, 0]
        for k in (1, 2):
            R = roots[i, k]
            st = R if R >= lo else R + ((lo - R + p2 - 1) // p2) * p2
            for x in range(st, hi, p2): mark[x - lo] = True
    while hpos < hits.shape[0] and hits[hpos] < hi:
        if hits[hpos] >= lo: mark[hits[hpos] - lo] = True
        hpos += 1
    MAXQ = 2048
    qel = np.zeros(MAXQ, dtype=np.int64); qpar = np.zeros(MAXQ, dtype=np.int64); qcls = np.zeros(MAXQ, dtype=np.int64)
    for j in range(n):
        if not mark[j]: continue
        x = lo + j
        stats[0] += 1
        ncol = state[0]
        r = x % 25
        if r == 7:
            s = ncol % W
            mate[s] = -1; nout[s] = 0
            freel[state[1]] = ncol; state[1] += 1
            state[0] += 1
            continue
        cmin = ncol - W if ncol > W else 0
        if r == 18:
            head = 0; tail = 1
            qel[0] = x; qpar[0] = -1; qcls[0] = -1
            done = -1; donecls = -1
            while head < tail and done < 0:
                e = qel[head]
                for t in range(state[1] + 24):
                    if t < state[1]:
                        c = freel[t]
                        if c < cmin or mate[c % W] >= 0: continue
                    else:
                        c = cmin + np.random.randint(0, ncol - cmin)
                    if c == qcls[head]: continue
                    stats[1] += 1
                    if not compat_free(e, c, W, mate, outs, nout, PRS, INV, LIM, True): continue
                    if mate[c % W] < 0:
                        done = head; donecls = c; break
                    on = False; h = head
                    while h >= 0:
                        if qcls[h] == c: on = True; break
                        h = qpar[h]
                    if on: continue
                    if tail < MAXQ:
                        qel[tail] = mate[c % W]; qpar[tail] = head; qcls[tail] = c; tail += 1
                    if head == 0 and tail >= 6: break   # enough children; go deeper
                head += 1
            if done < 0:
                return -x, hpos
            h = done; c = donecls
            while h >= 0:
                mate[c % W] = qel[h]
                c = qcls[h]; h = qpar[h]
            stats[2] += 1
            k = 0
            for t in range(state[1]):
                if freel[t] >= cmin and mate[freel[t] % W] < 0:
                    freel[k] = freel[t]; k += 1
            state[1] = k
            continue
        placed = False
        for t in range(512):
            c = cmin + np.random.randint(0, ncol - cmin)
            s = c % W
            if nout[s] >= CAP: continue
            stats[3] += 1
            if compat_free(x, c, W, mate, outs, nout, PRS, INV, LIM, False):
                outs[s, nout[s]] = x; nout[s] += 1; placed = True; break
        if not placed:
            stats[4] += 1
            for c in range(ncol - 1, cmin - 1, -1):
                s = c % W
                if nout[s] >= CAP: continue
                if compat_free(x, c, W, mate, outs, nout, PRS, INV, LIM, False):
                    outs[s, nout[s]] = x; nout[s] += 1; placed = True; break
        if not placed:
            return -x, hpos
    return 1, hpos

if __name__ == '__main__':
    np.random.seed(SEED)
    t0 = time.time()
    sq = math.isqrt(NMAX) + 1
    small = primes_upto(sq).astype(np.int64)
    roots = small_roots(small)
    hits = large_hits(NMAX, sq, small)
    print(f"N={NMAX}: {roots.shape[0]} small roots, {len(hits)} large hits, prep {time.time()-t0:.1f}s", flush=True)
    W = 1 << 20; CAP = 8
    mate = -np.ones(W, dtype=np.int64); outs = np.zeros((W, CAP), dtype=np.int64); nout = np.zeros(W, dtype=np.int64)
    freel = np.zeros(1024, dtype=np.int64)
    state = np.zeros(2, dtype=np.int64)
    stats = np.zeros(8, dtype=np.int64)
    hpos = 0
    SEG = 10**6
    lo = 1
    nextrep = 10**7
    @njit
    def _seed(s): np.random.seed(s)
    _seed(SEED)
    while lo <= NMAX:
        hi = min(NMAX + 1, lo + SEG)
        st, hpos = process(lo, hi, roots, hits, hpos, state, mate, outs, nout, freel, PRS, INV, LIM, W, CAP, stats)
        if st < 0:
            print(f"FAILED at x={-st}; ncol={state[0]}; stats={list(stats)}", flush=True); sys.exit(2)
        assert state[0] == (hi - 1 + 18) // 25, (state[0], hi - 1)
        lo = hi
        if hi - 1 >= nextrep or hi > NMAX:
            print(f"OK through N={hi-1}: ncol={state[0]}=M(N); |V|={stats[0]}; tests18={stats[1]} aug={stats[2]} testsOut={stats[3]} scans={stats[4]}; {time.time()-t0:.0f}s", flush=True)
            nextrep += 10**7
    print("CERTIFIED: omega(G_N) <= M(N) for all N <=", NMAX, flush=True)
    if len(sys.argv) > 3 and sys.argv[3] == 'dump':
        assert state[0] <= W
        with open(f'dump_{NMAX}.txt', 'w') as f:
            for c in range(state[0]):
                s_ = c % W
                f.write(' '.join(map(str, [7 + 25 * c, int(mate[s_])] + [int(v) for v in outs[s_, :nout[s_]]])) + chr(10))
        print('dumped', state[0], 'classes')

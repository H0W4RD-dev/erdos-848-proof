"""Erdős #848, small-N side: an evolving "squarefree-compatible" cover, logged for replay.

Notation. V(N) = {a <= N : a^2+1 is NOT squarefree}. Admissible A (ab+1 never squarefree,
a=b allowed) satisfy A ⊆ V(N), and A is a clique for "xy+1 not squarefree".
Call x,y COMPATIBLE if xy+1 IS squarefree.

Certificate. If V(N) is partitioned into g groups of pairwise-compatible numbers, then
|A ∩ group| <= 1, so |A| <= g. The bound B(N) = floor((N+18)/25) is constant on the segments
[25t+7, 25t+31]. Restricting a valid partition of V(25t+31) to V(N), for N in that segment,
leaves at most the same number of groups, so it is enough to have a partition of V(25t+31)
into <= t+1 groups for every t >= 0. (V(N) is empty for N <= 6, where B(N)=0.)

We process V in increasing order, keeping a partition and logging every operation:
  N g        new group g containing the vertex (only allowed when a new segment begins, i.e.
             the vertex is ≡ 7 mod 25; such vertices are pairwise incompatible anyway)
  A v g      add vertex v to group g
  M u g h    move u from group g to group h
  E t        end of segment t (verifier checks #nonempty groups <= t+1 here)
The verifier (verify_cover.py) replays the log with an independent squarefree test.

Usage: py -3.13 build_cover.py M [log.txt]
"""
import sys, math, random
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def vertex_set(M, P):
    """a in [1,M] with p^2 | a^2+1 for some prime p (only p ≡ 1 mod 4 can occur)."""
    mark = np.zeros(M + 1, dtype=bool)
    for p in P:
        p = int(p)
        if p % 4 != 1: continue
        q = p * p
        if q > M * M + 1: break
        g = 2
        while pow(g, (p - 1) // 2, p) != p - 1: g += 1
        r = pow(g, (p - 1) // 4, p)
        r = (r - (r * r + 1) * pow(2 * r, -1, q)) % q   # Hensel lift
        for root in (r, q - r):
            assert (root * root + 1) % q == 0
            if root <= M: mark[root::q] = True
    mark[0] = False
    return [int(x) for x in np.nonzero(mark)[0]]

class SqfreeTester:
    """Exact test for n < limit: no p^2 | n for p <= B = cbrt(limit); then the cofactor m
    (all primes <= B removed) has at most two prime factors, so n is squarefree iff m is
    not a perfect square > 1."""
    def __init__(self, limit):
        B = int(round(limit ** (1 / 3))) + 2
        while B ** 3 < limit: B += 1
        self.p = primes_upto(B).astype(np.int64)
        self.p2 = self.p * self.p
        self.limit = limit
        self.cache = {}
    def sqf(self, n):
        c = self.cache.get(n)
        if c is not None: return c
        assert n < self.limit
        if np.any(n % self.p2 == 0): r = False
        else:
            m = n
            for d in self.p[n % self.p == 0]: m //= int(d)
            s = math.isqrt(m); r = not (m > 1 and s * s == m)
        if len(self.cache) < 5_000_000: self.cache[n] = r
        return r
    def compatible(self, v, members, skip=None):
        return all(self.sqf(v * x + 1) for x in members if x != skip)
    def conflicts(self, v, members, cap=2):
        out = []
        for x in members:
            if not self.sqf(v * x + 1):
                out.append(x)
                if len(out) >= cap: break
        return out

def build(M, log, seed=1):
    rng = random.Random(seed)
    V = vertex_set(M, primes_upto(M + 1))
    T = SqfreeTester(M * M + 2)
    groups = []                     # list of member lists
    where = {}                      # vertex -> group id
    seg = -1
    repairs = 0
    def end_segments_before(v):
        nonlocal seg
        t_now = (v - 7) // 25       # segment containing v
        while seg < t_now - 1:
            seg += 1
            assert sum(1 for g in groups if g) <= seg + 1
            log.write(f"E {seg}\n")
    for v in V:
        end_segments_before(v)
        if v % 25 == 7:
            groups.append([v]); where[v] = len(groups) - 1
            log.write(f"N {len(groups)-1} {v}\n"); continue
        # 1) direct placement: recent groups, then random sample, then all
        order = list(range(len(groups) - 1, max(-1, len(groups) - 60), -1))
        order += rng.sample(range(len(groups)), min(60, len(groups)))
        placed = False
        for gi in order:
            if T.compatible(v, groups[gi]):
                groups[gi].append(v); where[v] = gi
                log.write(f"A {v} {gi}\n"); placed = True; break
        if not placed:
            for gi in range(len(groups)):
                if T.compatible(v, groups[gi]):
                    groups[gi].append(v); where[v] = gi
                    log.write(f"A {v} {gi}\n"); placed = True; break
        if placed: continue
        # 2) one-move repair: a group where v conflicts with exactly one non-anchor member u,
        #    and u fits into some other group h
        for gi in range(len(groups)):
            c = T.conflicts(v, groups[gi], cap=2)
            if len(c) != 1: continue
            u = c[0]
            for hi in range(len(groups)):
                if hi == gi: continue
                if T.compatible(u, groups[hi]):
                    groups[gi].remove(u); groups[hi].append(u); where[u] = hi
                    log.write(f"M {u} {gi} {hi}\n")
                    groups[gi].append(v); where[v] = gi
                    log.write(f"A {v} {gi}\n"); placed = True; repairs += 1; break
            if placed: break
        if not placed:
            raise RuntimeError(f"could not place vertex {v} (groups={len(groups)})")
    # close the final segment only if it is complete within [1, M]
    t_last = (M - 31) // 25
    while seg < t_last:
        seg += 1
        assert sum(1 for g in groups if g) <= seg + 1
        log.write(f"E {seg}\n")
    return V, groups, repairs, seg

if __name__ == "__main__":
    M = int(float(sys.argv[1]))
    out = sys.argv[2] if len(sys.argv) > 2 else f"log_{M}.txt"
    with open(out, "w") as log:
        log.write(f"# Erdos 848 cover log, M={M}\n")
        V, groups, repairs, seg = build(M, log)
    sizes = [len(g) for g in groups if g]
    print(f"M={M} |V|={len(V)} groups={len(sizes)} repairs={repairs} last_segment={seg} "
          f"(covers every N <= {25*seg+31}) max_group={max(sizes)} -> {out}")

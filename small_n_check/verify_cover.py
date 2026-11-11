"""Independent verifier for the Erdős #848 cover logs (build_cover.py).

It shares no code with the builder:
  * V(M) is recomputed by testing a^2+1 directly (a full factorisation by trial division
    with a precomputed prime list), not by the Hensel-root sieve;
  * squarefreeness of x*y+1 is tested by sympy.factorint, not the cube-root trick;
  * the log is replayed op by op, checking
      - vertices arrive exactly once each, in increasing order, and they are exactly V(M);
      - a new group is opened only by a vertex ≡ 7 (mod 25);
      - every add or move keeps the target group pairwise compatible (xy+1 squarefree);
      - at "E t" all vertices <= 25t+31 have arrived and #nonempty groups <= t+1.
Conclusion when it passes: for every N <= 25*t_last+31, any A ⊆ [N] with ab+1 never
squarefree (a=b allowed) has |A| <= floor((N+18)/25).

Usage: py -3.13 verify_cover.py log.txt M
"""
import sys
from sympy import factorint

def is_squarefree(n):
    return all(e == 1 for e in factorint(n).values())

def main(path, M):
    V = [a for a in range(1, M + 1) if not is_squarefree(a * a + 1)]
    Vset = set(V)
    groups = {}
    where = {}
    arrived = []
    last_t = -1
    checked_pairs = 0
    def add(v, g):
        nonlocal checked_pairs
        for x in groups[g]:
            assert is_squarefree(v * x + 1), f"incompatible {v},{x} in group {g}"
            checked_pairs += 1
        groups[g].append(v); where[v] = g
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            op, *a = line.split()
            a = list(map(int, a))
            if op == "N":
                g, v = a
                assert g not in groups and v % 25 == 7 and v in Vset and v not in where
                assert not arrived or v > arrived[-1]
                groups[g] = [v]; where[v] = g; arrived.append(v)
            elif op == "A":
                v, g = a
                assert v in Vset and v not in where and (not arrived or v > arrived[-1])
                add(v, g); arrived.append(v)
            elif op == "M":
                u, g, h = a
                assert where.get(u) == g and g != h and h in groups
                groups[g].remove(u); del where[u]
                add(u, h)
            elif op == "E":
                t, = a
                assert t == last_t + 1
                cut = 25 * t + 31
                need = [x for x in V if x <= cut]
                assert arrived[:len(need)] == need, f"segment {t}: missing vertices"
                assert all(x > cut for x in arrived[len(need):])
                ng = sum(1 for g in groups.values() if g)
                assert ng <= t + 1, f"segment {t}: {ng} groups > {t+1}"
                last_t = t
            else:
                raise ValueError(line)
    assert sorted(arrived) == V, "log does not cover exactly V(M)"
    print(f"OK: |V(M)|={len(V)}, groups={sum(1 for g in groups.values() if g)}, "
          f"pairs checked={checked_pairs}, bound proved for all N <= {25*last_t+31}")

if __name__ == "__main__":
    main(sys.argv[1], int(float(sys.argv[2])))

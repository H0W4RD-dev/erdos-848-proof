"""Exact small-N check for Erdos 848.
For each vertex v (v <= NMAX, mu(v^2+1)=0) verify omega(earlier neighbours of v) <= M(v)-1,
which is equivalent to omega(G_N) <= M(N) for all N <= NMAX.
Max clique via bitset branch and bound with greedy colouring bound (Tomita MCQ style),
searching only for a clique of size >= target."""
import sys, time
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000

def nonsqfree_sieve(L):
    ns = np.zeros(L + 1, dtype=bool)
    p = 2
    sieve = np.ones(int(L**0.5) + 2, dtype=bool); sieve[:2] = False
    for i in range(2, len(sieve)):
        if sieve[i]:
            sieve[i*i::i] = False
            ns[i*i::i*i] = True
    return ns

L = NMAX * NMAX + 1
NS = nonsqfree_sieve(L)
M = lambda n: (n + 18) // 25
V = [a for a in range(1, NMAX + 1) if NS[a*a + 1]]
idx = {a: i for i, a in enumerate(V)}
n = len(V)
adj = [0] * n
for i, a in enumerate(V):
    m = 0
    for j, b in enumerate(V):
        if j != i and NS[a*b + 1]:
            m |= 1 << j
    adj[i] = m

def colour_order(P):
    """greedy colouring of vertex set P (bitset); returns list of (v, colour) in increasing colour."""
    order = []
    col = 0
    U = P
    while U:
        col += 1
        Q = U
        while Q:
            low = Q & -Q
            v = low.bit_length() - 1
            Q &= ~low
            Q &= ~adj[v]
            U &= ~low
            order.append((v, col))
    return order

class Found(Exception):
    pass

nodes = 0
def expand(size, P, target):
    global nodes
    nodes += 1
    order = colour_order(P)
    for v, c in reversed(order):
        if size + c < target:
            return
        newP = P & adj[v]
        if size + 1 >= target:
            raise Found
        if newP:
            expand(size + 1, newP, target)
        P &= ~(1 << v)

def has_clique(P, target):
    if target <= 0:
        return True
    try:
        expand(0, P, target)
    except Found:
        return True
    return False

t0 = time.time()
worst = []
for i, v in enumerate(V):
    earlier = adj[i] & ((1 << i) - 1)
    k = M(v)  # clique of size k in earlier nbrs would give clique of size k+1 > M(v)
    if has_clique(earlier, k):
        print("VIOLATION at v =", v, "M(v) =", k, flush=True)
        sys.exit(1)
print(f"OK: omega(G_N) <= M(N) for all N <= {NMAX}; |V|={n}; nodes={nodes}; {time.time()-t0:.1f}s")

# positive control: the 7-class clique of size M(v)-1 must be found among earlier neighbours
bad = 0
for i, v in enumerate(V):
    if v % 25 == 7 and v > 7:
        earlier = adj[i] & ((1 << i) - 1)
        if not has_clique(earlier, M(v) - 1):
            bad += 1
print("positive control failures:", bad)
# brute-force cross-check for tiny N with networkx
if NMAX <= 400:
    import networkx as nx
    for N in range(1, NMAX + 1):
        G = nx.Graph(); W = [a for a in V if a <= N]; G.add_nodes_from(W)
        for x in range(len(W)):
            for y in range(x + 1, len(W)):
                if NS[W[x]*W[y] + 1]: G.add_edge(W[x], W[y])
        w = max((len(c) for c in nx.find_cliques(G)), default=0)
        assert w <= M(N), (N, w)
        if N % 25 == 7: assert w == M(N), (N, w, M(N))
    print("networkx cross-check OK to", NMAX)

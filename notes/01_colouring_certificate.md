# 01: Colouring certificates for small and medium N

**Graph.** V = {x ≥ 1 : μ(x²+1) = 0}. Join x ≠ y iff μ(xy+1) = 0.
Valid sets A ⊆ [N] are exactly the cliques of G_N = G[V ∩ [N]] (the condition at a = b is x ∈ V).

**Lemma C1 (PROVED, trivial).** ω(G_N) ≤ χ(G_N).
So if G_N has a proper colouring with M(N) = ⌊(N+18)/25⌋ colours, then the conjecture holds at N.

**Lemma C2 (PROVED).** Suppose we insert the vertices of V in increasing order, keep a proper colouring
after each insertion, and the colour count after inserting x is ≤ M(x). Then the conjecture holds for all N ≤ the last vertex.
Proof: for N between consecutive vertices x ≤ N < x', the current state is a proper colouring of G_N with ≤ M(x) ≤ M(N) colours.
Recolouring earlier vertices is allowed; each N needs only its own colouring.

**Structure (PROVED).**
- The 7-class and the 18-class mod 25 are each cliques (25 | xy+1).
- With exactly M(N) colours, every colour class contains exactly one 7-class element (its "anchor").
- At N = 18+25j the 18-class elements below N number j+1, the same as the anchors.
  So the 18-class must be **perfectly matched** to anchors. Each 18-element b is paired with an anchor a such that ab+1 is squarefree.
  This forces augmenting-path repairs (first failure of naive first-fit: N = 43).
- An outsider (x ≢ ±7 mod 25) needs a class whose anchor, 18-mate and other outsiders all have squarefree xy+1.

**Results so far (COMPUTER)**
- small_clique.py: exact max-clique check ω(G_N) ≤ M(N) for all N ≤ 5000, cross-checked with networkx find_cliques for N ≤ 400.
  The positive control (finding the (M−1)-clique) passes.
- colour_cert.py (anchored classes plus 1-step repair, exact squarefree test by trial division to the cube root, then a perfect-square test of the cofactor):
  proper M(N)-colouring for all N ≤ 1e5. That run used 1380 repair moves. It is O(N²) because of full scans; a matching-based rewrite is pending.

**Squarefree test (PROVED).** Let n ≥ 1. Trial-divide by primes p in increasing order while p³ ≤ m, where m is the current cofactor
(divide p out once; if p still divides, n is not squarefree). At exit every prime factor of m exceeds the last p tried, and m < (that p)³.
So m ∈ {1, q, q·q', q²}, and n is squarefree iff m is not a perfect square.

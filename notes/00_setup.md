# Erdős #848 (JSP-000701): notes 00, setup

**Problem.** For N ≥ 1, let A ⊆ [N] = {1,...,N} satisfy μ(ab+1) = 0 for all a, b ∈ A (a = b allowed).
Is |A| ≤ M(N) := #{n ≤ N : n ≡ 7 (mod 25)} = ⌊(N+18)/25⌋ for **every** N?
The original source [Er92b, p.239] has no "sufficiently large".

**Status (2026-11-10)**
- CITED: Sawhney (note at math.columbia.edu/~msawhney/Problem_848.pdf) proves it for N ≥ N0, with N0 ineffective.
- CITED (unrefereed, AI-generated): N. Sothanaphan, "An Explicit Threshold in Erdős Problem #848" (Mar 24 2026), proves N0 = 2.64e17. Saved copy: src/g_1ujhm4...txt.
  - Structure: A7, A18, and A* = everything else.
  - Case 1: A* contains an even element.
  - Case 2a/2b: A* ≠ ∅ and all of A* is odd.
  - Case 3: A* = ∅ and both A7, A18 are nonempty.
  - Inputs:
    - Cor 1: single-witness squarefree count, error κ N^{3/4}.
    - Cor 2: |A*| ≤ #{x ≢ ±7 (mod 25) : μ(x²+1) = 0} ≤ (23/25) C_quad N + 0.001973 N.
  - The binding case is Case 1, at 0.0376 + errors against 0.04.
- Small N: Alexeev reports checks to a few hundred. A forum claim of N ≤ 1e7 (hjyuh/MalekZ) was flagged as likely incorrect.
- **Gap:** N from about 1e3 to 2.64e17 is uncovered.

**Basic facts (PROVED, elementary)**
- F1. If a ∈ A, then μ(a²+1) = 0. So p² | a²+1 for some prime p ≡ 1 (mod 4).
- F2. 25 | a²+1 ⟺ a ≡ ±7 (mod 25). 169 | a²+1 ⟺ a ≡ ±70 (mod 169).
- F3. For a ≡ 7 and b ≡ 18 (mod 25): ab+1 ≡ 2 (mod 25), so no cross pair uses 5.
  For b ∉ {7, 18} (mod 25) and a ≡ ±7: 25 ∤ ab+1.
- F4. The classes 7 and 18 mod 25 are cliques. M(N) = |class 7 ∩ [N]| ≥ |class 18 ∩ [N]|.

**Plan**
1. Exact small-N verification (max clique), as far as possible.
2. Audit Sothanaphan's proof and numerics.
3. Lower the analytic threshold. Ideas:
   - (i) least-witness slack: A* ⊆ [b, N].
   - (ii) joint constraint: a ∈ A* also needs μ(ab+1) = 0 for the witness b.
   - (iii) sharper root counts using the true ω(b) range.
4. Record whatever gap remains as the barrier.

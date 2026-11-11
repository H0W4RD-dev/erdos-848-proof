# Erdős #848 (JSP-000701): a complete computer-assisted proof for every N

**Theorem.** For every N ≥ 1 and every A ⊆ [N] = {1,…,N} with ab+1 not squarefree for all a, b ∈ A (a = b allowed),
  |A| ≤ M(N) := #{n ≤ N : n ≡ 7 (mod 25)} = ⌊(N+18)/25⌋.
Equality holds for A = [N] ∩ (7 + 25ℤ), because 25 | ab+1 for all a, b ≡ 7 (mod 25).
So the answer to Erdős's question ("is this sequence maximal?") is **YES, for every N**.

Context:
- An independent, reproducible proof. See README.md, section "Relation to JSP PR #4667", for how it relates to the earlier claim in TheJustinSunPrize/awards#4667.
- Produced 2026-11-10 with substantial assistance from Anthropic's Claude (several Claude Code sessions that cross-audited each other); not yet human-refereed.
- File paths below refer to the original working folder; in this repository the certificate files are in certificate/, the analytic scripts and logs in analytic/, and the notes in notes/.

Status labels: PROVED (a proof is written here or in the cited note), COMPUTER (an exact or interval-certified computation), CITED (literature).

---

## 0. Notation and the reduction (PROVED)
- A7 = A ∩ (7 + 25ℤ), A18 = A ∩ (18 + 25ℤ), A* = A \ (A7 ∪ A18).
- Basic facts:
  - (F1) 25 | x²+1 ⟺ x ≡ ±7 (mod 25).
  - (F2) 25 | xy+1 with x ≡ 7 holds iff y ≡ 7; with x ≡ 18, iff y ≡ 18.
  - (F3) 4 | xy+1 ⟺ x, y are odd and x ≢ y (mod 4).
- If A ⊆ A7, then |A| ≤ M(N). If A ⊆ A18, then |A| ≤ ⌊(N+7)/25⌋ ≤ M(N).
- Otherwise it suffices to show |A| < N/25 − 1, since M(N) ≥ (N−6)/25.

## 1. Small N: colouring certificate (COMPUTER; method PROVED in 01_colouring_certificate.md)
- **Graph.** V = {x : μ(x²+1) = 0}. Join x ≠ y iff μ(xy+1) = 0. Valid sets A ⊆ [N] are exactly the cliques of G_N = G[V ∩ [N]].
- **Lemma C.** If G_N has a proper colouring with M(N) colours, then |A| ≤ M(N).
- **Construction.** Vertices are inserted in increasing order, keeping a proper colouring with exactly M(x) classes after each insertion:
  - every 7-class vertex is a new class's anchor;
  - 18-class vertices are kept perfectly matched to anchors by augmenting paths;
  - outsiders join any class with which they are pairwise compatible.
  For N between consecutive vertices, the current state colours G_N with M(N) colours.
- **Certified range: 1 ≤ N ≤ 2·10⁷, independently replayed.**
  - Generator: colour4_log.py (numba, exact squarefree test, bounds-checked) wrote the op log ops_20000000.bin.
  - Verifier: replay_verify2.py shares no code with the generator. It uses sympy.factorint for squarefreeness and sympy.sqrt_mod to rebuild V.
  - After every one of the 2,103,229 insertions it checks four things:
    1. every touched class is pairwise compatible;
    2. mates form a bijection onto the 18-class vertices ≤ x;
    3. #classes = M(x);
    4. the insertion order is exactly V.
  - Log: replay_2e7.log, "REPLAY VERIFIED N=20000000".
  - The verifier was mutation-tested: corrupted logs are rejected.
- Additionally, colour4.py certified N ≤ 7.1·10⁸ with in-code checks only (cert_1e9.log). It is **not needed** for the chain.

## 2. Lemmas for the analytic ranges (PROVED; full proofs in 02_lemmas.md and threshold/01_sharpened_lemmas.md; audited in 04_audit_threshold_S123.md)
- **R.** For b = 2^v m with m odd, #{y mod b : y² ≡ u} ≤ c_v ∏_{p|m}(1 + (u/p)), where (c_0, …, c_3) = (1, 1, 2, 4).
- **S2 (large-square tail).** #{a ≤ N : ∃ r > P1, r² | ab+1} ≤ ρ_st (2(N+1)/P1 + (N²+1)/P1²).
  A character-refined form, using J ⊆ {small primes of b}, replaces 2 by (1 − 1/p) per prime of J.
  Proof: write ab+1 = k r² and k r² ≡ 1 (mod b), so r lies in ≤ ρ_k classes mod b, where ρ_k vanishes unless k is a QR mod every p ∈ J.
- **S1 (one witness).** For a witness b and classes on which 4 ∤ ab+1 and 25 ∤ ab+1:
  - primes ≤ 29 are handled by inclusion–exclusion inside each class;
  - 29 < p ≤ Qs gives one class mod qp² per class;
  - **Qs < p ≤ P1 is counted once in total**, because {a : p² | ab+1} is a single class mod p²;
  - p > P1 is handled by S2.
- **Q / E1 / Pell.** Outsiders x need p² | x²+1 with p ≡ 1 (mod 4), p ≥ 13.
  - Large p are counted through x² + 1 = d u².
  - The admissible d < D (3 ∤ d, 4 ∤ d) number at most D/2 + 1.
  - Each d has ≤ 1 + log N / log(3+2√2) solutions with x ≤ N.
- **J / C1′ (Case 1, b ∈ A* even).** The joint count of x with μ(x²+1) = μ(bx+1) = 0. This includes the genuine root-class term: x ≡ b (mod p²) when p² | b²+1.
- **Case analysis (S3).** Assume A ⊄ A7, A18.
  - C1: A* has an even element.
  - C2: A* ≠ ∅ is all odd, split by whether A* uses both residues mod 4.
  - C3: A* = ∅, with sub-cases by parity and residues mod 4.
  - Each sub-case bounds |A| by a sum of S1/Q/J bounds (explicit list in threshold/01 §S3).
- **Monotone bound.** Every term is nondecreasing in N, and parameters range over fixed grids. So B(N1) < N0/25 − 1 certifies all N ∈ [N0, N1].

## 3. Proof chain (every range certified; overlaps shown)
| Range | Certified by | Arithmetic | Status |
| --- | --- | --- | --- |
| 1 ≤ N ≤ 2·10⁷ | §1 colouring certificate; replay_verify2.py (replay_2e7.log) | integers, sympy | **DONE** |
| 1.6·10⁷ ≤ N ≤ 1.387721·10⁹ | B3 (bound3.py), 40 intervals; exact: threshold/verify_exact.py → verify_exact_part1/2.log "ALL INTERVALS VERIFIED EXACTLY" | exact rationals | **DONE** (verifier audited in 04; spot re-check below) |
| 1.3·10⁹ ≤ N ≤ 1.219·10¹⁹ | B3L (bound3L.py), 101 intervals (coverL_13e8_12e19.log); exact: threshold/verify_exactL.py → verify_exactL.log "ALL INTERVALS VERIFIED EXACTLY" (101/101 OK, contiguous) | exact rationals | **DONE** |
| N ≥ 10¹⁹ | closed form αN + β(N), β(N)/N decreasing (03_analytic_range.md); tail_bound.py; interval-certified by threshold/tail_exact.py "CERTIFIED OK" | mpmath intervals + exact prime sums | **DONE** |

The overlaps are 1.6·10⁷ < 2·10⁷, 1.3·10⁹ < 1.3877·10⁹ and 10¹⁹ < 1.219·10¹⁹. **Every N ≥ 1 is covered.** ∎

The exact logs were checked for tiling: the 40 intervals of [1.6e7, 1.387721e9] are contiguous; 16 of them were verified in part1 and resumed in part2. The 101 intervals of [1.3e9, 1.218524e19] are contiguous with no failures.

Independent spot re-check (spot_exact.py; my own exact implementation of S1–S3, sharing no code with threshold/):
- [1.6e7, 1.616e7]: 0.03974, OK.
- [1.826993e11, 2.512115e11]: 0.03861, OK.
- Their tightest interval [5.473266e8, 7.29738e8] gives 0.04002 with my coarser grids. Re-covered as 4 sub-intervals of ratio 1.08, all OK (≤ 0.0327).
- So the range is independently covered (spot_exact.log; 04 §spot).

## 4. Cited inputs
- Dusart (2010): π(x) ≤ (x/ln x)(1 + 1/ln x + 2.51/ln²x) for x ≥ 355991.
- Robin (1983): ω(n) ≤ 1.38402 ln n / ln ln n for n ≥ 3. Used only for N ≥ 10¹⁹.
- Elementary facts: square roots modulo prime powers, and the negative Pell equation (proved in 02).
- Asymptotic idea: Sawhney's mod-25 decomposition (Problem_848.pdf). The structure of the explicit-threshold method follows Sothanaphan's note; it was audited, and no error was found.

## 5. Audit trail
- 04_audit_threshold_S123.md: Lemmas S1–S3, C1′ and S5 **CORRECT**. bound3/bound3L, verify_exact(L) and tail_exact **CORRECT** (code review plus brute-force checks).
- test_sqfree.py: the numba squarefree routine agrees with sympy on 5700 adversarial inputs.
- lemma_check.py, j_check.py, audit_S12.py, audit_roots.py: lemma bounds dominate brute-force counts; root counts are exact.
- **Remaining caveat:** the large computations are machine-checked (exact or interval arithmetic, independently replayed or re-implemented where stated) but not human-refereed.

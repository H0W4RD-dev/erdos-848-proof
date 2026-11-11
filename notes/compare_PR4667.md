# Erdős #848: mathematical comparison of PR TheJustinSunPrize/awards#4667 with our independent proof

Written 2026-11-10 by the threshold-side session (c7f777-thr), for the parent session. **Nothing has been posted.**

**Source.** github.com/q1lin570/JSP-lean @ 2cfd2281. I fetched these files read-only through the GitHub API into `pr4667_src/`:
- solution_overview.md, solution.md, closing_proof.md, threshold_2p28.md, threshold_2p41.md;
- completion_audit.md, check_constants.py, and the JSON results.

I ran **none** of their code. jsp-elementary-nt is auditing their Lean separately.
Their notes say they were produced by "Codex GPT-6-Astra" and claim no independent review.

Our side: `FULL_PROOF.md`, `threshold/PROOF_CHAIN.md`, `threshold/01_sharpened_lemmas.md`, `02_lemmas.md`.

## 1. Range decomposition

| PR #4667 | Range | Method | Ours | Range |
| --- | --- | --- | --- | --- |
| Finite certificate | 1 ≤ N ≤ 2^25 = 33,554,432 | Prefix colouring, at most ⌊(N+18)/25⌋ colours; C++ generator plus an independent C++ checker | Colouring certificate (colour4) | 1 ≤ N ≤ 7.1e8 generated; independently replayed (sympy only) to 2e7 |
| "Bridge" | [2^25, 2^28] | 5 cases, 3 dyadic intervals, 5,366 factor profiles, exact rationals | bound3 (S1–S3) | [1.6e7, 1.388e9], 40 intervals, exact rationals |
| "Intermediate" | [2^28, 2^41] | 4 cases, 13 dyadic intervals, 33,425 profiles, exact rationals | bound3L | [1.3e9, 1.22e19], 101 intervals, exact rationals |
| Closed form | N ≥ 2^41 ≈ 2.2e12 | Symbolic bounds, monotone, exact check at 2^41 | tail_bound | N ≥ 1e19, mpmath interval check at 1e19 |

Both are certificate + interval cover + closed-form tail. Their certificate is shorter (2^25 vs our 7.1e8), their analytic range starts at 2^25 ≈ 3.4e7 (ours at 1.6e7), and their closed-form tail starts much earlier (2.2e12 vs our 1e19).

## 2. Step-by-step mapping

| Step | PR #4667 | Ours | Same idea? | Problems found? |
| --- | --- | --- | --- | --- |
| Reduction to cliques and colouring | Clique ≤ χ; prefix colouring; colours ≤ ⌊(N+18)/25⌋ after each event | Lemmas C1, C2 of 01 | **Same** | None. **Independent check:** I recomputed #{x ≤ 2^25 : x²+1 not squarefree} with my own sieve and got **3,528,561, exactly their vertex count** |
| Squarefree test in the certificate | Trial division to the cube root, then a square test on the cofactor | Same cube-root lemma | **Same** | None |
| Case split (A7, A18, A*) | Sawhney/Sothanaphan 4 cases, refined in the bridge to 5 (I, IIa, IIb, III, IV) | C1, C2(i), C2(ii), C3 (S3) | **Same idea.** Their IIa ("A* meets both odd classes mod 4, two witnesses") is our C2(i). Their IIb ("A* in one class mod 4", main term d/4) is our Q_u case | None; exhaustive. Our C2(ii) goes further (P7/P18 control each other, so no 2·N/100 term); theirs pays 2(N/100+1) in IIb, which costs them nothing since IIb is not binding |
| Small-prime off-diagonal factors | ≥ 2^41: union bound Σ1/p² (C < 0.1624). Bridge: exact inclusion–exclusion over the allowed primes < 19, error 2^{\|F\|} | Exact inclusion–exclusion over primes ≤ 29 (02 Lemma W), error 2^k − 1 | **Same** (theirs uses fewer primes) | None |
| ρ(b) | ≥ 2^41: ρ(b) ≤ 18 b^{1/5}. Middle ranges: factor profiles (2-adic part, J ⊆ {3..19}, k large primes), h = c_v 2^{\|J\|+k}, l = lower bound for b | ρ ≤ c_v 2^{\|T\|+w} over structures (v, T ⊆ {3..29}, w) | **Same** profile idea. 18 b^{1/5} is a Sothanaphan-style power bound | Checked: max of ρ/b^{1/5} = (4/8^{1/5}) ∏_{p≤31} 2/p^{1/5} = **17.04 < 18**. OK |
| Large off-diagonal tail | a·b+1 = m p², p in ρ(b) classes mod b: ρ(2X/(bR) + X/R²) | Lemma W tail (02) | **Same** | None |
| Shared tail per witness | "Restrictions using the same b share one large-prime tail" | **S1** | **Same idea**, found independently. They share only the q > R tail; we also share the middle prime range (Qs, P1] | None |
| Character-refined tail | δ = φ(c0)/(c0 ρ(c0)), κ = c0 δ: the cofactor m must be an invertible square mod c0, so U ≤ ρ[δ(...) + κ(2√X/b + 1)] (2^28 note, (6)–(8)) | **S2**: k must be a nonzero QR mod each p ∈ J ⊆ T, so a factor \|G\|/L, + \|G\| + √((N+1)/base)·Σg^{−1/2} | **Same idea**, found independently. Theirs also includes the **2-adic part** in c0 (squares mod 8: δ gains 1/8); ours uses odd primes only | None. I re-derived (6)–(8): Σ_{m in class} m^{−1/2} ≤ 2δ√M + 2κ holds |
| Unknown-prime refinement | Bridge (10): if b has k ≥ 2 unknown large primes, the smallest q ≤ (b/B)^{1/k} gives δ_{c0 q} ≤ δ0/2 | not used | **Different idea** (theirs only) | Checked: q depends on b, but each b uses its own q. The term √((N+1)/b)·b^{1/k} decreases in b for k ≥ 2, so b = l is the worst case. Valid |
| Joint constraint for A* (Case 1) | Bridge and intermediate: only p = 13, 17. One diagonal root class mod p² is compatible with b, the other needs a second square factor, saving G = Σ_{13,17}(1 − H_p)/p². ≥ 2^41: **none** (A* bounded by the diagonal count alone) | **C1′** (J-part): full joint (p, q) split for all p ≤ R0 = 1000, plus the "same" root-class term | **Same idea**; ours is much more extensive | None. This is the main structural difference in Case 1; see §3 |
| Diagonal count (A*) | Joint-root discrepancy: two classes ±r mod q, minus 2 classes mod 25, ≤ 46N/(25q) + 69/25 | Lemma Q (02): per class 2·23(N/(Qp²) + 1) | Same purpose; **their rounding is sharper** | Their lemma was checked by their script; **I brute-forced it** for all p ≡ 1 (4), 13 ≤ p < 400, q ∈ {p², 2p²}, over a full period: max excess 2.576 ≤ 69/25 = 2.76. OK |
| Large diagonal factors (Pell) | x² − dy² = −1. ≥ 2^41: all d ≤ Y, ratio > 3. Bridge: d ≡ 1, 2, 5 (mod 8), gcd(d, 21) = 1 (3Y/14 + 45/14), and "only the minimal solution can have prime y unless y0 = 1" | Lemma P (02): admissible d ≤ D/2 + 1, ratio 3+2√2; counts all integer u > R1 | **Same core.** Their bridge refinements (prime denominator, residues mod 168) are **different and stronger** | Checked: the mod-168 offsets are exactly 45/14 and 15/7 (attained, so "≤" is right). The prime-denominator argument: higher odd powers have y-coefficient a multiple of y0 greater than y0, so they are composite if y0 > 1. Correct |
| Prime counting | ≥ 2^41: elementary π(x) ≤ 3x/log x for x ≥ 2^16 (via θ and central binomials, proved). Middle: Rosser–Schoenfeld (3.2). π₁ bound checked on [2^15, 2^21] | Exact sieve to 6e7, then Dusart 2010 | Different citations, same role | Checked their θ-derivation: the constants work (θ < 3x/2 needs (ln²x/ln2 + ln x)/x < 0.114, which is 0.003 at 2^16) |
| Constants | D = Σ_{p≡1(4), p≥13} 2/p² < 0.0277; C = Σ_{p≠2,5} 1/p² < 0.1624 | Products 1 − ∏(1 − p⁻²) (exact finite) + Σ tails | Same quantities; theirs are sums (union bound), ours products where possible | **Recomputed:** D = 0.0276275 + (tail < 2e−8) < 0.0277 (margin 7.2e−5); C = 0.1622474 < 0.1624. OK |
| Monotonicity / interval logic | Profiles generated to the right endpoint, errors checked at the left endpoint (errors decrease) | Bound nondecreasing in N, evaluated at the right endpoint against N0/25 − 1 at the left | Same logic, mirrored | None |
| Exact arithmetic | Float search, then rational recheck with integer-sqrt intervals and rational log series | Same (verify_exact*, tail_exact) | **Same** | None |

## 3. The tight 0.04 margin: theirs vs ours
Both proofs have cases within ~1e−4 of 1/25 somewhere. Both handle this the same way: an exact rational (or outward-rounded) recheck, with a strict "< 1/25" test.

| Where | PR #4667 worst bound | Ours worst bound |
| --- | --- | --- |
| Start of the analytic range | Bridge [2^25, 2^26]: Case III **0.039522** | [1.6e7, 2.6e7]: **0.03994** (C1) |
| Middle | [2^28, 2^29]: Case 3 **0.039877** | [5.47e8, 7.30e8]: **0.039972** (tightest overall) |
| Large N | **N = 2^41, Case 1: 0.039926**, from constants and errors (independently reproduced: 0.0399263) | Large N: ~0.0297 (point values); the cover's worst 0.039961 comes only from wide intervals |
| Closed-form tail | 0.039926 at 2^41 (decreasing errors) | 0.031496 at 1e19 |

Why their 2^41 Case 1 is tight:
- With no joint constraint for A*, the main term is (23/25)D + (2/25)C = 0.038476 (rounded-up constants), against Sothanaphan's 0.037612 (products). The union-bound constants lose about 8.6e−4.
- The remaining errors at 2^41 are Q + 2P + T ≈ 3.85e−4 + 8.1e−4 + 2.6e−4 = 1.45e−3, giving 0.039926. With the true D and C it would be 0.039847.
- So their threshold of 2^41 is essentially forced by Case 1. At 2^40 the same formula gives 0.04027 and at 2^39 0.04069 (computed), so the bound fails just below 2^41.
- **This is not an error.** The margin of 7.4e−5 is real and checked exactly, and every input constant reproduces (D, C, 18, the error formulas).

Ours: the tight cases are intervals of the monotone cover (the bound at the right endpoint against N0/25 at the left), not a structural limit. At fixed N our worst case at large N is ~0.0297, because C1′ (the full joint constraint) removes the Case 1 bottleneck that both Sothanaphan and #4667's ≥ 2^41 step keep.

## 4. Summary verdict (mathematical, document level)
- **Architecture:** the same as ours, since both follow Sawhney/Sothanaphan. Several key refinements were **found independently by both**: shared tail per witness (S1), the square-residue restriction on the cofactor (S2), the mod-4 split with two odd witnesses (S3: IIa/IIb), and the root-class argument mod p² (their G for p = 13, 17; our C1′ for all p ≤ 1000).
- **Ideas only in #4667:**
  - 2-adic δ;
  - the smallest-unknown-prime refinement;
  - the prime-denominator and mod-168 Pell count;
  - the joint-root discrepancy rounding 46/25 + 69/25;
  - elementary prime counting.
- **Ideas only in ours:**
  - the full J-part for all p ≤ 1000 (so Case 1 asymptotically ≈ 0.027 vs their 0.0385);
  - sharing of the middle prime range;
  - P7/P18 mutual control.
- **Problems found in #4667's written mathematics:** none.
  - I spot-verified: D, C, the ρ ≤ 18b^{1/5} constant, the 2^41 Case 1 value, the mod-168 counts, the root-discrepancy lemma (brute force), and the vertex count to 2^25 (exact match).
  - I also re-derived (6)–(8) and (10) of their notes and found them correct.
- **Not checked by me:**
  - their Python and C++ scripts (deliberately not run);
  - the 33,425-profile and 5,366-profile computations themselves;
  - the Lean formalisation (jsp-elementary-nt's job);
  - their prose claim that "Sawhney's classification" is used only as the case split. I read it that way; no asymptotic theorem is used as an input.
- **Our proof's value relative to #4667:** an independent, different implementation of the same strategy with different parameter choices. Its tightest margins sit in different places, so the two cross-check each other. It also confirms their **answer**: |A| ≤ ⌊(N+18)/25⌋ for all N.

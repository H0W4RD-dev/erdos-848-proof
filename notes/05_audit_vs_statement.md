# Erdős #848: audit of FULL_PROOF.md against the literal statement (session c7f777, 2026-11-11)

> Note for this repository: paths such as `threshold/`, `lean/`, `leanprep/` and `src/` refer to the authors' working tree. In this repository the analytic scripts are in `analytic/`, the lemma notes in `notes/`, and the certificate in `certificate/`.


This is a fresh, case-by-case pass over every range and every cited result.

**Verdict.**
- There is **no logical gap.**
- I found 6 documentation-level slips, all repaired below (G1–G6), none of which changes a certified number.
- Two cited results remain: Dusart 2010 and Robin 1983. Both are used only for N > 6·10⁷ (Dusart) or N ≥ 10¹⁹ (Robin), and §5 shows how to remove both.

Reproduced today:
- `threshold/tail_exact.py 1e19` → `C1=0.031496 C2a=0.025934 C2b=0.029080 C3=0.020996 CERTIFIED OK`.
- The interval logs were re-parsed for tiling by my own script (a short ad-hoc script), with the results in §3.

## 0. The statement being proved
- **Source.** Erdős [Er92b, p.239]; Bloom's page #848. In Sothanaphan's quotation: ab+1 "never squarefree for all (not necessarily distinct) a,b ∈ A". The original has no "N sufficiently large".
  (erdosproblems.com returned 403 to my fetch; the wording was checked against the saved sources in `src/`.)
- **Formal version** (formal-conjectures, mirrored in `lean/E848Basic.lean`):
  - `Erdos848For N := ∀ A ⊆ range N, NonSqProd A → |A| ≤ #{n < N : n ≡ 7 (25)}`;
  - `NonSqProd A := ∀ a ∈ A, ∀ b ∈ A, ¬Squarefree (ab+1)`, so a = b is included.
- **Literal claim proved:** for every N ≥ 1 and every A ⊆ {1..N} with μ(ab+1) = 0 for all a, b ∈ A (a = b allowed), |A| ≤ M(N) := ⌊(N+18)/25⌋.
  - **M(N) is the class count.** #{1 ≤ n ≤ N : n ≡ 7 (25)} = ⌊(N−7)/25⌋ + 1 = ⌊(N+18)/25⌋ for N ≥ 7, and 0 = ⌊(N+18)/25⌋ for N ≤ 6. ✔
  - **The bound is attained.** 7·7 ≡ 49 ≡ −1 (mod 25), so 25 | ab+1 on the class 7 mod 25. ✔
  - **Formal form.** It follows at N−1 (`fc_of_upper`, sorry-free): 0 ∉ A because 0·0+1 = 1 is squarefree, and #{n < N : n ≡ 7} = ⌊(N+17)/25⌋. ✔
- **Not covered: the distinct-pairs variant** (a ≠ b only; Sothanaphan's separate note `src/g_15yo…txt`).
  The diagonal condition is used essentially: A ⊆ V = {x : μ(x²+1) = 0}, and the outsider counts (Lemma Q) rely on it.
  This is correct for the literal problem. It must not be advertised as covering the variant.

## 1. Reduction (FULL_PROOF §0): checked
- **Trivial cases.** A ⊆ A7 gives |A| ≤ M(N); A ⊆ A18 gives |A| ≤ ⌊(N+7)/25⌋ ≤ M(N); A = ∅ is trivial.
- **Otherwise,** |A| < N/25 − 1 ≤ (N−6)/25 ≤ M(N), because ⌊(N+18)/25⌋ ≥ (N+18−24)/25.
- **The three cases are exhaustive.** If A ⊄ A7 and A ⊄ A18, then A* ≠ ∅, or A* = ∅ with A7 ≠ ∅ ≠ A18:
  - C1: A* has an even element;
  - C2: A* ≠ ∅ and A* is all odd;
  - C3: A* = ∅. ✔

## 2. Range [1, 2·10⁷]: colouring certificate. Checked (verifier logic re-read line by line)
`replay_verify2.py`, run at N = 2·10⁷, ended "REPLAY VERIFIED … 2103229 insertions", EXIT 0.

1. **Vertex set is exact.** μ(x²+1) = 0 ⟺ p² | x²+1 for a prime p ≡ 1 (4): 4 ∤ x²+1, and q ∤ x²+1 for q ≡ 3 (4).
   The script marks x ≡ r (mod p²) for every root r of r² ≡ −1 (sympy `sqrt_mod`, all roots), over every prime 5 ≤ p ≤ N with p ≡ 1 (4). So every vertex x ≤ N is marked. ✔
2. **Every vertex ≤ x is coloured after inserting x.**
   - 7-class vertices are anchors of their own class.
   - Outsiders are appended to a class and never removed.
   - 18-class vertices are covered through the mates. Every touched mate has count 1, untouched counts are unchanged, all counts are ≤ 1, and #classes-with-mate = n18 = #18-vertices ≤ x. So the mate map is a bijection. ✔
3. **Classes are independent.** Each touched class is re-checked pairwise (xy+1 squarefree, `sympy.factorint`); untouched classes are unchanged. ✔
4. **Class count.** The script asserts #classes = ⌊(x+18)/25⌋ after each insertion.
   - Between consecutive vertices, M is nondecreasing, so this gives a proper M(N)-colouring of G_N for every N ≤ 2·10⁷. The loop runs to N itself.
   - N ≤ 6 has no vertices, and M = 0. ✔
5. A clique (that is, a valid A) meets each independent class at most once, so |A| ≤ M(N). ✔

## 3. Ranges [1.6·10⁷, 1.387721·10⁹] and [1.3·10⁹, 1.218524·10¹⁹]: exact covers. Checked
**Tiling.** I re-parsed every `exact B3/a = … OK` line:

| Logs | Intervals | Contiguous | Endpoints | Worst B3/a |
| --- | --- | --- | --- | --- |
| verify_exact_part1+2 | 40 | yes | 16 000 000 → 1 387 721 000 | 0.039972 |
| verify_exactL | 101 | yes | 1 300 000 000 → 12 185 240 000 000 000 000 | 0.039961 |

No FAIL lines appear, and both logs end "ALL INTERVALS VERIFIED EXACTLY".

**Lemmas re-derived again today** (as well as 04_audit_threshold_S123.md):
- **Lemma R:** Hensel lifting for odd p^e; mod 2^e there are 1, 2, 4 roots for e = 1, 2, ≥ 3; then CRT.
- **S1:** the inclusion–exclusion moduli ∏_{U}p² (p ∈ P0 ⊆ {3,7,…,29} \ T) are coprime to q | 100. An AP of length L hits each class ⌊L/m⌋ or ⌈L/m⌉ times, so the error is < 1 per term. "p² | ab+1, p ∤ b" is one class mod p².
- **S2 tail**, including the J-refinement: the map a ↦ (k, r) is injective, and ρ_k ≤ ρ_st·1[k mod L ∈ G].
- **Witness structure:** base = 2^v∏T | b; w primes > 29 give m ≥ ∏T·(31·37·…).
- **C1′:** every route for x ∈ A* is covered:
  - p > R0 → E1;
  - q = 5 → s_b;
  - q ∈ {3,7,11,13} \ {p} → inclusion–exclusion (oth);
  - 13 < q ≤ Qs, q ≠ p → pairs;
  - q = p → same: x ≡ b (p²), since x² ≡ −1 and bx ≡ −1 force b ≡ −x⁻¹ ≡ x;
  - q > Qs → the shared S1 range.
- **S3 sub-cases:** each checked against the residues mod 4 of the 8 classes mod 100:
  - A7: 7≡3, 32≡0, 57≡1, 82≡2;
  - A18: 18≡2, 43≡3, 68≡0, 93≡1.

**Monotonicity.** Each summand is nondecreasing in N, the parameter grids are fixed, and the structure set grows with N. The check B3(b_i) < a_i/25 − 1 covers all integers in [a_i, b_i]. ✔

## 4. Range N ≥ 10¹⁹: closed-form tail. Checked term by term (`threshold/tail_exact.py` vs 02_lemmas.md)
| Code term | Lemma | Check |
| --- | --- | --- |
| `Wc(q,P0)` | Lemma W | T = ∅ maximises both the density and 2^{\|P0\|}−1; Σ_{29<p≤P1}p⁻² ≤ SIG29; π(P1)−π(29) with π(29) = 10; (bN+1)/(bP1) ≤ (N+1)/P1, (bN+1)/P1² ≤ (N²+1)/P1² ✔ |
| `rho` = 4·N^{1.38402 ln2/lnln N} | Lemma R + Robin | ρ(b) ≤ c_v 2^{ω(odd part)} ≤ 4·2^{ω(b)}; see G2 for b ≤ N ✔ |
| `E1` | Lemma Q/E1 + Pell | 2N·Σ_{p>1000, p≡1(4)}p⁻² + 2π(R1) + (D/2+2)(1 + ln N/ln(2+√3)); (D/2+2) ≥ D/2+1 admissible d ✔; per-d count: G4 ✔ |
| `qodd` | Lemma Q, 23 odd classes mod 50 | ✔ (S50o has 23 elements, kernel-checked as `S50o_card`) |
| `jb`, `med`, `same`, `E2` | Lemma J | oth uses n_o = 22 (5 ∤ b) and +15 ≥ 2^{k_p}−1; med keeps the q = p pair (overcount); same over the full 𝒢; E2 = q > 1000 plus the ρ-tail ✔; the 5 \| b case: G1 |
| `c2a`, `c2b`, `c3` | 02 C2a/C2b/C3 | ✔; the C3 formula: G5 |
| target `(N/25 − 1).a` | reduction | the lower endpoint of the interval ✔ |

**Monotonicity (03).** Every N-dependent part of β(N)/N is one of:
- c/N;
- N^{−0.2}·h(0.8 ln N) with h(y) = (1 + 1/y + 2.51/y²)/y decreasing;
- ρ*N^{−c} for c ∈ {0.6, 0.8}, multiplied by the decreasing factors (1 + 1/N) and (1 + N^{−2});
- N^{−0.6}(1 + ln N/ln(2+√3)) and N^{−1}(1 + ln N/·).

So α + β(N)/N ≤ α + β(N*)/N* < 1/25 − 1/N* < 1/25 − 1/N for N > N*. ✔

## Gaps / slips found and repaired
- **G1 (tail, C1 with 5 | b, unstated).** `jb` = sb + 22·X is the 5 ∤ b form. For 5 | b the bound is 0 + 23·X.
  - Write X = Σ_p X_p with X_p = 2((N/(25p²)+1)δ_p + 15) and sb_p = 2(N/(25p²)+1).
  - Then X_p ≤ sb_p ⟺ (N/(25p²)+1)(1−δ_p) ≥ 15. For p ≤ 1000 and N ≥ 10¹⁹ the left side is ≥ 4·10¹¹·0.8.
  - So 23X ≤ sb + 22X, and the code's form dominates. **Repaired:** stated here, and the same remark was added to FULL_PROOF §2.
- **G2 (Robin applied to b ≤ N).** Robin bounds ω(b) in terms of b, not N. f(t) = ln t/lnln t is increasing for t ≥ e^e, so f(b) ≤ f(N) for 16 ≤ b ≤ N. For b ≤ 15, ω(b) ≤ 2 and ρ ≤ 16 ≤ ρ*(10¹⁹) ≈ 2.7·10⁵.
  - The constant is consistent: the primorial ratios k·lnln(p_k#)/ln(p_k#) peak at k = 9 with 1.3840127 < 1.38402.
  - **Repaired:** stated.
- **G3 (03 derivative condition).** 0.9593(1/ln L − 1/ln²L) − c < 0 holds for **all** L > 1, since max_u(1/u − 1/u²) = 1/4. The text's "ln L > 1.6" is unnecessary but harmless.
- **G4 (Pell count proof).** 02 asserts x′ > (2+√3)x for consecutive solutions of x² − du² = −1; 01 says the ratio is ≥ 3+2√2. An elementary proof, with no fundamental-unit theory:
  - Take two solutions α = x + u√d < β = x′ + u′√d. Then μ := β/α = −βᾱ lies in ℤ[√d], has norm 1 and μ > 1, so μ = s + t√d with s, t ≥ 1.
  - s = 2 forces d t² = 3, so d = 3, which is inadmissible because 3 | x²+1 is impossible. Hence s ≥ 3 and t√d = √(s²−1) ≥ √8.
  - Then x′ = s x + d t u > (s + t√d)x ≥ (3+2√2)x, using u√d = √(x²+1) > x.
  - With x₁ ≥ 1, the count of x ≤ N is ≤ 1 + ln N/ln(3+2√2) ≤ 1 + ln N/ln(2+√3). ✔
  - Lean has this as `pell_growth` (factor 3) plus `card_le_log_of_lacunary`.
- **G5 (02 C3 total).** The all-odd sub-case gives 2·max(2W100, cls), but the stated total lists only 2·cls. It is still valid: W100 ≤ W25′ gives 4W100 ≤ cls + 3W100 + W25′, so max(…, 2·max(2W100, cls)) = the stated max.
  (01 S3 lists 2·max(2W_o(1), cls) explicitly.) **Repaired:** note added.
- **G6 (FULL_PROOF table).** "1.3·10⁹ ≤ N ≤ 1.219·10¹⁹" should read **1.218524·10¹⁹** (the log's last endpoint). It still overlaps the tail at 10¹⁹. **Repaired** in FULL_PROOF.md.

## 5. Cited external results: what each is used for
| Result | Used where | Needed? |
| --- | --- | --- |
| Dusart 2010 (arXiv:1002.0442): π(x) ≤ (x/ln x)(1 + 1/ln x + 2.51/ln²x), x ≥ 355991 | bound3L/verify_exactL for P1 > 6·10⁷; tail_exact π(N^0.8) | Yes for the written chain. Hypothesis met: every use has x > 6·10⁷ |
| Robin 1983 (Acta Arith. 42): ω(n) ≤ 1.38402 ln n/lnln n, n ≥ 3 | tail_exact only | Yes for the written tail (with G2) |
| Sawhney (Problem_848.pdf), Sothanaphan (explicit threshold) | none (method origin only) | **No**: every lemma is proved in 02/01 |
| Hensel/CRT root counts, negative Pell | Lemma R, Lemma P | Elementary, proved (G4) |

**Removing both citations** (this is the Lean route, and it also works on paper):
- π(x) ≤ 48(⌊x/210⌋+1) + 4 (210-wheel; Lean `prime_count_le`).
- ρ(b) ≤ 4·2^{|T|}·2^{|W|} ≤ 2048·N^{1/4}, since 16^{|W|} ≤ 31^{|W|} ≤ N.
- P1 = R1 = ⌊N/1000⌋ (so the Pell D is O(1)).

With these, the tail is **linear** with no decreasing-β argument needed. `leanprep/tail_lean_model.py` (float model) gives worst B/N = 0.03314 for every N ≥ 10²⁴. The Lean cover already reaches 1.03·10²⁸.

## 6. Remaining caveat (unchanged)
The large computations are machine-checked but not human-refereed:
- the replay;
- exact rationals (141 intervals);
- the interval tail.

The Lean formalisation (LEAN_PLAN.md) is the route to removing this caveat.

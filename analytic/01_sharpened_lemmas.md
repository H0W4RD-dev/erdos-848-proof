# Threshold side 01: sharpened lemmas and the case analysis behind bound3.py

Notation follows `../02_lemmas.md`. A ⊆ [N] is valid, and A7, A18, A* are as there.
The goal is |A| < N/25 − 1 whenever A ⊄ A7 and A ⊄ A18. That suffices, because M(N) = ⌊(N+18)/25⌋ ≥ (N−6)/25.
A "class" means a residue class mod q ∈ {25, 50, 100} intersected with [N]. Its size is ≤ N/q + 1.

## Witness structure
For 1 ≤ b ≤ N, write b = 2^{v'} m with m odd. Let v = min(v', 3), c_0 = c_1 = 1, c_2 = 2, c_3 = 4.
Let T = {p ∈ {3,5,…,29} : p | b}, and let w be the number of prime factors of m that exceed 29.
Then w ≤ w_max(v, T), the largest w with 2^v·∏T·(31·37·…, w primes) ≤ N. Put base := 2^v ∏T ≤ b.
Write ρ_st := c_v 2^{|T| + w_max}.
**Lemma R (PROVED, standard).** For every unit u mod b, #{y mod b : y² ≡ u} ≤ c_v ∏_{p | m} (1 + (u/p)) ≤ ρ_st.
(Odd prime powers have 1 + (u/p) roots; modulo 2^{v'} there are ≤ c_v roots.)

## Lemma S2 (character-refined large-square tail, PROVED)
Let b ≤ N have structure st, let P1 ≥ 1, and let J ⊆ T, with L := ∏J and
G := {g ∈ [1, L] : g is a nonzero quadratic residue mod every p ∈ J}. Then
  #{a ≤ N : ∃ integer r > P1 with r² | ab+1}
  ≤ ρ_st · [ (|G|/L)(2(N+1)/P1 + (N²+1)/P1²) + |G| + √((N+1)/base) · Σ_{g∈G} g^{−1/2} ],
and also ≤ ρ_st (2(N+1)/P1 + (N²+1)/P1²), which is the J = ∅ form (Lemma W of 02).

*Proof.*
1. Put X := bN + 1. For each such a write ab + 1 = k r² with r > P1, so k ≤ K := X/P1². The map a ↦ (k, r) is injective.
2. From k r² ≡ 1 (mod b): gcd(k, b) = 1 and r² ≡ k⁻¹ (mod b).
3. So r lies in at most ρ_k := #{y mod b : y² ≡ k⁻¹} classes mod b. By Lemma R, ρ_k ≤ c_v ∏_{p|m}(1 + (k/p)).
   - For p ∈ J the factor is 0 unless k mod p is a nonzero QR, and then it is 2.
   - So ρ_k ≤ ρ_st · 1[k mod L ∈ G].
4. Since r ≤ √(X/k), the count is ≤ Σ_{k ≤ K, k mod L ∈ G} ρ_st (√(X/k)/b + 1).
5. Bound the two sums:
   - Σ_{k ≤ K, k ≡ g (L)} k^{−1/2} ≤ g^{−1/2} + (2/L)√K, by comparing the terms g + jL, j ≥ 1, with an integral.
   - #{k ≤ K : k mod L ∈ G} ≤ |G|(K/L + 1).
6. Simplify, using b ≥ base:
   - √X/b = √(N/b + 1/b²) ≤ √((N+1)/base);
   - (√X/b)·(2/L)√K = 2X/(bLP1) ≤ 2(N+1)/(LP1);
   - K ≤ (N²+1)/P1².

For J = ∅ the sharper Σ_{k≤K} k^{−1/2} ≤ 2√K gives the second form. ∎

*Effect.* A prime p ∈ J multiplies ρ by 2, but |G|/L by (p−1)/(2p), so its net cost is (1 − 1/p) instead of 2.
The constant |G| + √(N/base)·Σ g^{−1/2} limits how large J can be (the code uses J with L ≤ 3e5, minimised over J).

## Lemma S1 (witness group, PROVED)
Fix b ≤ N and a set 𝒞 of classes, where class C_i is mod q_i ∈ {25, 50, 100}. Assume that on every class, 4 ∤ ab+1 and 25 ∤ ab+1.
Let P0 = {3,7,11,13,17,19,23,29} \ T and k0 = |P0|. Then for any 29 ≤ Qs < P1:
  #{a ∈ ∪C_i : μ(ab+1) = 0}
  ≤ Σ_i [ (N/q_i + 1)(1 − ∏_{p∈P0}(1 − p⁻²)) + 2^{k0} − 1 ]   (primes ≤ 29, inclusion–exclusion per class)
  + Σ_i Σ_{29<p≤Qs} (N/(q_i p²) + 1)
  + Σ_{Qs<p≤P1} (N/p² + 1)
  + Tail_{S2}(b, P1).

*Proof.* Each such a has a prime p with p² | ab+1. Here p ≠ 2 and p ≠ 5 (by hypothesis), and p ∤ b (since ab+1 ≡ 1 mod p).
- p ≤ 29: inclusion–exclusion inside each class. Every term is one residue class mod a modulus coprime to q_i (it is odd and prime to 5), so its count is within 1 of (length)/(modulus), and there are 2^{k0} − 1 terms.
- 29 < p ≤ Qs: one class mod q_i p² per class.
- Qs < p ≤ P1: **{a ∈ ℤ : p² | ab+1} is ONE class mod p²**, so all classes together hold ≤ N/p² + 1 such a.
- p > P1: Lemma S2, which counts all a ≤ N.

The union bound over the ranges finishes the proof. ∎

The shared ranges (Qs, P1] and (P1, ∞) are paid **once per witness b**, however many classes b controls.
In 02_lemmas.md each class paid them separately (W100 was ~75% tail at 1e8).

## Lemma C1' (Case 1: A* has an even element b, PROVED)
All 8 classes mod 100 of A7 ∪ A18 (that is, 2 classes mod 25) have ab+1 odd and 25 ∤ ab+1.
For x ∈ A*, x needs:
- p² | x²+1 for some p ≡ 1 (mod 4), p ≥ 13 (𝒢 = such p ≤ R0, or p > R0), and
- q² | bx+1 for some odd prime q ∤ b.

Split by (p, q). The pieces of 02/J are unchanged except:
- **same** (q = p): needs p² | b²+1, and then x ≡ b (mod p²).
  - Class s_b (x ≡ −1/b mod 25) is already counted in full by the s_b piece.
  - So x lies in n_o := 22 classes mod 25 (23 if 5 | b), and the count is ≤ min(n_o(N/(25p²) + 1), N/p² + 1).
  - Maximise Σ_{p∈S} over prime sets S such that the CRT class of ±√−1 mod ∏_{p∈S} p² meets [1, N] (DFS as in 02).
- **q ∈ (13, Qs]**: 2(N/(p²q²) + 1) per pair, with q ∉ T.
- **q > Qs**: shared with A7 ∪ A18 through Lemma S1. Every a ≤ N, in A7, A18 or A*, with q² | ab+1 lies in one class mod q².

Total:
  |A| ≤ 2[(N/25 + 1)δ(P0) + 2^{k0} − 1] + 2Σ_{29<q≤Qs}(N/(25q²) + 1)
      + sb + oth + same + med(Qs) + E1(R0) + Σ_{Qs<q≤P1}(N/q² + 1) + Tail_{S2}(b, P1).

The terms are:
- sb = Σ_{p∈𝒢} 2(N/(25p²) + 1) if 5 ∤ b, else 0;
- oth = n_o Σ_{p∈𝒢} 2[(N/(25p²) + 1)δ(QS_p) + 2^{|QS_p|} − 1], with QS_p = {3,7,11,13} \ ({p} ∪ T);
- E1 = #{x ≤ N : ∃ p ≡ 1 (4), p > R0, p² | x²+1} ≤ 2N Σ_{R0<p≤R1, p≡1(4)} p⁻² + 2(π₁(R1) − π₁(R0)) + Pell(N, R1).

**Pell (PROVED).** The admissible d < D := (N²+1)/R1² number ≤ D/2 + 1 (3 ∤ d and 4 ∤ d; residues 1,2,5,7,10,11 mod 12).
Each d has ≤ 1 + log N / log(3+2√2) solutions with x ≤ N, because consecutive solutions grow by a factor ≥ v + w√d with v ≥ 3 (v = 2 forces d = 3).

## Case analysis S3 (PROVED)
Facts used:
- (i) 25 | xw+1 with x ≡ 7 (resp. 18) mod 25 iff w ≡ 7 (resp. 18). So A* elements are witnesses for both A7 and A18; A18 elements for A7; A7 elements for A18.
- (ii) 4 | xw+1 iff x, w are odd and x ≢ w (mod 4).

The classes mod 100 are:
- A7: 7 (≡ 3 mod 4), 32 (≡ 0), 57 (≡ 1), 82 (≡ 2);
- A18: 18 (≡ 2), 43 (≡ 3), 68 (≡ 0), 93 (≡ 1).

Below, W_o(c) and W_e(c) are Lemma S1 bounds for c classes mod 100, maximised over odd or even witness structures. W_e25 is one class mod 25 with an even witness. cls := N/100 + 1.

**C1.** A* has an even element: Lemma C1'.

**C2.** A* ≠ ∅ and A* is all odd.
- (i) A* has elements b ≡ 1 and b′ ≡ 3 (mod 4). Every x ∈ A7 ∪ A18 is controlled: even x by b, odd x by the witness of its own residue.
  - |A| ≤ Q_odd + W_o(6) + W_o(2), where Q_odd counts the 23 odd classes mod 50 by Lemma Q of 02.
- (ii) A* ⊆ {≡ u (mod 4)}. Then |A*| ≤ Q_u: 23 classes mod 100 (odd, ≡ u mod 4, ≢ ±7 mod 25).
  - b controls the 6 classes of A7 ∪ A18 that are even or ≡ u.
  - Let P7, P18 be the odd ≡ −u classes of A7 and A18. Elements of P7 and P18 control each other (both odd, same residue).
  - P7 ≠ ∅ ≠ P18: W_o(6) + 2W_o(1).
  - P7 ≠ ∅ = P18 with an even z ∈ A18: W_o(6) + W_e(1).
  - P7 ≠ ∅ = P18 with A18 ⊆ {odd ≡ u}: A18 lies in 1 class, so W_o(4) + cls.
  - Symmetric cases are the same; P7 = P18 = ∅ gives W_o(6).
  - |A| ≤ Q_u + max(W_o(6) + 2W_o(1), W_o(6) + W_e(1), W_o(4) + cls).

**C3.** A* = ∅, A7 ≠ ∅ ≠ A18.
- Both contain even elements: 2W_e25.
- Only A18 has an even z (or symmetrically):
  - A7 ⊆ {7, 57} via z: W_e(2).
  - Pick x0 ∈ A7 with residue u0. It controls A18's 2 even classes and its odd ≡ u0 class: W_o(3).
  - A18's odd ≡ −u0 class: W_o(1) via an A7 element ≡ −u0, if one exists; else cls.
- Neither contains even elements: for each u ∈ {1, 3}, either both A7^u and A18^u are nonempty (2W_o(1)) or one is empty (≤ cls).
- |A| ≤ max(2W_e25, W_e(2) + W_o(3) + max(W_o(1), cls), 2·max(2W_o(1), cls)).

The cases C1, C2, C3 are exhaustive for A ⊄ A7, A18. ∎

## Monotonicity and the interval check (PROVED)
Every bound is a max over witness structures (a set that grows with N) of a min over **fixed** grids (Qs ∈ QSGRID, P1 ∈ PGRID, J ⊆ T with ∏J ≤ 3e5) of terms that are nondecreasing in N.
So B3 is nondecreasing, and B3(N1)(1 + 1e−9) < N0/25 − 1 certifies [N0, N1].

**Float rounding.** All summands are positive:
- prime sums carry +1e−9 absolute slack (the true accumulated error is < 2e−10);
- π is exact (sieve to 6e7);
- the 1e−9 relative factor covers the remaining rounding.
An exact-rational re-evaluation is in `verify_exact.py` (see SUMMARY).

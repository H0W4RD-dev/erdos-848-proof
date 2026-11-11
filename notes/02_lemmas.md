# 02: Explicit lemmas and the refined case analysis (target: N ≥ N_A ≈ 1e9)

Throughout, N ≥ 1 and A ⊆ [N] is valid. A7, A18 and A* are as in Sawhney / Sothanaphan.
It suffices to show |A| < N/25 − 1 whenever A ⊄ A7 and A ⊄ A18, since M(N) ≥ (N−6)/25.
For a residue class S = {a ≤ N : a ≡ t (mod q)}, |S| ≤ N/q + 1.

## Lemma R (root count, PROVED, standard)
ρ(b) := max over units u mod b of #{y mod b : y² ≡ u}.
For b = 2^v m with m odd, ρ(b) ≤ c_v 2^{ω(m)}, where c_0 = c_1 = 1, c_2 = 2 and c_v = 4 for v ≥ 3.

## Lemma W (single witness, PROVED)
Setup:
- q ∈ {25, 100}, t mod q, and 1 ≤ b ≤ N.
- For every prime ℓ | q, ℓ² ∤ ab+1 for all a ≡ t (mod q).
- 𝒫0 = {primes p ≤ 29 : p ∤ qb}. Let P1 ≥ 29.

Then
#{a ∈ S : μ(ab+1) = 0} ≤ |S| [1 − ∏_{p∈𝒫0}(1−p⁻²) + Σ_{29<p≤P1} p⁻²] + (2^{|𝒫0|} − 1) + (π(P1) − π(29)) + ρ(b) [2(bN+1)/(b P1) + (bN+1)/P1²].

Proof.
- Every a with μ(ab+1) = 0 has p² | ab+1 for some prime p ∤ b. Primes p | b are impossible because ab+1 ≡ 1 (mod p).
- p ≤ 29: inclusion–exclusion over 𝒫0. Each condition (∏U)² | ab+1 is one residue class mod (∏U)², coprime to q.
  An AP of length |S| with difference coprime to the modulus m hits a class ⌊|S|/m⌋ or ⌈|S|/m⌉ times, so each term has error < 1.
- 29 < p ≤ P1: each p hits S at most |S|/p² + 1 times.
- p > P1: write ab+1 = kp² with 1 ≤ k ≤ K := (bN+1)/P1². Then kp² ≡ 1 (mod b), so p lies in one of ≤ ρ(b) classes mod b, and p ≤ √((bN+1)/k).
  The count is ≤ Σ_{k≤K} ρ(b)(√((bN+1)/k)/b + 1) ≤ ρ(b)(2√(K(bN+1))/b + K). Note √(K(bN+1)) = (bN+1)/P1.
  Each a is counted at least once. ∎

The bound is increasing in b, so for worst-case purposes take b ≤ N. Witnesses with many small prime factors lose those primes from 𝒫0, which shrinks the density term; the calculator maximises over the structure of b.

## Lemma P (Pell tail, PROVED; same argument as Sothanaphan Prop 2)
If x² + 1 = d u² with x, u ≥ 1, then d ≥ 2 is not a square, 4 ∤ d, and every odd prime factor of d is ≡ 1 (mod 4).
Consecutive positive solutions satisfy x' > (2+√3) x, and the smallest satisfies x₁² ≥ d − 1.
So #{x ≤ N : x² + 1 = d u²} ≤ 1 + log(N/√(d−1)) / log(2+√3).
Hence T(N, R1) := #{x ≤ N : ∃ prime p > R1, p² | x²+1} ≤ Σ over admissible d with 2 ≤ d < (N²+1)/R1² of that count.

## Lemma Q (outsider counts, PROVED)
Let 𝒢(R0) = {p ≡ 1 (mod 4) prime : 13 ≤ p ≤ R0}. Define
E1(N, R0) := #{x ≤ N : ∃ p ≡ 1 (4), p > R0, p² | x²+1} ≤ Σ_{R0<p≤R1, p≡1(4)} 2(N/p² + 1) + T(N, R1).
For a set 𝒮 of residue classes mod Q (Q ∈ {25, 50, 100}) avoiding ±7 mod 25:
#{x ≤ N : x mod Q ∈ 𝒮, μ(x²+1) = 0} ≤ Σ_{p∈𝒢(R0)} 2|𝒮| (N/(Q p²) + 1) + E1(N, R0).
(For x ≢ ±7 mod 25, 25 ∤ x²+1. Also 4 ∤ x²+1, and q² ∤ x²+1 for q ≡ 3 mod 4.)

## Lemma J (joint bound for Case 1, PROVED)
Setup: b even, b ≢ ±7 (mod 25). We bound J(b) := #{x ≤ N : x ≢ ±7 (25), μ(x²+1) = 0, μ(bx+1) = 0}.
- Since bx+1 is odd, the prime q with q² | bx+1 is odd, and q ∤ b.
- If 5 ∤ b, then 25 | bx+1 exactly on one class s_b mod 25 (s_b ≢ ±7).

J(b) ≤ sum of:
- [class s_b] Σ_{p∈𝒢} 2(N/(25p²) + 1);
- [other classes s, p ∈ 𝒢, root ε] Σ 2[(N/(25p²) + 1) δ_p + 2^{k_p} − 1], where δ_p = 1 − ∏_{q∈QS_p}(1−q⁻²), QS_p = {3, 7, 11, 13} \ ({p} ∪ {q | b}), and k_p = |QS_p|;
- [q = p] Σ over p ∈ 𝒢 with p² | b²+1 of (N/p² + 1). Here x ≡ −1/b ≡ b (mod p²) is a root.
  The worst case over b ≤ N is found by searching for root classes of ∏p² that contain an integer ≤ N.
- [13 < q ≤ Q1, q ≠ p] Σ_p Σ_q 2(N/(p²q²) + 1);
- [p > R0] E1(N, R0);
- [q > Q1] Σ_{Q1<q≤P1}(N/q² + 1) + ρ(b)(2(N+1)/P1 + (N²+1)/P1²).

## Refined case analysis (PROVED; all W below are Lemma W with the stated exclusions)
W25' means q = 25 with 2 and 5 excluded; W100 means q = 100 with 2 and 5 excluded; Q*odd is Lemma Q with the 23 odd classes mod 50 avoiding ±7.

- **C1.** A* has an even b. Then |A| ≤ W25'(b; 7) + W25'(b; 18) + J(b).
- **C2.** A* ≠ ∅ and all of A* is odd; fix b ∈ A*. Then |A*| ≤ Q*odd. For x ∈ A7 ∪ A18:
  - 4 | xb+1 ⟺ x is odd and x ≡ −b (mod 4).
  - Write P7 = A7 ∩ {odd, ≡ −b mod 4} and Q7 = A7 \ P7, and similarly P18, Q18. Each P is one class mod 100.
  - **C2a** (A7 ∪ A18 all odd):
    - Q7 and Q18 are each one class mod 100, with xb+1 ≡ 2 (mod 4), so each is ≤ W100(b).
    - If P7 and P18 are both nonempty, then for x ∈ P7 and z ∈ P18: xz+1 ≡ 2 (mod 4) and ≡ 2 (mod 25). So |P7| ≤ W100(z0) and |P18| ≤ W100(x0).
    - Total: |A| ≤ max(4W100, N/100 + 1 + 2W100) + Q*odd.
  - **C2b** (some even element, WLOG b' ∈ A7 even; the A18 case is symmetric):
    - |A18| ≤ W25'(b'; 18), or ≤ W100(b') if A18 lies in one class mod 100.
    - The A7 classes other than P7 (two even classes and one odd ≡ b mod 4) give ≤ 3W100(b).
    - P7 ≤ N/100 + 1; or, if A18 has an element z that is even or ≡ −b (mod 4), then P7 ≤ W100(z).
    - Total: |A| ≤ max(4W100 + W25', N/100 + 1 + 4W100) + Q*odd.
- **C3.** A* = ∅, A7 ≠ ∅ ≠ A18.
  - If some element is even (WLOG x0 ∈ A7), then |A18| ≤ W25'(x0). Let z0 ∈ A18: if z0 is even, |A7| ≤ W25'(z0); if odd, |A7| ≤ N/100 + 1 + 3W100.
  - If all elements are odd: split by mod-4 residue u. A7^u is controlled by A18^u (≤ W100) when A18^u ≠ ∅, else by its class (≤ N/100 + 1).
  - Total: |A| ≤ max(N/100 + 1 + 3W100 + W25', 2W25', 2(N/100 + 1)).

Asymptotic densities (W100 ≈ 0.00156, W25' ≈ 0.00623, Q*odd ≈ 0.0126):
- C1 ≈ 0.0125 + J (≈ 0.005 plus the root-class term ≈ 0.011);
- C2a ≈ 0.0257; C2b ≈ 0.0288; C3 ≈ 0.0209.
All of these are well below 0.04. Sothanaphan's corresponding margins were 0.0024 (C1) and 0.0043 (C2a).

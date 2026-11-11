# Threshold side 00: audit of Sothanaphan, "An Explicit Threshold in Erdős Problem #848" (Mar 24 2026, N0 = 2.64e17)

Source: `../src/g_1ujhm4_WYpgRV_rd1rJXIfHyvx16COEKe.txt`. This is an unrefereed, AI-generated paper.
I checked it line by line on 2026-11-10 and **found no mathematical error**. Every numerical constant reproduces (checks below).
The threshold is large only because of lossy error terms, not because of any flaw in the method.

## Line-by-line

| Item | Claim | Verdict | Where it is lossy |
| --- | --- | --- | --- |
| Lemma 1 | ρ(q) ≤ λ q^{1/4}, λ = 128/(8·15015)^{1/4} < 6.876 | PROVED (checked). The max of 2^k/(3·5·…·p_k)^{1/4} is at k = 5 (2.89 vs 2.85 at k = 6). 17^{1/4} > 2 is OK. λ = 6.8755 | Replaces the true ρ(b) ≤ c_v 2^{ω(b)} by a q^{1/4} power. For b ≤ 1e9 the gap is moderate (1024 vs ~2700), but the key cost is the error *shape* below |
| Prop 1 | \|S2(x;q,a) − (x/q)∏_{p∤q}(1−p⁻²)\| ≤ D + x/(qD) + λq^{1/4}(2x/(qD) + x/D²) | PROVED. Standard Möbius split at D; tail via n = md², d² ≡ a·m̄ | Two-sided estimate; only the upper side is needed |
| Cor 1 | #{a ≤ N : a ≡ t (q0), μ(ab+1) = 0} ≤ (N/q0)(1 − ∏_{(p,q0b)=1}(1−p⁻²)) + κ N^{3/4}, κ25 = 4.7, κ100 = 5.3 | PROVED. The reduction to d1 n + c1 is checked: g is squarefree and divides d1, so μ² is unchanged. Each displayed inequality is checked. E25(2.64e17) = 3.125 + 1.574 = 4.699; E100 = 3.5 + 1.775 = 5.275 (recomputed) | **Main loss.** The modulus d1 ≈ q0·b ≈ 25N, so the error is κN^{3/4} with κ → α + λq0^{1/4}/α² ≈ 4.7 (constant). At N = 1e9, κN^{−1/4} = 0.026: unusable below ~1e16 |
| Lemma 2 | C_quad < 0.02734669 | PROVED (recomputed): ∏_{13 ≤ p ≤ 1e7, p ≡ 1 (4)} (1 − 2/p²) = 0.97265351066, giving C_quad ≤ 0.02734668 | n/a |
| Prop 2 | #{x ≡ t (q) : μ(x²+1) = 0} ≤ X·C_quad + R(1 + log R) + 2X(log R + 2)/R + ((N²+1)/R²)(1 + log_{2+√3} N) | PROVED. A_r has error ≤ 2^{ω(r)}; Σ_{r≤R} 2^ω ≤ Σ d(r) ≤ R(1+log R). The Pell growth step x2 > (2+√3)x1 is checked | **Second loss.** A full two-sided Möbius inversion pays Σ_{r≤R} 2^{ω(r)}, i.e. R log R with R ≈ N^{2/3}, about 28·N^{2/3}. An upper bound only needs the union bound over primes (Σ_p 2(X/p² + 1)), which costs ~π(R)·2 instead. Base 2+√3 can be 3+2√2, since the only d with v = 2 is d = 3, which is inadmissible |
| Cor 2 | (a) \|A*\|/N ≤ (23/25)C_quad + 0.001973; (b) odd version (23/50)C_quad + 0.001973 | PROVED. The error 8.577e−5 per class at N = 2.64e17 is recomputed ≈ 5.66e−5 (R log R) + 2.9e−5 (Pell) | Uses only μ(x²+1) = 0 for A*. It ignores the joint constraint μ(bx+1) = 0 against the witness b, worth 0.025 → ~0.016 in density |
| (1), (2) | C_{2,5} = 0.1556568030, C_5 = 0.3667426022 | PROVED (recomputed) | n/a |
| Case 1 | 23/25·C_quad + 2/25·C_{2,5} = 0.0376115 < 0.04 − errors | PROVED; recomputed 0.0376115 | **Binding**: margin 0.00239 = 0.001973 + 0.000415 exactly at N0. This is the line that fixes 2.64e17 |
| Case 2a | 23/50 C_quad + 1/50 + C_{2,5}/50 = 0.0356927 | PROVED. The mod-4 class logic is checked | Two whole classes mod 100 are counted trivially |
| Case 2b | 23/50 C_quad + C_5/25 + C_{2,5}/25 = 0.0334755 | PROVED. Witness b (odd) for A7: 25 ∤ 7b+1 since b ≢ 7. Witness b′ ≡ 7 even for A18 | Odd witness gives C_5, which loses the factor at 2 |
| Case 3 | 2C_5/25 = 0.0293395 | PROVED | n/a |
| Final step | \|A\| < N/25 ⟹ \|A\| ≤ ⌊(N+18)/25⌋ | PROVED. ⌈N/25⌉ − 1 ≤ (N−1)/25 | n/a |

## Where N0 = 2.64e17 comes from
Case 1 has asymptotic density 0.0376 against 0.04, so the margin is only 0.0024. Two error terms must fit inside it:
- Cor 1 (κN^{−1/4}): needs N ≳ 1e17 on its own (9.4·N^{−1/4} < 0.0004 ⟺ N > 3e17).
- Cor 2: ≈ 23·28·N^{−1/3}·(1.3).

No single fix to Sothanaphan's own inputs can get below ~1e12, because the κN^{3/4} shape is intrinsic to Prop 1 with modulus ≈ N.
The route that works is the one in `../02_lemmas.md`. It sieves small primes exactly, uses union bounds for the middle primes, uses the k-parametrised ρ(b)-tail (error ≈ (ρN²)^{1/3}), and uses the joint constraint for A*.
My sharpenings on top of that are in `01_sharpened_lemmas.md`.

## Flags (not errors)
- F-a. Cor 1 is stated for b ≤ N and N ≥ 2q0. Both hold in every application.
- F-b. In Case 2b the paper writes "4.7N^{−1/4}" for each of the two Cor 1 applications (q0 = 25). That is consistent.
- F-c. In Case 2a, 10.6 = 2·5.3 (q0 = 100). That is consistent.
- F-d. Lemma 2's tail Σ_{p>1e7} 2/p² < 2/(1e7 − 1) is fine (crude but valid).

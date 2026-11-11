# 04: Independent audit of threshold/01_sharpened_lemmas.md and threshold/bound3.py (session c7f777, 2026-11-10)

Method: re-derived each lemma from scratch, brute-force checks with my own code (audit_S12.py, audit_roots.py; no imports from threshold/), and a line-by-line code review of bound3.py.

| Item | Verdict | Notes |
|---|---|---|
| Lemma R (root count with Legendre factors) | **CORRECT** | Exact test, audit_roots.py: 48,040 (b, k) pairs, b < 4000, 0 violations of #{y : y² ≡ k⁻¹ (b)} ≤ c_v ∏_{p\|m}(1 + (k/p)). |
| **S2** (character-refined tail) | **CORRECT** | Re-derived.<br>- a ↦ (k, r) is injective; gcd(k, b) = 1; ρ_k ≤ ρ_st·1[k ∈ G].<br>- Σ_{k≡g (L), k≤K} k^{−1/2} ≤ g^{−1/2} + (2/L)√K, by an integral comparison over j ≥ 1.<br>- #{k ≤ K : k mod L ∈ G} ≤ \|G\|(K/L + 1).<br>- √X/b ≤ √((N+1)/base), since base \| b; X/(bP1) ≤ (N+1)/P1; K ≤ (N²+1)/P1².<br>Brute force at N = 3e4 (25 witnesses, J = ∅ and J ⊆ T): all bounds hold. |
| **S1** (witness group, shared ranges) | **CORRECT** | The key step is that {a : p² \| ab+1} is a single class mod p² when p ∤ b, so ranges > Qs are counted once over all a ≤ N.<br>Inclusion–exclusion moduli are odd and prime to 5, hence coprime to q ∈ {25, 50, 100}. Per-class counts N/(qp²)+1 are correct.<br>Brute force at N = 3e4: holds. |
| C1′ (Case 1 with shared ranges, same-term restricted to 22/23 classes) | **CORRECT** | Every route for x ∈ A* is accounted for: q = 5 → s_b; q ∈ {3,7,11,13} → I-E; 13 < q ≤ Qs → pairs; q = p → same; q > Qs → shared; p > R0 → E1.<br>min(n_o(N/(25p²)+1), N/p²+1) is valid. |
| S5 Pell (≤ D/2+1 admissible d; ratio ≥ 3+2√2) | **CORRECT** | Residues mod 12 coprime to 3 and ≢ 0 mod 4 are {1,2,5,7,10,11}; checking every partial block gives count(1..n) ≤ n/2+1.<br>v = 2 forces d = 3, which is inadmissible because 3 ≡ 3 mod 4. |
| **S3** (case split C2(i), C2(ii), C3) | **CORRECT** | Facts (i)/(ii) re-proved: 25 \| xw+1 ⟺ w is in the same 7/18 class; 4 \| xw+1 ⟺ both odd with different residues mod 4. The residues mod 4 of the 8 classes mod 100 are right.<br>Each subcase re-derived, including W_o(4) + cls when P18 = ∅ and A18 has no even element (then b only needs to control 4 classes). Q_u counts 23 classes mod 100. Exhaustive given A ⊄ A7, A18. |
| bound3.py vs lemmas | **faithful** | WG, tailvec, Qcount/E1/pell_T, C1, B3 checked term by term. Grids are fixed (PGRID, QSGRID, J with ∏J ≤ 3e5), so B3 is nondecreasing.<br>Minor robustness note: `structures` would IndexError if w exceeded 70; this is impossible for N ≤ 1e100.<br>Float slack argument is plausible. The exact re-check in verify_exact.py is theirs; I did not audit verify_exact.py itself. |

**Overall: S1, S2, S3: CORRECT** (no fixes needed).

Caveat: the brute-force checks run at small N, where all bounds are very loose (actual/bound ≤ 0.02). They catch sign/structure errors, not constant-level slips. The constant-level check is the line-by-line comparison above.

## Addendum: audit of the exact verifiers and bound3L (requested by parent)

| Item | Verdict | Notes |
|---|---|---|
| bound3L.py vs bound3.py | **CORRECT** | Diff only:<br>- Dusart π(x) beyond 6e7 (valid for x ≥ 355991; log lowered by 1e-12, ceil + 1);<br>- π restricted to p ≡ 1 (mod 4), beyond 6e7: bounded by π(6e7) for those primes + (π(x) − π(6e7)), an upper bound;<br>- prime sums beyond 6e7 get + 1/(2(PLIM−2)) ≥ Σ_{odd n>PLIM} n⁻² ≥ Σ_{p>PLIM} p⁻²;<br>- the P1 grid extended to 1.4e16 (still fixed);<br>- full-G same-term for N > 1e12, an upper bound.<br>s2(a, ·) is only called with a ≤ PLIM; otherwise it would IndexError, not return a wrong value. |
| verify_exact.py | **CORRECT** | - **Exactness:** every summand is a Fraction or an int. Prime sums use ceil(10¹⁸/p²)/10¹⁸ (rounded up). √ is upper-rounded via isqrt + 1 and g^{−1/2} via S/isqrt(gS²). log N is float × (1 + 1e-12), above the true value. log(3+2√2) ≥ 1.762747174039 (true value 1.76274717403908…). π comes from the sieve; densities are exact products.<br>- **Parameters:** chosen by float search; legitimate, since any parameter choice gives a valid bound.<br>- **Expressions:** WGx, C1x, Qx, E1x, same_x and B3x match bound3/S1–S3 term by term. med_x drops the q = p pair's +2, which is valid because that pair is in the same-term.<br>- **Tiling:** endpoints are the log values rounded to integers and chained (pts[i+1] is the next interval's left end), so there are no gaps. The check is B3x(b_i) < a_i/25 − 1 with B3 nondecreasing. |
| verify_exactL.py | **CORRECT** | Same code; imports bound3L; exact tail constant; full-G same-term for N > 1e12. |
| tail_exact.py | **CORRECT** | mpmath iv with outward rounding (dps 40). Prime sums are exact rationals with the PLIM tail. Dusart is evaluated in intervals, and the upper endpoints are compared with the lower endpoint of N/25 − 1.<br>It evaluates exactly the tail_bound.py formulas (P1 = R1 = N^0.8; +15 ≥ 2^{\|QS_p\|} − 1; med without the q = p subtraction, an upper bound).<br>It relies on the 03 monotonicity argument (β(N)/N decreasing); I re-checked the derivative computation. |

## Independent spot reproduction (spot_exact.py: my own exact implementation of S1–S3; shares no code with threshold/)
| Interval | My exact worst B(b)/a | Verdict |
| --- | --- | --- |
| [16000000, 16160000] (the first interval) | 0.03974 | OK |
| [547326600, 729738000] (their tightest, 0.03997) | 0.04002 | My coarser grids (P1 ∈ 1.2^k, \|J\| ≤ 2) fail by 2e-5 |
| … re-covered as 4 sub-intervals (ratio 1.08) | 0.0327, 0.0326, 0.0325, 0.0318 | OK: the range is independently covered |
| [182699300000, 251211500000] (from the bound3L cover) | 0.03861 | OK |

Conclusion: the S1–S3 bound independently reproduces their certification. On the tightest stretch it needed a finer interval split, not a different bound.

# Erdős #848, threshold side: SUMMARY (2026-11-10, session c7f777-thr)

**RESUMED and FINISHED (computations) 2026-11-10 as an independent check of PR #4667. See PROOF_CHAIN.md.** Steps 1–3 of the parent's plan are done; the lemma audit and the certificate replay are pending (jsp-elementary-nt).

(Earlier: paused by the parent session; all jobs stopped and the sub-claim released.)
Reason: a prior claim exists. JSP PR TheJustinSunPrize/awards#4667 (q1lin570, open since 2026-10-02) claims a complete all-N Lean proof with the same architecture.
So this work can only serve as an independent check of #4667; whether to finish it is the user's decision.

## Result
- **Best analytic threshold: N0 = 1.6e7.** COMPUTER (float cover) + PROVED lemmas (self-audited); exact re-check partial (see next item).
  - Sothanaphan's threshold was 2.64e17.
  - The other session's `analytic2.py` threshold was ≈ 6e8 (single points) and 1e9 (cover).
- `bound3.py cover 1.6e7 1.3e9` covers **[1.6e7, 1.388e9]** in 40 monotone intervals. Log: `cover3_16e6_13e8.log`.
- `verify_exact.py` re-evaluates every interval bound in exact rationals (rational upper bounds for √, log and the prime sums). Log: `verify_exact_16e6_13e8.log`.
  - **15/40 intervals were verified exactly before the stop, covering [1.6e7, 2.63e7]. No failure.** The remaining 25 intervals passed in float only.
- Worst point: B3(N1)/N0 ≈ 0.0399 near N = 1.7e7; elsewhere ≤ 0.0398.

## Where the gain came from (see 01_sharpened_lemmas.md)
| Change | Effect |
| --- | --- |
| S1: share the q > Qs ranges and the ρ-tail once per witness, across all classes it controls, including A* in Case 1 | Removes 2–4 duplicate N^{2/3} tails per case |
| S2: character-refined tail. The k with r² ≡ 1/k (mod b) solvable must be QRs mod b's small primes | Each small prime p of b costs (1 − 1/p) instead of 2, so worst-case ρ effectively drops from 1024 to ~128 |
| S3: new mod-4 case split for C2/C3 | At most 3 witnesses per case; Q_u (density 0.0063) when A* lies in one residue mod 4 |
| S4: 'same' term restricted to 22 (or 23) classes mod 25 | −12% on the 0.011 root-class term |
| S5: Pell growth ratio 3+2√2 (d = 3 is inadmissible) | Small |

## Audit of Sothanaphan (00_audit_sothanaphan.md)
- No error found; every constant reproduces.
- N0 = 2.64e17 is forced by Case 1's margin of 0.0024 against two lossy errors:
  - the κN^{3/4} of Cor 1 (modulus ≈ 25N, so it is intrinsic to the method);
  - the R log R from the full Möbius inversion in Cor 2.

## How the ends meet (combined chain)
1. N ≤ 7.1e8: the colouring certificate of jsp-elementary-nt (stopped at 7.1e8).
2. [1.6e7, 1.388e9]: this side's B3 cover (Lemmas S1–S3 + 02 base lemmas).
3. [1e9, 7.96e13]: analytic2 cover (stopped at 7.957575e13).
4. N ≥ 1e19: `../tail_bound.py` closed form + monotonicity (`../03_analytic_range.md`). I re-ran it: C1 = 0.0315 at 1e19, 1e21 and 1e30. I checked the β(N)/N decreasing argument.

**Remaining gap: [7.96e13, 1e19] is uncovered by any interval run.** Every point tested there passes at ≈ 0.030.
Closing it is about 1 h of `analytic2.py cover 7.9e13 1e19 0.2`. bound3 would also do it, but its sieve cap PLIM = 6e7 limits P1, so it would need extending.

Other open items:
- the exact re-check of the 25 remaining intervals of [2.63e7, 1.388e9];
- an independent audit of S1–S3 (the subagent was stopped before reporting; only the brute-force lemma_test.py at N = 6e4 passed);
- a third-party check of the 02 base lemmas (I audited them and found no error).

Status: **NOT a complete proof.** OPEN pieces as listed above, and a prior claim exists in any case (PR #4667).

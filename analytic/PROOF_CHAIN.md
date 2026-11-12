# Erdős #848 (JSP-000701): proof chain for the independent check of PR #4667

**Claim.** For every N ≥ 1 and every A ⊆ [N] with μ(ab+1) = 0 for all a, b ∈ A (a = b allowed), |A| ≤ ⌊(N+18)/25⌋.

**Reduction (PROVED, 02 and 01 here).**
- If A ⊆ A7 or A ⊆ A18, the claim is immediate.
- Otherwise it suffices to show |A| < N/25 − 1, since M(N) ≥ (N−6)/25.
- Every analytic step below proves this for all A ⊄ A7, A18, on the stated range of N.

Status labels: PROVED (lemma with proof), COMPUTER (exact or interval-certified computation), PENDING (not yet done or not yet audited).

| N range | Covered by | Lemmas | Script / log | Arithmetic | Status |
| --- | --- | --- | --- | --- | --- |
| 1 ≤ N ≤ 7.1e8 | Proper M(N)-colouring of the compatibility graph, built by insertion (Lemmas C1, C2 of ../notes/01) | ../notes/01_colouring_certificate.md | ../certificate/colour4.py; ../certificate/cert_1e9.log ("OK through N=710000000") | Integer (exact squarefree test, cube-root trial division + square test) | COMPUTER. Independent replay with sympy: done to N ≤ 1e6 (replay_verify.py, superseded and not included). Independent sympy-only replay to N = 2e7 PASSED (../certificate/replay_2e7.log, 2,103,229 insertions) |
| 1.6e7 ≤ N ≤ 1.387721e9 | Monotone bound B3(N1) < N0/25 − 1 on 40 intervals | 02 base lemmas (W, Q, J, P) + S1, S2, C1′, S3 (01_sharpened_lemmas.md) | bound3.py cover → cover3_16e6_13e8.log; verify_exact.py → verify_exact_part1.log + verify_exact_part2.log ("ALL INTERVALS VERIFIED EXACTLY") | Exact rationals; rational upper bounds for √, log and the prime sums | COMPUTER done; S1–S3 audited CORRECT (../04_audit_threshold_S123.md) |
| 1.3e9 ≤ N ≤ 1.2e19 | Same bound, extended for large N (Dusart π, prime-sum tail beyond 6e7, full-G 'same' term) | as above | bound3L.py cover → coverL_13e8_12e19.log; verify_exactL.py → verify_exactL.log | Exact rationals | COMPUTER done: 101/101 intervals, "ALL INTERVALS VERIFIED EXACTLY" (worst 0.039961); audited CORRECT |
| N ≥ 1e19 | Closed-form bound αN + β(N), with β(N)/N decreasing on [1e19, ∞) | 02 base lemmas; ../notes/03_analytic_range.md (monotonicity); Robin 1983 ω-bound; Dusart 2010 | tail_bound.py (float); tail_exact.py (mpmath interval arithmetic + exact prime sums): "CERTIFIED OK", C1 = 0.031496 at 1e19 | Interval-certified | COMPUTER done; monotonicity argument checked by me (threshold side) |

Every exact log ends with "ALL INTERVALS VERIFIED EXACTLY" and contains no FAIL line: 40/40 intervals for the first cover and 101/101 for the second.

The ranges overlap at 1.6e7 ≤ 7.1e8, 1.3e9 ≤ 1.387e9 and 1e19 ≤ 1.2e19, so there is **no gap**. The certificate is needed only up to 1.6e7, and that range is independently replayed.

## Cited inputs
- Dusart (2010): π(x) ≤ (x/ln x)(1 + 1/ln x + 2.51/ln²x) for x ≥ 355991.
- Robin (1983): ω(n) ≤ 1.38402 ln n / ln ln n for n ≥ 3. Used only for N ≥ 1e19.
- Elementary facts: roots of y² ≡ u mod prime powers; growth of solutions to the negative Pell equation (proved in 02/01).

## Not float-only anywhere
The float scripts (bound3*.py, analytic2.py, tail_bound.py) only *search* parameters and intervals.
Every interval and the tail are re-evaluated by verify_exact*.py or tail_exact.py, whose arithmetic is exact rationals or outward-rounded intervals.
The other session's float-only cover [1e9, 7.96e13] (../cover_1e9_1e21.log) is **not used** in this chain; it is superseded by the bound3L exact cover.

## Status (2026-11-10)
- Both audits are done by jsp-elementary-nt: the lemma and code audit is in ../04_audit_threshold_S123.md, and the certificate replay is in ../replay_2e7.log.
- **Complete independent all-N proof chain:** PROVED lemmas + COMPUTER (exact or interval-certified), independently audited.
- Combined write-up: ../FULL_PROOF.md.
- Nothing has been posted. The prior claim is PR TheJustinSunPrize/awards#4667; any use of this work is the user's decision.

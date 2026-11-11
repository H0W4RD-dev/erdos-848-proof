# Erdős Problem #848: a complete, reproducible proof for every N

**Theorem.** For every N ≥ 1 and every A ⊆ {1,…,N} such that ab+1 is not squarefree for all a, b ∈ A (a = b allowed),

  |A| ≤ ⌊(N+18)/25⌋ = #{n ≤ N : n ≡ 7 (mod 25)}.

Equality holds for A = {n ≤ N : n ≡ 7 (mod 25)}, since 25 | ab+1 for all a, b ≡ 7 (mod 25). So the residue class 7 mod 25 is extremal for **every** N, not only for sufficiently large N. This answers [Erdős Problem #848](https://www.erdosproblems.com/848) (Justin Sun Prize catalogue entry JSP-000701) in full.

Sawhney proved the statement for all sufficiently large N, with an ineffective threshold. Sothanaphan made the threshold explicit (N ≥ 2.64·10¹⁷, unrefereed note, Mar 2026). This repository closes the remaining finite range and lowers the analytic threshold. Every step is either a written proof or an exact / interval-certified computation that you can rerun.

## How the proof is split

| Range of N | Method | Certified by |
|---|---|---|
| 1 ≤ N ≤ 2·10⁷ | Explicit proper colouring of the "ab+1 not squarefree" graph with ⌊(N+18)/25⌋ colours, maintained by augmenting-path repairs | `certificate/replay_verify2.py` replays the full operation log `certificate/ops_20000000.bin` (2,103,229 insertions). It uses sympy only and shares no code with the generator. Log: `certificate/replay_2e7.log` |
| 1.6·10⁷ ≤ N ≤ 1.388·10⁹ | Explicit counting bound (lemmas S1–S3, below 1/25 of N) | `analytic/verify_exact.py`: 40 contiguous intervals in exact rational arithmetic. Logs: `analytic/verify_exact_part1.log`, `analytic/verify_exact_part2.log` |
| 1.3·10⁹ ≤ N ≤ 1.219·10¹⁹ | Sharper explicit bound (bound3L) | `analytic/verify_exactL.py`: 101/101 intervals, exact rationals (worst case 0.039961 < 0.04). Log: `analytic/verify_exactL.log` |
| N ≥ 10¹⁹ | Closed-form tail bound | `analytic/tail_exact.py`: mpmath interval arithmetic, "CERTIFIED OK" (C₁ = 0.031496) |

The ranges overlap, so every N is covered. The full argument is in **[FULL_PROOF.md](FULL_PROOF.md)**. Lemma proofs are in `notes/02_lemmas.md` and `analytic/01_sharpened_lemmas.md`, and the chain of ranges is in `analytic/PROOF_CHAIN.md`.

## Reproducing it

Requires Python 3.13 with `numpy`, `numba`, `sympy` and `mpmath`.

```bash
# 1. Small-N certificate: independent replay of the stored operation log (~30 min)
cd certificate && python replay_verify2.py 2e7          # expect "REPLAY VERIFIED N=20000000"

# 2. Analytic ranges in exact arithmetic
cd analytic && python verify_exact.py  <cover log> <out>   # see PROOF_CHAIN.md for the exact commands per range
cd analytic && python verify_exactL.py <cover log> <out>
cd analytic && python tail_exact.py                     # expect "CERTIFIED OK"

# 3. Optional: a second, independent small-N certificate (N ≤ 1e5 in about 2 minutes)
cd small_n_check && python build_cover.py 1e5 log.txt && python verify_cover.py log.txt 1e5
```

To regenerate the operation log from scratch rather than replaying the stored one, run `python certificate/colour4_log.py 2e7 1 replay`.

## Checks performed

- The certificate replay verifier was mutation-tested: deliberately corrupted logs are rejected.
- A second certificate method (`small_n_check/`, a different builder and a sympy-only verifier) independently confirms the bound for every N ≤ 99,981.
- Lemmas S1, S2, S3, C1′ and S5 and the exact verifiers were each audited by a session that did not write them: `notes/04_audit_threshold_S123.md`.
- Sothanaphan's explicit-threshold note was audited, and no error was found: `analytic/00_audit_sothanaphan.md`.

## Relation to JSP PR #4667

[TheJustinSunPrize/awards#4667](https://github.com/TheJustinSunPrize/awards/pull/4667) (opened 2026-10-02) claims an all-N Lean proof with a similar overall structure. We compared the two (`notes/compare_PR4667.md`) and found no error in its written mathematics. However, at commit `2cfd228` its Lean modules read 2,788 certificate JSON files that are not in the repository, its manifest or any generator script, so its build cannot currently be reproduced (`notes/audit_PR4667.md`). This repository was developed independently and does not use that code.

## Credits and disclosure

- Problem: P. Erdős (see erdosproblems.com/848 for the original source).
- Asymptotic solution: M. Sawhney. Explicit threshold: N. Sothanaphan (unrefereed note, 2026). Both are cited, not reproduced.
- This proof: **[SOLVER NAME]**. It was produced with substantial assistance from Anthropic's Claude (several Claude Code sessions, 2026-11-10), which generated the lemmas, code and certificates and cross-audited each other's work. It has not yet been refereed by a human mathematician.

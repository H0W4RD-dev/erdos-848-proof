# Audit of JSP PR TheJustinSunPrize/awards#4667: Lean side (session c7f777, 2026-11-10)

Repo: github.com/q1lin570/JSP-lean @ 2cfd2281a215c5b662186e5fe5bd99eeedb857d2 (5.3 GB, 804 files).
Method: read-only. All 587 files under 2 MB were downloaded and grepped. The 217 generated `ReplayFullBlock*/Factors.lean` files (~24 MB each) and the 100 MB `.cert` file were not downloaded.
**Nothing was built or executed**; no lake or lean was run on their code.

## Verdict: **CAN'T TELL WITHOUT BUILDING, and the published repository cannot be built as it stands**
The design is sound in principle. But the proof depends on certificate data that is not in the repository, so no third party can currently re-check the claim.

## (1) Forbidden constructs (grep over all 555 downloaded `.lean` files)
- `sorry`, `admit`, `axiom` declarations, `native_decide`, `unsafe`, `implemented_by`, `@[extern]`, `ofReduceBool`, `opaque`, `#eval`: **none** in the proof sources.
  - The only `sorry` / `native_decide` hits are in a comment of `references/erdos848_asymptotic.lean`, which is a third-party file not imported.
- `run_cmd` appears only in `verification/Solution*Audit.lean`.
  - These collect the transitive axioms of `JSP848.upper_bound` and `JSP848.solution` from the **kernel-checked** environment (`env.checked`) and reject anything outside {propext, Classical.choice, Quot.sound}. That is a reasonable audit.
- **6 custom `elab` term macros** (`eventChunk%`, `eventSlice%`, `diagonalChunk%`, `diagonalSlice%` in CertificateInput.lean; `primeHintChunk%` in PrimeCountInput.lean; `replaySnapshot%` in SnapshotInput.lean).
  - At elaboration time they read JSON from disk (`IO.FS.readFile`) and emit **only data constructors** (`mkNatLit`, list literals of `BucketMove.mk` and similar).
  - So no proof enters through them: every property of the data must still be proved and kernel-checked (`theorem valid : RawReplay.Valid snapshot := …`). This is a legitimate technique, not a soundness hole.

## (2) Statement fidelity
`JSP848.UpperBound := ∀ N A, A ⊆ Finset.Icc 1 N ∧ (∀ a ∈ A, ∀ b ∈ A, ¬Squarefree (a*b+1)) → A.card ≤ (N+18)/25`.
- The diagonal case is included, and `solution` also proves the extremal class 7 mod 25. **Faithful to Erdős #848.**
- Versus formal-conjectures `Erdos848.erdos_848 : ∀ N, Erdos848For N`, which uses `A ⊆ Finset.range N` and the bound `#{n < N : n ≡ 7 (mod 25)}`:
  - 0 can never lie in A (0·0+1 = 1 is squarefree), so FC at N is exactly UpperBound at N−1.
  - **Mathematically equivalent**, but the PR does not prove FC's `erdos_848` literally; a short bridging lemma would be needed.

## (3) Are the finite and interval parts proved in Lean?
- Top level: `upper_bound` = `upper_bound_finite` (N ≤ 2^25) or `upper_bound_from_finite_cutoff` → `upper_bound_middle` (2^25..2^41) / `upper_bound_large`.
  - These are theorems, with no extra hypotheses or assumption parameters.
- **Gap: missing data.**
  - 330 modules in their own declared proof closure (`verification/proof_modules.list`; 330/330 are listed) read **2,788 distinct JSON files** from three directories:
    - `certificates/33554432-v1/` (diagonal and event chunks);
    - `certificates/prime-count33554432-v1/`;
    - `certificates/prime-one-count2097152-v1/`.
  - **None of these files is in the repository, and none is in `proof_source.sha256`.** No included script (check_*.py, *.cpp, verify_release.sh) generates them.
  - The README says raw replay snapshots are "not source dependencies of the selected theorem target and are not included". That is **contradicted by the source**: e.g. `ReplayFullBlock000/Replay.lean` does `replaySnapshot% "../../certificates/replay-snapshots33554432-v1/snapshot-0.json"`, and it is on the proof path.
  - Import-closure analysis (my script over their imports):

    | Module | Modules in closure | Data-reading |
    | --- | --- | --- |
    | FiniteRangeUpperBound | 557 | 249 |
    | MiddleRangeUpperBound | 217 | 81 |
    | LargeRangeUpperBound | 146 | 81 |
    | Solution | 766 | 330 |

    So **all three ranges**, not just the finite one, depend on the missing data.
  - Consequence: `lake build JSP848.Solution` on the published commit should fail at the first `…%` elaborator ("file not found").
  - Their claim that it "passed local kernel checking" can therefore not be reproduced from what was published.
  - Because the data only feeds kernel-checked statements, a successful build with *some* data would still be sound. The issue is **reproducibility**, not a demonstrated error.
- Build cost: the README notes the axiom audit needed a Lean memory limit of 262,144 MiB (256 GB).
- The 217 committed `Factors.lean` files (~5.3 GB) were not inspected.

## (4) PR thread
- OPEN, no labels, **no comments, no reviews**, last update 2026-10-02 (its creation).
- No maintainer engagement in over five weeks.
- The material is self-described as produced by "Codex GPT-6-Astra"; no independent review is claimed.

## Bottom line
- **Soundness: no defect found** in what is visible. No sorry, axiom or native_decide; the data-only elaborators are sound; the statement is faithful.
- **Verifiability: gap.**
  - The certificates that the Lean proof reads are absent and unpinned, so the formalisation is currently not reproducible by anyone else.
  - A complete check needs (a) the missing `certificates/` data from the author, and (b) a machine with very large memory.
- **Does this make our proof "first valid"?** Not demonstrably.
  - The PR's claim is unverified, not refuted. Its written-proof architecture (certificate to 2^25 plus explicit ranges) matches ours, and we found no mathematical error.
  - Priority is the maintainers' call. A maintainer would reasonably ask the author for the missing data rather than reject the claim.

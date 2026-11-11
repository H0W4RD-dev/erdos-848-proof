# 03: Covering every N ≥ 1e9 analytically

All lemmas are in 02_lemmas.md. Code: analytic2.py (intervals) and tail_bound.py (N ≥ 1e19).

## Monotone bound B(N) (PROVED, given the lemmas)
analytic2.B(N) returns, for each case C1, C2a, C2b, C3, a sum of the Lemma W/Q/J/P bounds.
Every free parameter (P1, R1) is minimised over a **fixed** grid GRID = {round(1.04^k)} ∩ [30, 1e22]; R0 = Q1 = 1000.

Every summand is nondecreasing in N:
- linear terms N/q, N/p²;
- the ρ-tails;
- the Dusart-bounded π(P1);
- the Pell count (increasing in D = (N²+1)/R1² and in N);
- the root-class term (more prime sets become feasible as N grows);
- the set of witness structures (it grows with N).

A minimum over a fixed grid of nondecreasing functions is nondecreasing, and so is a maximum over a growing set.
Hence |A| ≤ B(N1) for every N ∈ [N0, N1] and every valid A not contained in A7 or A18.
The interval check is B(N1)(1+1e-9) < N0/25 − 1. That gives |A| < N/25 − 1 < M(N) on the whole interval.

## Explicit inputs (CITED)
- π(x) ≤ (x/ln x)(1 + 1/ln x + 2.51/ln²x) for x ≥ 355991 (Dusart 2010). Below 2e6 the sieve gives π exactly.
- Σ_p p⁻² = P(2) = 0.452247420041065498506543…
- ω(n) ≤ 1.38402 ln n / ln ln n for n ≥ 3 (Robin 1983). Only tail_bound.py uses this.

## Tail N ≥ 1e19 (PROVED modulo the numerical evaluation)
Take P1 = R1 = N^0.8, and bound ρ(b) ≤ 4·2^{ω(b)} ≤ ρ*(N) := 4 N^{1.38402 ln 2 / ln ln N} for every b ≤ N.
Then use the worst density over structures (T = ∅) and the worst ρ. Each case bound has the form αN + β(N), where β(N)/N is a finite sum of:
- c/N;
- π(N^0.8)/N, which with Dusart is ≈ N^{−0.2}/(0.8 ln N): decreasing;
- ρ*(N)(2N^{−0.8}(1+1/N) + N^{−0.6}(1+N^{−2})).
  Here ρ*(N) N^{−0.6} = 4 exp(L(0.9593/ln L − 0.6)) with L = ln N.
  The derivative in L is 0.9593/ln L − 0.9593/ln²L − 0.6 < 0 for ln L > 1.6, so this is decreasing; the same holds with −0.8;
- (D/2 + 2)(1 + ln N / ln(2+√3))/N with D = N^{0.4} + N^{−1.6}: decreasing for N > e^{2/0.6}.

So β(N)/N is decreasing on [1e19, ∞), and it suffices to evaluate at N* = 1e19.
tail_bound.py gives C1 = 0.03150, C2a = 0.02593, C2b = 0.02908, C3 = 0.02100, all < 0.04 − 1e−19. It is also OK at 1e21 and 1e40.

## Intervals [1e9, 1e21] (COMPUTER, analytic2.py cover)
- Log: cover_1e9_1e21.log. Adaptive ε, starting at 3%.
- Tightest ratio B(N1)/N0 ≈ 0.0398 near N0 ≈ 1.2e9 (margin ≈ 2e5 elements).

## Point values B(N)/N (worst case)
| N | 5e8 | 1e9 | 2e9 | 1e10 | 1e12 | 1e18 |
| --- | --- | --- | --- | --- | --- | --- |
| B(N)/N | 0.0407 (FAIL) | 0.0381 | 0.0364 | 0.0332 | 0.0299 | 0.0301 |

The binding case is C2b near 1e9 and C1 (root-class term) for large N.

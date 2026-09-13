# Elkies family search

The Explorer search uses only the Elkies family. Polynomial data are transcribed
from [Elkies, arXiv:2608.25406v1, Theorem 4](https://arxiv.org/html/2608.25406v1#S2.SS1).
The implementation is in `explorer/Models/ElkiesSearchFamily.cs`; coefficient arrays
are in descending powers. It bundles S of degree 8, T of degree 12, and the 17
published x-coordinate polynomials of degree 4. No download is needed at runtime.

## Integral specialization

The published equation is `y² = x³ − 27 S(t)x + (27/4) T(t)`. For a reduced
`t = a/b`, b positive, let S_h, T_h and X_i,h be their homogeneous evaluations.
The change `x' = 4 b⁴ x`, `y' = 8 b⁶ y` yields

```
y'² = x'³ − 432 S_h(a,b) x' + 432 T_h(a,b).
x'_i = 4 X_i,h(a,b)
```

The y-coordinate is recovered as an exact nonnegative square root. A sign choice
does not change independence. Every supplied section must satisfy the specialized
equation. Denominator clearing uses integer arithmetic and does not factor the
discriminant or compute a minimal model. Singular specializations have no rank
claim and no section list.

`t = −9529/5471` is the published rank-28 specialization. The family package gives
only its 17 generic sections, so loading this parameter alone does not reproduce
the 28-point lower-bound certificate. It also does not make nearby parameters
inherit rank 28.

## Candidate score

The score is the finite sum `Σ log(#E(F_p)/p)` over primes `3 < p ≤ B` at which
the chosen integral model is nonsingular. Other primes contribute zero; no
minimalization is done to recover additional good reductions. For every prime,
the worker precomputes point counts for parameters in P¹(F_p), then looks them
up for each coprime pair a,b. Parameters with `p | b` use the point at infinity.
The tests compare these table scores with direct counts on specialized curves.
This is a finite heuristic, not an analytic-rank computation or a proof.
The general family-scoring approach is described by
[Elkies and Klagsbrun](https://arxiv.org/html/2003.00077v1).

Pairs are processed in denominator order in batches of 256. Only the highest
scores are retained, with deterministic ties by denominator and numerator.
New shortlisted curves receive exact section checks and the bounded additional
point search described in the Explorer guide. Concurrency is limited by the
chosen worker count. The implementation does not include a SIMD sieve or
reduced 2-covering search; it is an experimental home-computer starting point.

## What the rank certificate proves

`EllipticCurveQ.GetRankLowerBound(points, reductionPrimeBound, token)` returns
`PointRankCertificate`. It validates every supplied coordinate and clears rational
coefficient denominators by a rational model isomorphism. For an integral model,
the substitution `X=4x, Y=8y+4a1*x+4a3` gives

```
Y² = X³ + b2 X² + 8 b4 X + 16 b6.
```

At each tested good odd prime, each simple cubic root gives a quadratic-character
homomorphism into F₂. Binary elimination computes the image dimension of the
supplied points. Good reduction also bounds the rational 2-torsion dimension;
subtracting this upper bound gives a rigorous rank lower bound. A separate
infinite-order check preserves rank ≥1 when a non-torsion point has trivial
images in the tested quotients. No BSD, GRH, parity assumption, discriminant
factorization or full descent is used.

`IndependenceCertified` is true only when the lower bound equals the number of
supplied points (including duplicates and infinity in that count). These checks
neither prove saturation nor an upper bound. Increasing the prime limit can help
an inconclusive result, but cannot resolve every possible relation modulo 2.
For example, the bundled tests certify ≥15 at t=0,1,−1 and ≥17 at t=2/3 and
t=−9529/5471 using primes up to 1009. The smaller values are not exact-rank claims.

## Checkpoints and verification

`.ecsearch` stores version 1 JSON: family identifier, options, next pair index,
tested count and retained reports. Writes use a temporary file and atomic replace;
files and requests are size bounded. The most recent completed reported batch
survives Pause. Resume reconstructs and rechecks retained candidates. Saved rank
values are reports, not proof inputs. Range completion means the selected pair
box was scored, not that all rational points on those curves were found.

Regression coverage includes exact section specialization, dependent points,
rational model changes, cancellation, direct point-count comparisons, checkpoint
resume, invalid saved data, real worker execution, pause/restart controls and
compiled WPF layout/menu/import/undo checks.

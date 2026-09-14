# Equation-only search for a rank lower bound

`EllipticCurveQ.SearchRankLowerBound` takes a curve and resource limits. It does
not accept known witness points and does not query a database. It finds rational
points and returns an unconditional lower bound, with those points in the
original input coordinates. It uses managed C# only.

For verification of already available points, use `GetRankLowerBound(points)`.
Searching and verification have very different costs.

## API

```csharp
using System;
using System.Numerics;
using EllipticCurves;

var curve = new EllipticCurveQ(1, 0, 0,
    new BigRational(BigInteger.Parse("-1375414938269729933430")),
    new BigRational(BigInteger.Parse("20780863582673042643051322404516")));

var result = curve.SearchRankLowerBound(new RankLowerBoundSearchOptions
{
    NumeratorRadius = BigInteger.Parse("1200000000"),
    DenominatorRootBound = 1,
    TimeLimit = TimeSpan.FromSeconds(10),
    TargetLowerBound = 17
});
Console.WriteLine(result.LowerBound);
Console.WriteLine(result.Points.Count);
Console.WriteLine(result.StopReason);
```

Search coordinates have the form `x=m/d²`, with `gcd(m,d)=1`. The numerator
interval is inclusive: `NumeratorCenter-NumeratorRadius` through
`NumeratorCenter+NumeratorRadius`. Both center and radius are arbitrary-size
integers. Denominator roots `d=1,...,DenominatorRootBound` are visited in order;
numerators are visited in ascending order. A time limit can therefore stop
before reaching later denominators. No claim is made about unvisited coordinates.

For rational input coefficients, the method clears denominators using
`x'=scale²*x, y'=scale³*y`. Bounds apply to the resulting integral model and
`IntegralScale` reports this scale. Returned points are mapped back. This does
not require a global minimal model or integer factorization, but can produce
large search coordinates.

Defaults are numerator radius 2,000,000,000, center 0, denominator root bound 8,
10 seconds, 1,000,000 exact square tests, 1,024 retained points, and certificate
primes up to 1,009. A target lower bound is optional. These are resource settings,
not a promise that a particular positive rank will be detected.

The result distinguishes exhaustion of the specified box, elapsed time, square
test budget, retained point limit, and success at the requested target. Limits
return the already proved bound and its witnesses. Cancellation throws
`OperationCanceledException`, consistently with the library's existing APIs;
an optional `IProgress<RankLowerBoundSearchResult>` receives snapshots when the
bound increases. The elapsed-time limit is checked cooperatively, so a single
arithmetic operation or progress callback can overrun it.

## Sieve and exact proof

On an integral Weierstrass equation, define the usual invariants `b2,b4,b6`.
For `x=m/d²`, complete the square after clearing denominators:

```
F(m,d) = 4m³ + b2*m²*d² + 2b4*m*d⁴ + b6*d⁶
       = (2d³*y + a1*m*d + a3*d³)².
```

For each small odd prime through 251, precompute which residues of `m` make
`F` a quadratic residue, including zero. A 64-bit mask tests 64 consecutive
numerators at once. Blocks contain at most 65,536 candidates. Only candidates
surviving every mask undergo a coprimality check and exact integer square test.
The sieve also uses primes of bad reduction: quadratic-residue testing remains
a necessary condition there. The independent rank certificate only uses good
odd reductions, as required by its proof.

For an exact square, recover
`y=(sqrt(F)-a1*m*d-a3*d³)/(2d³)`. One point from each pair `P,-P` suffices for
rank lower bounds. Every point is checked by the certificate engine and mapped
to the original model. This does not assert that all returned points are independent.

The existing `RationalPointRank` engine incrementally builds a binary matrix of
good-reduction Kummer characters. Its image rank, minus an upper bound on the
rational 2-torsion dimension, is a rigorous rank lower bound. The engine also
retains its independent infinite-order check for a first point invisible to the
chosen characters. See [the certificate explanation](rank31-verification.md).

## Measurements on the two supplied curves

Observed on the local Windows machine with .NET 8.0.31, Release, on 2026-09-14.
Both equation-only searches used `|m| <= 1,200,000,000`, `d=1`, and 30 seconds.
The first stopped at target 17. Timings exclude process startup and build time.

| Curve | How points were obtained | Proved lower bound | Observed time |
|---|---|---:|---:|
| `[1,0,0,-1375414938269729933430,20780863582673042643051322404516]` | New C# search from the equation; 18 points found | 17 | 0.962 s |
| Same curve | Verify 17 previously found witnesses | 17 | 1.061 ms |
| ICARM #302, exactly the second supplied equation | New C# search from the equation; no integer points in this box | 0, inconclusive for the actual rank | 1.339 s |
| ICARM #302 | Verify its 31 published witnesses | 31 | 1.137 ms |

The first search prepared 25,880 sieve blocks and needed only 18 exact square
tests before proving the requested bound. It received no witness coordinates.
The second prepared 36,622 blocks and rejected every candidate before a square
test. This proves only the absence of integer points in the specified interval,
not absence of rational points elsewhere. These are individual timings, not
performance guarantees. Earlier PARI/GP exploration found 42 signed integer
points in the first box in 94 ms of reported CPU time.

The second curve has a 67-digit `a4`, a 99-digit `a6`, and a characteristic
coordinate scale `sqrt(|a4|)` of approximately `1.13e33`. Published witnesses
have enormous coordinates. Scanning up to that scale directly is impractical.
Its 31-point data already exist in `tests/Fixtures/icarm-302.json`; no known
rank, torsion or factorization metadata are trusted by the verification.
The independent Python verifier also obtained image rank 31 using primes up to
503 and a rootless good reduction at 31 to exclude rational 2-torsion.

Source and attribution: [ICARM curve #302](https://elliptic-rank.icarm.cloud/curve/302),
attributed there to Claude, Levent Alpöge and Ava Howell. ICARM is supported by
NSF Grant DMS 2425401. Witness verification here proves only rank at least 31.

## What a general large-curve solver still needs

The implemented method is a bounded first stage. It does not recover the 31
published points from the second equation. Increasing the original-coordinate
box or adding processors cannot by itself make a scan on the `1e33` scale practical.

A stronger equation-only solver can add exact coordinate changes with recorded
point maps; reduced 2-coverings and, where feasible, higher coverings; searches
with small denominators on these alternative models; and height-based subgroup
reduction. Coverings can make a large point on the original curve correspond to
a smaller point on the auxiliary curve. Computing useful coverings can itself
require costly arithmetic, so this does not give a uniform runtime guarantee.

Every successful search stage should feed the same incremental certificate
engine. A failed finite character test is inconclusive: independent even
multiples can be invisible modulo 2. Certified height-pairing matrices can
provide an additional independence method; an unqualified floating-point
determinant is not a rigorous substitute.

The covering and height-pairing search stages are not implemented by this API.
The current method provides exact, useful partial results and explicit stopping
states. It never promotes a heuristic score or an empty search to a proved rank.

For searching *new curves* rather than analysing an arbitrary input equation,
families with known sections and Mestre–Nagao scoring can guide expensive point
searches. This is a different input model from the equation-only method requested
here. See [Elkies and Klagsbrun, sections 2–5](https://arxiv.org/html/2003.00077v1)
for modular scoring, sieving and searching on coverings.

## Reproduce

The saved local experiment directory is `artifacts/rank-lower-bound-comparison`.
It contains two witness JSON files, a standard-library Python verifier, the
original PARI search, and a C# benchmark. The JSON for curve #302 explicitly
attributes its points to ICARM; the first JSON contains points found here.

```powershell
python artifacts/rank-lower-bound-comparison/verify.py
dotnet run --project artifacts/rank-lower-bound-comparison/Benchmark.csproj -c Release
dotnet test tests/EllipticCurves.Tests.csproj -c Release --filter FullyQualifiedName~RankLowerBoundSearchTests
```

The numerical tests compare the sieve with independent exact rational-point
enumeration, exercise rational coefficient maps, enormous coordinate centers,
word/block endpoints, torsion, inconclusive empty searches, limits and progress.

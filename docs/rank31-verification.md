# Verifying 31 independent points on the record curve

The library's existing `RationalPointRank` computation proves an unconditional
lower bound of 31 for [ICARM curve #302](https://elliptic-rank.icarm.cloud/curve/302).
The experiment uses the published equation and rational points, credited to
Claude, Levent Alpöge and Ava Howell on 2026-08-23. It does not search for those
points or compute an upper bound for the rank.

## Run the experiment

From the repository root, with the .NET SDK and the project's restored packages:

```powershell
dotnet test tests/EllipticCurves.Tests.csproj -c Release --filter FullyQualifiedName~RecordRankTests --logger "console;verbosity=detailed"
```

The [saved JSON](../tests/Fixtures/icarm-302.json) makes the verification offline;
no Sage, PARI, external arithmetic service or new library API is required. The
[test implementation](../tests/RecordRankTests.cs) uses the library's internal
certificate engine through the test assembly's existing access. No production
code changes were needed.

## What is proved

1. All 31 supplied points are distinct and satisfy the equation exactly, using
   `BigInteger` and `BigRational` arithmetic. The curve is nonsingular and its
   coefficients are integral.
2. There is no rational 2-torsion. With `X=4x` and
   `Y=8y+4a1*x+4a3`, the equation becomes
   `Y^2 = X^3 + b2*X^2 + 8*b4*X + 16*b6`.
   At the good prime 31 this cubic has no root, as checked by trying every residue.
   Since rational 2-torsion injects under good reduction at an odd prime,
   `E(Q)[2]` is zero. No factorization or database torsion claim is used.
3. The library evaluates Kummer characters at good odd primes. For every root
   `r` of the cubic modulo a prime `p`, the character of a point is the Legendre
   symbol of `X-r`; the value at `(r,0)` is given by the cubic's derivative,
   and the value at infinity is 1. Writing these signs as bits produces an exact
   binary matrix, whose rank is computed by elimination over `F_2`.
4. This matrix has rank 31. Its rank is at most the Mordell-Weil rank plus the
   rational 2-torsion dimension. Together with step 2, this proves `rank >= 31`
   and independence of the 31 supplied points modulo torsion.

The [two-descent notes](two-descent.md#lower-bounds-from-rational-points) explain
the character certificate used by the library. This experiment calls only that
part of the implementation. It never calls `GetRankBounds`, minimalization,
conductor computation, saturation or the full Selmer enumeration.

## Observed result

On 2026-09-13, with .NET 8.0.31 on the local Windows machine:

| Largest prime allowed | Binary matrix rank | Proved lower bound |
| --- | ---: | ---: |
| 101 | 10 | 10 |
| 251 | 25 | 25 |
| 503 | 31 | 31 |
| 1009 | 31 | 31 |

One run took approximately 19 ms for reading the data and all four verification
passes, excluding test discovery and .NET startup. Each test has a 30-second
cancellation limit. This is a recorded run, not a performance guarantee.

All four test cases passed. Replacing the last point with either a duplicate or
the sum of the first two points yields a bound of 30. Incrementing a supplied
y-coordinate by one is rejected because the point is no longer on the curve.

The conclusion is **rank at least 31**. This does not prove that the exact rank
is 31, that the subgroup is saturated, or that the supplied points form a full
Mordell-Weil basis. Neither BSD nor GRH is used in this verification.

Data attribution: the record and witness points are maintained by the NSF
Institute for Computer-Aided Reasoning in Mathematics (ICARM), supported by
NSF Grant DMS 2425401. The JSON snapshot was retrieved on 2026-09-13.

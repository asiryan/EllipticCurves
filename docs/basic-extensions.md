# Coefficients, CM, division, isogenies, prime fields and real logarithms

All native APIs below use managed C# and run offline. They do not call PARI, Sage,
Magma or LMFDB. Exact algebraic operations use integer/rational arithmetic.
The elliptic logarithm is explicitly a numerical result.

## Fourier coefficients and reduction counts

`GetFourierCoefficients(count, maxPointCountingWork, token)` returns a read-only
list indexed by n: the sentinel at index zero is 0, followed by a_1=1 through
a_count. Count zero is valid. `GetFourierCoefficient(n, ...)` requires n>=1 and
computes the list through n. All coefficients are exact signed 64-bit integers;
the implementation uses checked arithmetic.

`GetFrobeniusTrace(p, ...)` returns a_p. At a good prime,
`CountPoints(p, ...) = p+1-a_p`, including the point at infinity. At bad primes,
`GetFrobeniusTrace` returns the local L-series coefficient (0 or +/-1), while
`CountPoints` rejects the argument because reduction is not an elliptic curve.

These operations use the global minimal model, so rational changes of the input
equation preserve the output. Prime arguments are validated. Counting is direct,
with a quadratic-residue sieve for odd primes and a separate characteristic-two
case. Coefficients at composite indices use multiplicativity and the good/bad
prime recurrences. The default work limit is 20000000. Exceeding it throws rather
than returning a truncated list. Storage grows with the requested coefficient
count and the largest counted prime; raising work limits also permits larger
allocations. Limits do not bound minimalization or factorization time.

## Complex multiplication over Q

`HasComplexMultiplication` and `CmDiscriminant` use the complete table of thirteen
rational CM j-invariants. The discriminant refers to the endomorphism order over
the algebraic closure, not to endomorphisms all being defined over Q. Zero means
non-CM. Quadratic twists and rational model changes preserve the answer.
The API does not construct endomorphism maps or handle general number fields.

The finite-table method and its scope are explained in
[Cremona, Computing the endomorphism ring of an elliptic curve over a number field](https://lucant.org/papers/2023/230123-Cremona.pdf).

## Exact rational division points

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0);
var p = new EllipticCurvePoint(0, 0);
var options = new PointDivisionOptions { MaxWork = 2000000, MaxDivisionDegree = 1024 };
var preimages = e.GetDivisionPoints(e.Multiply(p, 6), 2, options);
bool divisible = e.TryDividePoint(p, 2, out var quotient, options); // false
```

`GetDivisionPoints(P,n)` finds all Q in E(Q) such that [n]Q=P, including torsion
translates, on the original model. Signed n is supported, except zero and
`int.MinValue`. For P=O it filters the full rational torsion subgroup. For finite
P it solves the exact multiplication equation on a short model, enumerates all
rational roots and validates both ordinate signs with the exact group law.
An empty result proves nondivisibility. `TryDividePoint` selects a preimage when
one exists; its output is an unused infinity sentinel when it returns false.

The degree limit applies to the n^2 multiplication equation. Counted polynomial
and root-isolation work is limited separately. An incomplete search throws
`ArithmeticException`, so callers cannot mistake exhaustion for nondivisibility.
Cancellation is cooperative. Existing torsion preparation is outside the counter
and is not internally interruptible.

## Velu isogenies and dual 2-isogenies

`CreateIsogeny(kernelGenerators)` accepts generators of a finite subgroup of
E(Q), including noncyclic subgroups. It computes the full kernel, then applies
Velu's formulas after an exact change to a short equation with u=1. Infinite-order
generators are rejected by exact addition and the rational torsion-order bound.
The trivial kernel gives a degree-one map to the short model.

`RationalIsogeny` exposes `Source`, `Target`, `Degree`, all `Kernel` points and
`Map(P)`. The target is a short equation, not necessarily globally minimal.
Kernel points map to infinity; invalid source points are rejected. Point maps
use exact rational arithmetic and each output is checked against the target.

For y^2=x^3+Ax+B, summing over all nonzero kernel points Q gives

    t = sum (3*x(Q)^2+A)
    w = sum (2*y(Q)^2+x(Q)*(3*x(Q)^2+A))
    A' = A-5*t, B' = B-7*w.

The coordinate maps are the sums x(P)+sum(x(P+Q)-x(Q)) and similarly for y.
The formulas and normalization follow
[Sutherland's MIT elliptic-curve lectures, Lecture 5](https://math.mit.edu/classes/18.783/2023/LectureSlides5.pdf).

`CreateTwoIsogeny(T)` requires a nonzero rational point of order two and returns a
`TwoIsogenyPair` with `Forward` and `Dual`. Both compositions equal multiplication
by two on their respective models, including the exact original input model.
For a short-model kernel root alpha, the dual kernel root is -2*alpha. A second
normalized Velu quotient has coefficients 16*A and 64*B; scaling by u=2 and
undoing the original short-model change gives the dual with the correct sign.
`GetTwoIsogenies()` returns one such pair for each rational nonzero 2-torsion point.

This API requires a pointwise rational kernel. It does not accept a kernel
polynomial with nonrational roots, discover all rational isogenies, or construct
a full isogeny graph. Only degree-two duals are exposed.

## Curves over prime fields

```csharp
var e = new EllipticCurveFp(7, 0, 0, 0, -1, 0);
var p = e.CreatePoint(0, 0);
var twice = e.Double(p); // infinity
var count = e.CountPoints();
var points = e.Points(); // lazy complete enumeration, infinity first
```

`EllipticCurveFp` supports general Weierstrass equations over a proved prime
modulus, including 2 and 3. Coefficients and point coordinates are canonical
residues. Singular equations and composite moduli are rejected. Primality testing
uses the existing exact native machinery and can be costly beyond 64 bits.

The public operations are `CreatePoint`, `IsOnCurve`, `Negate`, `Add`, `Subtract`,
`Double`, `Multiply`, `CountPoints`, `Points` and `GetPointOrder`, together with
the discriminant, c4, c6 and j. Affine `EllipticCurvePointFp` values carry their
field modulus; points from a different field are rejected. The universal infinity
is `EllipticCurvePointFp.Infinity`; a default-initialized struct is not a valid point.

Counting and enumeration visit every x-coordinate; the default limit is 1000000
x-coordinate steps. Odd-characteristic enumeration uses Tonelli-Shanks square
roots. `GetPointOrder` first counts the group, then factors its order and removes
prime factors by exact scalar multiplication. The counter is not a bound on
individual modular operations, elapsed time or factorization. These are basic
algorithms, with no SEA or extension-field support.

On `EllipticCurveQ`, `ReduceModuloPrime(p)` reduces the global minimal model at a
good prime. `ReducePointModuloPrime(P,p)` first maps P onto that model, then reduces
it. Poles of the minimal x-coordinate map to infinity. Both reject bad reduction.

## Numerical elliptic logarithm of a rational point

`RealEllipticLogarithm(P, options, token)` returns `RealEllipticLogarithmResult`.
It uses the same minimal model, differential and period convention as `GetPeriods`.
It accepts rational points, which are real points; general complex input is not
supported. `RealComputationOptions` controls certified root/period preparation and
iteration limits. It does not give the final double logarithm a certified error bound.

- `RealPart` lies in [0,w1), subject to floating-point rounding.
- `PrimitiveRealPeriod` is w1, the least positive real period, not the BSD period
  multiplied by the number of real components.
- `ComponentIndex` is zero on the component containing infinity and one on the
  bounded component when it exists.
- `ImaginaryPart` is zero on the identity component and Im(w2)/2 on the bounded
  component. Thus the two coordinates describe a period-class representative.
- Infinity has logarithm zero. The sign of 2*y+a1*x+a3 selects the branch.

The method isolates the roots of the completed-square cubic exactly, and evaluates
real incomplete elliptic integrals through Carlson R_F. Duplication and the
degree-seven expansion are from [DLMF 19.26.18](https://dlmf.nist.gov/19.26.E18)
and [DLMF 19.36.1](https://dlmf.nist.gov/19.36.E1). Exact rational interval preparation
avoids subtracting nearly equal floating-point roots. Insufficient relative root
or period precision, double range failures and exhausted iterations throw explicitly.
The numerical logarithm itself is not used to certify ranks, heights or saturation.

## Additional stored LMFDB metadata

These properties read the same cached aggregate snapshot, with no extra HTTP call:

| Property | Stored field |
|---|---|
| `CremonaLabel`, `IsogenyClassLabel`, `IsogenyClassSize` | `Clabel`, `lmfdb_iso`, `class_size` |
| `CmDiscriminant` | `cm` (0 means non-CM; null means missing) |
| `IsogenyDegrees`, `IsogenyMatrix` | `isogeny_degrees`, `ec_classdata.isogeny_matrix` |
| `FourierCoefficients`, `PrimeFourierCoefficients` | `ec_classdata.anlist`, `aplist` |
| `ModularDegree`, `ManinConstant`, `TorsionOrder` | `degree`, `manin_constant`, `torsion` |
| `FaltingsHeight`, `StableFaltingsHeight` | `faltings_height`, `stable_faltings_height` |
| `BsdShaOrder`, `AnalyticShaOrder` | `sha`, `ec_mwbsd.sha_an` |
| `LeadingLValue` | `ec_mwbsd.special_value`, the Taylor coefficient L^(r)(1)/r! |
| `IntegralPointXCoordinates` | `ec_mwbsd.xcoord_integral_points` on the stored minimal model |

Optional missing fields stay null; arrays and matrix rows are read-only.
Coefficient indexing, class labels, matrix shape and consistency are checked.
`PrimeFourierCoefficients` follows consecutive primes 2,3,5,...; `FourierCoefficients`
includes the index-zero sentinel. Stored isogeny matrices do not supply point maps.

The Sha fields represent the stored analytic/BSD value, not a native calculation
or unconditional proof of the group order; see
[LMFDB's analytic Sha convention](https://www.lmfdb.org/knowledge/show/ec.q.analytic_sha_order).
Likewise the stored integral-point coordinates do not certify native search completeness.

## Verification

Tests compare coefficients and CM metadata with offline LMFDB records, enumerate
small prime-field curves independently in both coordinates, check reduction and
group laws, and check every division point by exact multiplication. Isogeny tests
include cyclic/noncyclic kernels, both dual compositions and rational model changes.
Additional PARI fixtures contain 99 real logarithms and five isogeny quotients;
logarithm tests also check addition modulo periods and both real components.
See [fixture generation instructions](../tests/Fixtures/README.md).

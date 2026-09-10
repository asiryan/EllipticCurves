<p align="center"><img width="25%" src="docs/png/ec_logo_v3b.png" /></p>

# About
**EllipticCurves** is a small C# library for studying and working with elliptic curves. It provides functionality to compute and explore:
* coefficients and group structure,  
* discriminant and j-invariant,  
* torsion/rational/integral points,  
* short and minimal Weierstrass model,  
* proved algebraic rank bounds by general 2-descent or descent by 2-isogeny,
* native analytic rank estimates, with rigorous certificates for ranks 0 and 1,
* exact model isomorphisms and point maps,
* local reduction types, Kodaira symbols and Tamagawa numbers,
* certified canonical/local heights, height pairings and subgroup regulators,
* certified real/complex periods and period-lattice area,
* exact subgroup saturation at explicitly requested primes,
* public Frobenius traces, Fourier coefficients and good-reduction point counts,
* exact CM recognition over the rationals,
* exact rational division points and Velu isogenies with rational kernels,
* explicit 2-isogenies and their dual maps,
* curves and point arithmetic over prime fields, including characteristics 2 and 3,
* numerical elliptic logarithms of rational points on both real components,
* algebraic/analytic ranks from optional LMFDB metadata,
* LMFDB label/url,  
* conductor, etc.  

# Version
You can build **EllipticCurves** from sources or install to your own project using nuget package manager.
| Assembly | Specification | OS | Download | Package |
|-------------|:-------------:|:-------------:|:--------------:|:--------------:|
| [EllipticCurves](sources) | .NET Standard 2.0 | Cross-platform | [Release](https://github.com/asiryan/EllipticCurves/releases/) | [NuGet](https://www.nuget.org/packages/EllipticCurves/) | 

# Installation
C# interface  
```c#
using EllipticCurves;
```
To get started with **EllipticCurves** it is recommended to take a look at the [example project](examples).  
Here are some results for the [elliptic curve](https://arxiv.org/abs/2510.11768): **Y^2 = X^3 - 17X^2 + 72X**.
```
E: y^2 = x^3 - 17*x^2 + 72*x
Short Weierstrass: y^2 = x^3 - 73/3*x + 1190/27
Torsion: Z/2Z x Z/4Z
b2 = -68
b4 = 144
b6 = 0
b8 = -5184
D  = 82944
c4 = 1168
c6 = -38080
j  = 1556068/81
Torsion points:
O
(0, 0)
(8, 0)
(9, 0)
(6, 6)
(6, -6)
(12, 12)
(12, -12)
LMFDB: 48.a3
Url: https://www.lmfdb.org/EllipticCurve/Q/48.a3/
Minimal Weirstrass model: y^2 = x^3 + x^2 - 24*x + 36
Torsion: Z/2Z x Z/4Z
Rank(E) = 0
Analytic rank(E) = 0
Cond(E) = 48
Isomorphic to E: True
Native minimal Weierstrass model: y^2 = x^3 + x^2 - 24*x + 36
Native rank bounds(E) = 0
Exact native rank proved: True
Native Cond(E) = 48
Native minimal model matches LMFDB: True
Native conductor matches LMFDB: True
LMFDB rank is within native bounds: True
```

## Rank bounds and conductor

`EllipticCurveQ` computes minimal models, conductors and rank bounds in C# using
exact integer/rational arithmetic. No database, native binary, Sage or PARI installation
is required. `LmfdbEllipticCurve` remains available separately for optional online metadata.

```csharp
using System.Numerics;
using EllipticCurves;

var e = new EllipticCurveQ(0, -17, 0, 72, 0);
var minimal = e.GlobalMinimalModel;
BigInteger conductor = e.Conductor;       // 48
RankBounds bounds = e.GetRankBounds();
Console.WriteLine(bounds.LowerBound);     // 0
Console.WriteLine(bounds.UpperBound);     // 0
Console.WriteLine(bounds.ExactRank);      // 0 (null unless proved)
```

The conductor is computed on a global minimal model, with Tate's algorithm handling
wild reduction at 2 and 3. Rational coefficients and non-minimal input models are supported.

Rank bounds use descent by 2-isogeny when a rational point of order 2 is available,
and general binary-quartic 2-descent otherwise. Local checks decide solubility over
the reals and the necessary p-adic fields. Rational points, distinct soluble covering
classes and exact good-reduction characters prove lower bounds, including ranks above 1.
These computations do not assume BSD, GRH, parity, or finiteness of Sha.

```csharp
var e = new EllipticCurveQ(0, 0, 1, -7, 6);
var rank = e.GetRankBounds();
Console.WriteLine(rank.ExactRank);           // 3
Console.WriteLine(rank.UsedGeneralTwoDescent); // True
Console.WriteLine(rank.TwoSelmerDimension);  // 3
```

A complete descent gives a proved upper bound, which can exceed the rank because
of Sha. An unsuccessful point search never proves rank zero. A work limit during
descent gives `UpperBound = null` and an explanation in `Reason`; any proved lower
bound is retained. `ExactRank` is populated only when both bounds agree.

`GetRankBounds(searchBound: 64)` increases the point search. `maxSquareClasses` defaults
to 65,536 per isogeny; exceeding it throws rather than silently truncating the descent.
Factorization and searches can be expensive. Use the cancellable methods when needed:

```csharp
using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
var conductor = e.GetConductor(timeout.Token);
var minimal = e.GetGlobalMinimalModel(timeout.Token);
var rank = e.GetRankBounds(searchBound: 64, cancellationToken: timeout.Token);
```

For more control, pass `RankComputationOptions`:

```csharp
var rank = e.GetRankBounds(new RankComputationOptions
{
    SearchBound = 64,
    MaxDescentWork = 10000000,
    MaxPointSearchWork = 2000000
}, timeout.Token);
Console.WriteLine(rank.Reason);
```

General descent enumerates a complete reduction region; its cost can grow rapidly
with the curve invariants. The work limits bound counted steps, not elapsed time or
integer factorization. Point-search exhaustion preserves a completed upper bound.
`PreferGeneralTwoDescent` also enables the general method for curves with 2-torsion.
See [the descent construction and proof conditions](docs/two-descent.md).

## Analytic rank and BSD

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0);
var result = e.EstimateAnalyticRank();
Console.WriteLine(result.EstimatedRank);  // 1
Console.WriteLine(result.Status);         // Certified
Console.WriteLine(result.ProvenRank);     // 1; null for an unproved numerical estimate
Console.WriteLine(result.Derivatives[1]); // L'(E,1), approximately 0.3059997738340523
Console.WriteLine(e.RootNumber);          // -1 (exact)
```

`EstimateAnalyticRank` computes Fourier coefficients and central L-function derivatives
in C#, without HTTP or external mathematical software. BSD predicts that the order
of vanishing equals the algebraic rank. Results distinguish three cases:

- `Certified`: a rigorous interval proves analytic rank 0 or 1. Known theorems then
  prove the same algebraic rank, so this status does not assume BSD.
- `NumericalEstimate`: small derivatives have only been recognized numerically.
  BSD alone does not turn numerical zero recognition into a proof.
- `Inconclusive`: a work limit, numerical ambiguity, or the derivative-order limit
  prevents a rank estimate. `EstimatedRank` and `ProvenRank` are null.

`AnalyticRankOptions` controls the derivative order (default 4, maximum 8), zero
threshold, number of coefficients, point-counting work, integration evaluations and
interval-certificate length. Numerical derivatives use `double`; the separate
rank 0/1 certificates use exact outward-rounded dyadic intervals and rigorous tails.
`EstimatedErrors` are numerical diagnostics, not proof bounds. Large conductors are
expensive; this implementation uses direct point counting, not SEA. Heights,
regulators, periods and Tamagawa numbers have separate native APIs below.
The Tate–Shafarevich group order is not computed.

```csharp
var analytic = e.EstimateAnalyticRank(new AnalyticRankOptions
{
    MaxDerivativeOrder = 6,
    MaxTerms = 30000,
    MaxPointCountingWork = 50000000,
    ZeroTolerance = 1e-9
}, cancellationToken: timeout.Token);
```

## Local data, heights, periods and saturation

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0); // 37.a1
var generator = new EllipticCurvePoint(0, 0);
var modelMap = e.GetMinimalModelIsomorphism();
var local = e.GetLocalData(37); // I1, nonsplit multiplicative, c_37 = 1
var height = e.CanonicalHeight(generator); // approximately 0.05111140824
var regulator = e.Regulator(new[] { generator });
var periods = e.GetPeriods(); // RealPeriod approximately 5.98691729246

var saturation = e.Saturate(
    new[] { e.Multiply(generator, 6) }, new[] { 2, 3 });
// On completion: IndexGain = 6, CertifiedPrimes = [2, 3].
// IsComplete certifies only the requested primes, not a full Mordell-Weil basis.
```

Real results expose exact rational `LowerBound` and `UpperBound`; `Approximation`
is for display. `RealComputationOptions` controls absolute accuracy and work limits.
`Regulator(points)` refers to the supplied subgroup modulo torsion. Saturation
requires independent input generators and reports unresolved primes when limited.

The optional LMFDB adapter exposes matching stored data:

```csharp
var stored = await LmfdbEllipticCurve.FetchAsync(e, cancellationToken: timeout.Token);
var storedGenerators = stored.GetGeneratorsOnModel(e);
var storedLocalData = stored.LocalData;
var storedHeights = stored.GeneratorHeights;
var storedRegulator = stored.Regulator;
var storedRealPeriod = stored.RealPeriod;
var storedArea = stored.PeriodArea;
```

The synchronous constructor remains available. Missing optional fields are `null`;
stored decimal values retain their text and precision metadata, but are not certified
intervals. `FromStoredDataJson(json)` reads a downloaded curve-data snapshot offline.

See [height, period, saturation and LMFDB notes](docs/heights-and-saturation.md)
and [rank algorithm notes](docs/native-arithmetic.md) for conventions and limits.

## Coefficients, CM, division points and isogenies

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0); // 37.a1
var p = new EllipticCurvePoint(0, 0);
var ap = e.GetFrobeniusTrace(5);             // -2
var coefficients = e.GetFourierCoefficients(100); // a[n], with a[0]=0
var count = e.CountPoints(5);               // 8, including infinity
var divided = e.GetDivisionPoints(e.Multiply(p, 6), 2); // exactly {3P}
var logarithm = e.RealEllipticLogarithm(p); // numerical, on the minimal model

var cm = new EllipticCurveQ(0, 0, 0, 0, 1);
Console.WriteLine(cm.CmDiscriminant);        // -3; zero denotes non-CM
var threeIsogeny = cm.CreateIsogeny(new[] { new EllipticCurvePoint(0, 1) });
Console.WriteLine(threeIsogeny.Degree);      // 3
var two = cm.CreateTwoIsogeny(new EllipticCurvePoint(-1, 0));
// two.Dual.Map(two.Forward.Map(Q)).Equals(cm.Double(Q)) is true.

var finite = e.ReduceModuloPrime(5);
var reducedPoint = e.ReducePointModuloPrime(p, 5);
var order = finite.GetPointOrder(reducedPoint);
var allPoints = finite.Points();
```

Division returns every rational preimage or an empty list proving nondivisibility;
exhausting a work limit throws. Isogeny kernels in this API consist of rational
points. General isogeny-class discovery and global saturation are separate,
unimplemented tasks. Prime-field counting uses direct search, not SEA.

The LMFDB adapter additionally reads stored Fourier coefficients, CM discriminants,
isogeny degrees/matrices, modular degrees, Manin constants, Faltings heights,
analytic Sha values, leading L-values and integral-point x-coordinates. These are
database records; in particular the Sha fields do not assert a native proof.
See [the new API conventions and examples](docs/basic-extensions.md).

Run the example and regression tests with:

```sh
dotnet run --project examples/EllipticCurves.Example.csproj
dotnet test EllipticCurves.sln
```

The example retains its LMFDB lookup and checks the native results against it,
so running the example requires internet access. The native API and regression
tests work offline.

# License
**MIT**  

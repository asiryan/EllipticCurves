<p align="center"><img width="25%" src="docs/png/ec_logo_v3b.png" /></p>

# About
**EllipticCurves** is a C# library for studying elliptic curves over the rationals and finite fields. It provides functionality to compute and explore:
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
* certified Faltings and stable Faltings heights,
* exact subgroup saturation at explicitly requested primes,
* public Frobenius traces, Fourier coefficients and good-reduction point counts,
* exact CM recognition over the rationals,
* exact rational division points and Velu isogenies with rational kernels,
* explicit 2-isogenies and their dual maps,
* curves and point arithmetic over prime fields, including characteristics 2 and 3,
* finite extensions F_(p^k), with exact irreducibility checks and curve arithmetic,
* numerical elliptic logarithms of rational points on both real components,
* algebraic/analytic ranks from optional LMFDB metadata,
* LMFDB label/url,  
* conductor and its prime factorization.

# Version
Build **EllipticCurves** from source or install the NuGet package in your project.
| Assembly | Specification | OS | Download | Package |
|-------------|:-------------:|:-------------:|:--------------:|:--------------:|
| [EllipticCurves](sources) | .NET Standard 2.0 | Cross-platform | [Release](https://github.com/asiryan/EllipticCurves/releases/) | [NuGet](https://www.nuget.org/packages/EllipticCurves/) | 

The NuGet package contains the library. The [Console application](console/README.md)
and [desktop Explorer](explorer/README.md) are distributed separately. Their Windows
release archives include .NET and run without a separate runtime installation.
See [release preparation](docs/releasing.md) for version
settings, validation and packaging commands.

# Installation
```shell
dotnet add package EllipticCurves
```

This installs the latest published stable package. The examples below describe
the current source checkout, which may be newer. For changes not yet published
to NuGet, reference the library project or build a local package using
[the packaging instructions](docs/releasing.md#nuget-package).

Import the namespace in your C# code:
```csharp
using EllipticCurves;
```

See the [Console application](console/README.md) for command-line arguments, run instructions
and sample output.

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
    MaxPointSearchWork = 2000000,
    MaxDegreeOfParallelism = 4 // optional; the library defaults to 1 (sequential)
}, timeout.Token);
Console.WriteLine(rank.Reason);
```

General descent enumerates a complete reduction region; its cost can grow rapidly
with the curve invariants. The work limits bound counted steps, not elapsed time or
integer factorization. Point-search exhaustion preserves a completed upper bound.
`PreferGeneralTwoDescent` also enables the general method for curves with 2-torsion.
`MaxDegreeOfParallelism` parallelizes general binary-quartic enumeration only;
the 2-isogeny method and other rank stages remain sequential. All workers share
the same work limits. Cancellation and work exhaustion stop and join the workers
before returning. Parallel scheduling may change covering representatives, work
counts and partial lower bounds; an upper bound still requires complete descent.
See [the descent construction and proof conditions](docs/two-descent.md).

## Rank lower bounds from supplied points

`GetRankLowerBound` verifies supplied rational points and proves a lower bound
using exact good-reduction characters, without factorization, minimalization
or full descent:

```csharp
var e = new EllipticCurveQ(0, 0, 1, -7, 6);
var certificate = e.GetRankLowerBound(new[]
{
    new EllipticCurvePoint(0, 2),
    new EllipticCurvePoint(1, 0),
    new EllipticCurvePoint(2, 0)
});
Console.WriteLine(certificate.LowerBound);             // 3
Console.WriteLine(certificate.IndependenceCertified);   // True
Console.WriteLine(certificate.Reason);
```

Points must use the input model's coordinates. `reductionPrimeBound` controls
the tested primes (default 1009); cancellation is supported. The returned
`PointRankCertificate` reports the proved bound, point count, independence status,
character image dimension and reduction-prime diagnostics. A smaller bound does
not prove dependence, and zero does not prove rank zero. Neither an upper bound
nor saturation is asserted. See [the 31-point verification](docs/rank31-verification.md)
for the certificate construction and a larger example.

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
See [API conventions and examples](docs/basic-extensions.md).

## Faltings heights and finite extensions

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0);
var faltings = e.FaltingsHeight();
var stableFaltings = e.StableFaltingsHeight();
// Both return certified RealEnclosure values in the LMFDB normalization.
// FaltingsHeight uses the global minimal model, even for nonminimal input.
```

Finite extensions use a specified irreducible polynomial, with coefficients in
ascending order. Primality and irreducibility are proved when constructing the field.

```csharp
using System.Numerics;

var field = new FiniteField(3, new BigInteger[] { 1, 0, 1 }); // F_9, t^2+1
var alpha = field.Generator;
Console.WriteLine((alpha * alpha + 1).IsZero); // True

var e = new EllipticCurveFq(field,
    field.Zero, field.Zero, field.Zero, alpha, field.One);
var p = e.CreatePoint(0, 1);
var twice = e.Double(p);
var count = e.CountPoints();
var points = e.Points();
```

`EllipticCurveFq` supports general equations, including characteristics 2 and 3.
Counting and enumeration check all q^2 affine coordinate pairs; the default limit
is 1000000 pairs. Field presentations must agree before their elements can be mixed.
See [Faltings-height conventions and finite-extension limits](docs/faltings-and-finite-fields.md).

## Tests

From the repository root, with the .NET 8 SDK installed, run the portable
regression tests with:

```sh
dotnet test tests/EllipticCurves.Tests.csproj -c Release
```

The regression tests work offline after dependencies have been restored. The Windows-only WPF
regression check is a separate command documented in
[the Explorer development guide](explorer/README.md#development).

## Desktop Explorer (Windows)

The [WPF Explorer](explorer/README.md) provides a modern desktop interface on
.NET 8: full formula input for simple and general Weierstrass equations, optional
exploration sliders, an interactive real-locus plot, a linked period-lattice and
3D complex-torus view, exact invariants and bounded rational-point samples.
The Tools menu exposes the library's computations,
including torsion, ranks, Faltings heights, periods, isogenies and finite fields.
**Verify rank from supplied points** accepts exact coordinates and reports
a proved rank lower bound with its independence certificate.
Parameter windows feed a results panel with session history, progress, cancellation
and time limits. Native computations run locally; only the explicit LMFDB search and fetch
commands require internet access.

# License

**MIT**  

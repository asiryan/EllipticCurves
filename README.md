<p align="center"><img width="25%" src="docs/ec_logo.png" /></p>
<p align="center"> .NET library for studying and working with elliptic curves </p>    

# About
**EllipticCurves** is a small C# library for studying and working with elliptic curves. It provides functionality to compute and explore:
* coefficients and group structure,  
* discriminant and j-invariant,  
* torsion/rational/integral points,  
* short and minimal Weierstrass model,  
* proved algebraic rank bounds (exact when the bounds coincide),
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

# Offline arithmetic

`EllipticCurveQ` computes minimal models, conductors and rank bounds in C# using
exact integer/rational arithmetic. No database, native binary, Sage or PARI installation
is required. `EllipticCurveLMFDB` remains available separately for optional online metadata.

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

Rank bounds use descent by 2-isogeny when the curve has a rational point of order 2.
Quartic point searches prove lower bounds; real and modular obstructions prove upper
bounds. These are unconditional bounds, with no BSD or parity assumption. The finite
local sieve may leave an interval even when a full Selmer computation could resolve it.

For curves **without rational 2-torsion**, the current implementation proves only a
lower bound of 0 or 1 by point search; `UpperBound` and `ExactRank` are `null`.
General 2-descent and native analytic rank computation are not implemented. An unsuccessful
point search never proves rank zero.

`GetRankBounds(searchBound: 64)` increases the point search. `maxSquareClasses` defaults
to 65,536 per isogeny; exceeding it throws rather than silently truncating the descent.
Factorization and searches can be expensive. Use the cancellable methods when needed:

```csharp
using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
var conductor = e.GetConductor(timeout.Token);
var minimal = e.GetGlobalMinimalModel(timeout.Token);
var rank = e.GetRankBounds(searchBound: 64, cancellationToken: timeout.Token);
```

See [algorithm notes](docs/native-arithmetic.md) for the bounds and limitations.
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

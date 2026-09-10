# EllipticCurves

A managed C# library for elliptic curves over the rationals and prime fields,
targeting .NET Standard 2.0. Native computations run offline without PARI, Sage,
Magma or other mathematical software. LMFDB integration is optional.

```csharp
using EllipticCurves;

var e = new EllipticCurveQ(0, 0, 1, -1, 0); // 37.a1
var p = new EllipticCurvePoint(0, 0);

var twice = e.Double(p);
var conductor = e.Conductor;                // 37
var bounds = e.GetRankBounds();
var rank = bounds.ExactRank;                 // null unless proved
var local = e.GetLocalData(37);
var height = e.CanonicalHeight(p);           // certified rational enclosure
var periods = e.GetPeriods();                // on the global minimal model
var count = e.CountPoints(5);                // 8, including infinity
var finite = e.ReduceModuloPrime(5);
```

## Capabilities

- Exact rational arithmetic, invariants, point arithmetic, rational torsion,
  minimal models, isomorphisms and quadratic twists.
- Conductors, root numbers, local reduction data and Tamagawa numbers.
- General 2-descent and descent by 2-isogeny, with proved algebraic rank bounds.
- Numerical central L-function derivatives and analytic rank estimates, with
  separate rigorous certificates for analytic ranks 0 and 1.
- Certified canonical/local heights, subgroup regulators and periods.
- Subgroup saturation at explicitly requested primes and exact rational division.
- Fourier coefficients, Frobenius traces, rational CM recognition, Velu isogenies
  with pointwise rational kernels, and dual 2-isogenies.
- Prime-field curves, arithmetic, direct point counting and enumeration,
  including characteristics 2 and 3.
- Numerical elliptic logarithms of rational points on both real components.
- Optional cached LMFDB data, including generators, heights, periods and
  isogeny-class metadata, with offline snapshot parsing.

## Scope and limits

Point searches are bounded. Saturation certifies only the requested primes and
does not establish a full Mordell-Weil basis. Work limits can leave rank bounds
incomplete. Numerical estimates and stored LMFDB values are distinguished from
proofs; the library does not compute the Tate-Shafarevich group order.
Prime-field counting uses direct search, not SEA. General isogeny-class discovery,
higher descents and curves over general number fields are outside this release.

Version 3 uses `LmfdbEllipticCurve` in place of `EllipticCurveLMFDB`; the old name
has been removed. Callers using the old name must update it.

See the [repository README](https://github.com/asiryan/EllipticCurves/blob/main/README.md)
for examples and the algorithm notes for
[ranks](https://github.com/asiryan/EllipticCurves/blob/main/docs/native-arithmetic.md),
[heights, periods and saturation](https://github.com/asiryan/EllipticCurves/blob/main/docs/heights-and-saturation.md),
and [coefficients, division, isogenies and prime fields](https://github.com/asiryan/EllipticCurves/blob/main/docs/basic-extensions.md).

Licensed under MIT.

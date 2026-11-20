# Faltings heights and curves over finite extensions

These two additions use managed C# only. LMFDB data and PARI outputs are used for
comparison in tests, with no native-binary or network dependency at runtime.

## Certified Faltings heights over Q

`FaltingsHeight(options, token)` and `StableFaltingsHeight(options, token)` return
`RealEnclosure`, using the same `RealComputationOptions` as heights and periods.
Successful results have exact rational endpoints and width at most
10^(-DecimalDigits). `Approximation` is a double for display, not a proof bound.

Both methods first use the reduced global minimal model. Write A for the area of
its Neron period lattice, Delta for its discriminant and d for the positive
denominator of its j-invariant. The formulas are

    h_F = -log(A)/2
    h_F_stable = h_F + log(d/|Delta|)/12.

This matches the [LMFDB Faltings-height normalization](https://www.lmfdb.org/knowledge/show/ec.q.faltings_height).
The correction is also documented in
[Sage's Faltings-height API](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_rational_field.html#sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.faltings_height).
Some conventions use the area on the supplied equation, which changes on scaling
the equation. Our `FaltingsHeight` always uses the minimal model, consistently with
`GetPeriods` and LMFDB, and so both public methods are invariant under Q-model
changes. The stable height is also invariant under quadratic twisting. For a
semistable minimal curve d=|Delta|, so the two heights agree.

The area is enclosed by the existing root-isolation/AGM period algorithm. Logarithms,
the exact rational correction and all subsequent operations round outwards.
The final width is checked after taking the logarithm; an inadequate period
enclosure or insufficient precision throws instead of silently lowering accuracy.
Singular curves throw `InvalidOperationException`; invalid options throw
`ArgumentOutOfRangeException`; exhausted root/AGM work or inadequate accuracy
throws `ArithmeticException`. Cancellation throws `OperationCanceledException`
at cooperative checks. Null options select defaults: 12 decimal digits, 256
precision bits, 512 maximum iterations and 100000 root-isolation steps. The API
does not automatically increase precision or retry after a work-limit failure.
Neither height is constrained to be positive.

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0);
var options = new RealComputationOptions { DecimalDigits = 30, PrecisionBits = 512 };
var height = e.FaltingsHeight(options);
var stable = e.StableFaltingsHeight(options);
```

## Exact finite fields

`FiniteField(p, definingPolynomial, token)` represents F_p[t]/(f), including degree
one. The coefficient list is ascending: `[1,0,1]` means 1+t^2. Construction reduces
coefficients modulo p, removes trailing zeros and makes f monic. Constant or zero
polynomials and reducible polynomials throw `ArgumentException`; a nonprime
characteristic throws `ArgumentOutOfRangeException`. A null coefficient list
throws `ArgumentNullException`.

The characteristic is proved prime with the existing native number-theory routines.
Irreducibility uses the exact Rabin criterion: for degree k, t^(p^k)=t modulo f,
and gcd(f,t^(p^(k/r))-t)=1 for every prime r dividing k. Powers are calculated by
successive Frobenius steps modulo f; the potentially enormous polynomial t^(p^k)
is never materialized. See the finite-field algorithms and irreducibility criterion
in [Shoup, A Computational Introduction to Number Theory and Algebra](https://www.shoup.net/ntb/ntb-v2_5.pdf).
Construction is deterministic once the prime proof completes; it does not accept
probable irreducibility. There is no automatic search for a defining polynomial.

`Characteristic`, `Degree`, `Order` and `Modulus` describe the field. `Zero`, `One`
and `Generator` are elements; Generator is the residue class of t and is not
necessarily a generator of the multiplicative group.

`CreateElement(c0,c1,...)` reduces the given polynomial modulo p and f. A
`FiniteFieldElement` has immutable canonical coefficients of degree less than k,
without trailing zeros. Zero has an empty coefficient list. The field and element
coefficient collections are read-only and input arrays are copied.

Elements support addition, subtraction, negation, multiplication, division,
`Inverse(token)` and `Pow(signedExponent, token)`. Field methods expose the same
operations. Integer addition, subtraction and multiplication embed the integer
in the prime subfield; an integer is not interpreted as a base-p element encoding.
Zero division and negative powers of zero throw `DivideByZeroException`. The
convention for exponent zero is a^0=1, including 0^0. Multiplication is dense polynomial multiplication/reduction;
inversion uses a^(q-2). A default-initialized element is invalid and arithmetic
rejects it. Calling `Inverse` or `Pow` on a default element throws
`InvalidOperationException`; passing a default or incompatible element to a field
arithmetic method throws `ArgumentException`.

Fields compare by characteristic and normalized defining polynomial. Separately
constructed equal presentations can interoperate. Equal cardinality alone does
not permit mixing elements: automatic field isomorphisms and embeddings are not
part of this API. Equality and hashes include the presentation.

`Elements(maxElements, token)` lazily enumerates the whole field in increasing
base-p coefficient encoding, starting with zero. It rejects a field whose order
exceeds maxElements (default 1000000) with `ArithmeticException` before yielding
its first element; a negative limit throws `ArgumentOutOfRangeException`. As with
point enumeration, validation begins on the first `MoveNext`, not when the
enumerable is obtained.

## Curves and points over F_(p^k)

`EllipticCurveFq` takes a field and five general Weierstrass coefficients. One
constructor accepts field elements; another embeds integer coefficients into the
prime subfield. Singular equations and incompatible field presentations are
rejected. The general formulas apply in characteristics 2 and 3 as well.

```csharp
using System.Numerics;

var field = new FiniteField(3, new BigInteger[] { 1, 0, 1 });
var alpha = field.Generator;
var e = new EllipticCurveFq(field,
    field.Zero, field.Zero, field.Zero, alpha, field.One);
var p = e.CreatePoint(0, 1);
var twice = e.Double(p);
var inverse = e.Negate(p);
var multiple = e.Multiply(p, -7);
var count = e.CountPoints();
```

The public API exposes c4, c6, discriminant and j, membership checking, point
creation, negation, addition, subtraction, doubling and signed scalar multiplication.
`EllipticCurvePointFq.Infinity` is a universal infinity; a default point is invalid.
Both coordinates of affine points must belong to the curve's field presentation.
Group operations reject points outside the curve.

Point counting and enumeration use exact integer and finite-field arithmetic.
`CountPoints(maxWork, token)` includes infinity and counts the solutions of
y^2+b*y=rhs at each x without constructing the points. In odd characteristic,
the discriminant b^2+4*rhs gives one solution when zero, two when a nonzero square,
and none otherwise. Odd fields of order q<=2^20 use a square-root table (at most
4 MiB), with a cap independent of maxWork; larger fields use exact power tests.
Nonresidues for Tonelli-Shanks are sought in the whole field presentation. In
even degree, all nonzero prime-subfield elements are squares, so the search starts
at the generator without constructing or testing the p constants.

In characteristic two, b=0 gives one solution because squaring is bijective;
otherwise z=y/b satisfies z^2+z=rhs/b^2. The latter map has kernel {0,1}, so a
soluble equation gives exactly the two ordinates bz and bz+b. Squaring and
z->z^2+z are F_2-linear maps. Their exact echelon forms use O(k^2) bits for degree k
and solve the equations in both even and odd extension degrees. Solubility of the
Artin-Schreier equation is equivalent to vanishing absolute trace. Normalization
by the inverse square of b uses exact binary polynomial extended Euclid, followed
by squaring and reduction, avoiding a dense-polynomial exponentiation for every x.
For constant nonzero b, its inverse square is reused.
See the discussion of trace and quadratic equations in
[Fong et al., Field inversion and point halving revisited, section 4.2](https://cacr.uwaterloo.ca/techreports/2003/corr2003-18.pdf).

When all five coefficients belong to the prime subfield, counting first counts
over F_p and then uses the exact Frobenius trace recurrence to obtain the count
over F_(p^k). Writing t=p+1-#E(F_p), the recurrence is

    t_0=2, t_1=t, t_n=t*t_(n-1)-p*t_(n-2),
    #E(F_(p^k))=p^k+1-t_k.

This shortcut checks the coefficients, not just the j-invariant. Only completed
counts are cached on the curve instance; failed or cancelled attempts can be
retried. Pre-cancelled tokens throw even for cached results, and concurrent
waiters observe their own tokens.

The abscissa-based approach is described by PARI's
[`Flxq_ellcard_naive` in FlxqE.c](https://pari.math.u-bordeaux.fr/lcov-report/basemath/FlxqE.c.gcov.html).
PARI's `elltrace_extension` and `Fp_ffellcard` in
[FpE.c](https://pari.math.u-bordeaux.fr/lcov-report/basemath/FpE.c.gcov.html) use the
extension trace identity above. These are algorithmic references; the library
uses independently written C# and requires no PARI installation.

Both methods check q<=maxWork before starting, including before returning a
cached or prime-subfield count. The default limit is 1000000 x-coordinates.
`Points(maxWork, token)` visits each x and solves the quadratic in y, with
infinity yielded first. Odd fields beyond the square-table cap use exact
Tonelli-Shanks roots with bounded storage. The order of the two ordinates at a
fixed x is unspecified.
Both limits are checked using BigInteger to avoid overflow.
A work-limit exception does not report a partial count. Stopping lazy point
enumeration early does not establish completeness.
Both methods reject negative limits with `ArgumentOutOfRangeException` and
insufficient limits with `ArithmeticException`; `Points` checks its limit on the
first `MoveNext`, before yielding infinity.

Work limits count potential x-coordinates, not field operations or seconds.
One field operation can itself be expensive for large degree or characteristic.
Construction and basic arithmetic have no work counter. The `FiniteField`
constructor and `Multiply`, `Divide`, `Inverse` and `Pow` methods accept tokens.
`CreateElement`, field addition/subtraction/negation, the element operators and
`EllipticCurveFq` construction have no cancellation parameter. Curve `Multiply`,
`CountPoints` and `Points` propagate tokens into field multiplication, inversion
and powers. Cancellation is checked in those polynomial loops, table preparation,
linear solving, between point group operations, and at every enumerated x; a single
BigInteger operation or field addition is not internally interruptible. Large
fields remain usable for arithmetic even when counting or enumeration is disallowed. The existing
[`EllipticCurveFp` API](basic-extensions.md#curves-over-prime-fields) uses specialized
residue tables and cubic finite differences.

These APIs do not implement SEA, discrete logarithms, field embeddings,
extension-field isogenies, or curves over number fields.

## Verification

The tests use 46 independently generated 80-digit PARI Faltings-height rows, all
five existing LMFDB snapshots, precision up to 60 decimal places, rational model
changes, quadratic twists and semistable/additive examples. Reference decimal
comparisons have explicit tolerances and are not themselves interval proofs.

Finite-field tests compare the Rabin criterion with independent exhaustive trial
factorization of small monic polynomials. They check field laws, inverses,
Frobenius, immutable representations, linear extensions, large characteristics,
invalid inputs and limits. PARI fixtures contain 90 field-arithmetic rows and 35
curves over fields of orders 4,8,9,16,25,27,49,125, with genuine extension-field
coefficients, invariants, point counts and point operations. Base-change counts
are also checked against the Frobenius trace recurrence; degree-one group laws
are compared with the existing prime-field implementation, including nonzero
generators in degree-one presentations. All 768 nonsingular general models over
F_4 are checked against full coordinate-pair searches. Over F_(1031^2), with
t^2=-1, an independent prime-field Legendre test of Norm(a+b*t)=a^2+b^2 checks
both squares and nonsquares in the large-field fallback. Further cases cover
enumeration over F_(101^2), binary solving and inversion through degree 128,
cancellation and work limits. A curve over F_9 with j=0 and a nonconstant
coefficient guards against incorrectly treating a prime-subfield j-invariant
as proof that all coefficients lie in the prime subfield.

See [fixture regeneration instructions](../tests/Fixtures/README.md).

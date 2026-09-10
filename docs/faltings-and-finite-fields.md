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
Singular curves, invalid options, exhausted period work and cancellation follow
the existing period API's behavior. Neither height is constrained to be positive.

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
polynomials and composite characteristics are rejected.

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
Zero division and negative powers of zero throw. The convention for exponent zero
is a^0=1, including 0^0. Multiplication is dense polynomial multiplication/reduction;
inversion uses a^(q-2). A default-initialized element is invalid and arithmetic
rejects it.

Fields compare by characteristic and normalized defining polynomial. Separately
constructed equal presentations can interoperate. Equal cardinality alone does
not permit mixing elements: automatic field isomorphisms and embeddings are not
part of this API. Equality and hashes include the presentation.

`Elements(maxElements, token)` lazily enumerates the whole field in increasing
base-p coefficient encoding, starting with zero. It rejects a field whose order
exceeds maxElements (default 1000000) before yielding its first element.

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

`Points(maxWork, token)` visits every affine pair (x,y), with infinity yielded
first. `CountPoints(maxWork, token)` counts this complete enumeration. Both check
q^2<=maxWork before starting, using BigInteger to avoid overflow; the default is
1000000 coordinate-pair checks. Thus the default supports q<=1000. A work-limit
exception does not report a partial count. Stopping lazy enumeration early does
not establish completeness.

This intentionally simple O(q^2) search counts pairs, not field operations or
seconds. One field operation can itself be expensive for large degree or
characteristic. Construction and basic arithmetic have no work counter. Cancellation
is checked in polynomial power/reduction loops, between point group operations,
and at every enumerated pair; the operators and a single BigInteger operation
are not internally interruptible. Large fields remain usable for arithmetic even
when enumeration is disallowed. The existing `EllipticCurveFp` API retains its
more efficient direct counting for prime fields.

This addition does not implement SEA, discrete logarithms, field embeddings,
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
are compared with the existing prime-field implementation.

See [fixture regeneration instructions](../tests/Fixtures/README.md).

# Exact rational input and checked points

## Numeric input

```csharp
var fraction = BigRational.Parse("-7/12");
var decimalValue = BigRational.Parse("0.1");  // exactly 1/10
bool parsed = BigRational.TryParse("3e-2", out var value); // 3/100
var quotient = BigRational.Parse("1.5 / 0.25"); // 6
```

The grammar accepts signed ASCII integers, decimals with a dot, an optional
`e` or `E` exponent, and two such numbers separated by one `/`. Whitespace is
allowed around the input and each fraction component. Unicode minus `−` is
accepted as a sign, including in an exponent. Decimal commas, grouping separators,
NaN, infinity and expressions with operators are rejected. Conversion uses exact
integer arithmetic without a floating-point intermediate, independently of culture.

`Parse(null)` throws `ArgumentNullException`. Malformed text, a zero denominator
and resource-limit exhaustion throw `FormatException`. `TryParse` returns false
for these inputs, including null, and sets its output to `BigRational.Zero`.
Invalid parsing options throw `ArgumentOutOfRangeException` in either API,
before text validation or assignment to the `TryParse` output.

`RationalParsingOptions` can be supplied to `Parse(text, options)` or
`TryParse(text, options, out value)`. Null options select the defaults:

| Property | Default | Meaning |
| --- | --- | --- |
| `MaxTextLength` | 20,000 | Total input characters, including whitespace |
| `MaxExponent` | 4,096 | Absolute decimal exponent in each numeric component |
| `MaxValueBits` | 32,768 | Bit length of the final normalized numerator and denominator |
| `MaxIntermediateBits` | 131,072 | Bit length allowed while constructing components and the final fraction |

Length and bit limits must be positive; the exponent limit may be zero.
Intermediate limits are separate because large components can cancel to a small
result, such as `N/N`. Each component must still fit the intermediate budget;
for example, zero divided by an excessive component is rejected. The parser
checks length and exponents first, rejects obviously excessive digit strings and
decimal powers before allocation, and checks exact integer and product sizes
after construction. These are acceptance limits, not strict peak-memory or
elapsed-time bounds. Parsing has no cancellation-token overload.
The limits are snapshotted at the start of each call. A larger exact value can
be parsed with explicitly increased limits; round-tripping `value.ToString()`
requires enough space and bit allowance for that value.

The Explorer and console retain their existing input budgets through a thin
`RationalText` wrapper. Decimal display in the editor remains a presentation
choice; the library's canonical `ToString()` format remains an integer or `n/d`.

## Points and group operations

```csharp
var e = new EllipticCurveQ(0, 0, 1, -1, 0);
var p = e.CreatePoint(BigRational.Parse("1/4"), BigRational.Parse("-5/8"));
var q = e.DoubleChecked(p);
var sum = e.AddChecked(p, q);
var multiple = e.MultiplyChecked(p, -7, cancellationToken);
```

`CreatePoint(x,y)` checks nonsingularity and the exact curve equation. A singular
model gives `InvalidOperationException`; coordinates outside the curve give
`ArgumentException`. `TryCreatePoint(x,y,out point)` returns false for either
condition, placing an unused `Infinity` sentinel in the output. On success it
returns an affine point in the existing `EllipticCurvePoint` representation.

`AddChecked`, `DoubleChecked`, `NegateChecked`, `SubtractChecked` and
`MultiplyChecked` reject singular models and points outside their target curve.
Validation precedes identity shortcuts and a zero multiplier. The point at
infinity is accepted on any nonsingular model. `MultiplyChecked` honors
cancellation before validation and within scalar multiplication; the input is
validated once rather than on every internal group operation. Cancellation cannot
interrupt a membership check or an individual rational group operation.

Points are coordinate values without a stored curve identity. Membership is
checked against the receiving curve's equation. `default(EllipticCurvePoint)`
is the affine point `(0,0)`, whose validity depends on that equation. It is not
an invalid-point sentinel.

The raw point constructor and original `Add`, `Double`, `Negate`, `Subtract`
and `Multiply` retain their existing behavior and assume valid inputs. Singular
equations can still be represented and examined with `IsOnCurve`; checked group
operations reject them. Existing Explorer group actions already validate inputs,
so checked aliases do not create duplicate menu actions or change saved IDs.

## Arithmetic implementation

Multiplication and division cancel cross factors before forming products.
Addition and subtraction use the GCD of the denominators, then cancel remaining
common factors. Internal construction of a proven reduced result avoids a
second GCD, while the public fraction constructor continues to normalize arbitrary
inputs. Results keep a positive denominator, coprime numerator and denominator,
and the canonical zero `0/1`, including `default(BigRational)`.

The scalar multiplication loop skips its unused final doubling. Cross cancellation
can reduce temporary integer sizes when operands share factors; for small or
coprime inputs, the added GCD work may outweigh that benefit. Mathematical results
remain exact.

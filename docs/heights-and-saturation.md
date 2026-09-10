# Local data, heights, periods and saturation over Q

These native APIs run in managed C# using `BigInteger` and `BigRational`, without
network requests, subprocesses or mathematical-software dependencies. The source
is an independent implementation of mathematical formulas; PARI is used only to
generate reference test data. All methods require a nonsingular elliptic curve.

## Models and points

`ChangeModel(u,r,s,t)` uses

    x = u² x' + r,  y = u³ y' + s u² x' + t,  u != 0.

Its `WeierstrassIsomorphism` contains both models and `Map`, `MapBack`, `Inverse`.
`TryGetIsomorphism(target, out map)` verifies an exact rational isomorphism;
equal j-invariants alone do not suffice. `GetMinimalModelIsomorphism(token)` maps
to the reduced global integral minimal equation. A map need not select the same
sign of u as another isomorphism; point maps are always used consistently.

## Local data

`GetLocalData(token)` returns all bad primes, sorted increasingly.
`GetLocalData(prime, token)` also accepts good primes, returning I0, conductor
exponent 0, Tamagawa number 1 and root number +1. Composite arguments are rejected.
`LocalData` and `TamagawaProduct` are convenience properties.

`LocalReductionData` contains the prime, minimal discriminant and conductor
valuations, j-denominator valuation, Kodaira symbol, splitting type, Tamagawa number
and local root number. Tate's algorithm includes the wild branches at 2 and 3,
the I_n* refinement and exceptional starred types. Auxiliary cubic roots are
counted by polynomial gcd over F_p. See
[Cremona, Chapter III, §3.2](https://johncremona.github.io/book/fulltext/chapter3.pdf).

## Height conventions and guarantees

`NaiveHeight(P)` is log(max(|numerator(x)|, denominator(x))). `CanonicalHeight(P)`
uses the Cremona/LMFDB normalization lim 4^(-n) h_x(2^n P). Infinity has height 0.
The canonical value is invariant under the supplied rational model changes.

`LocalHeight(P,p)` and `ArchimedeanHeight(P)` refer to the reduced global minimal
model and are undefined at infinity. Their sum is the canonical height. For integral
minimal coefficients, finite contributions outside the bad primes and x-denominator
primes vanish. Local contributions can be negative.

Finite terms use exact valuations; the real term uses the two-chart Tate–Silverman
series with a rigorous truncation bound. The chart is changed only after an interval
inequality certifies its validity. See
[Cremona, Chapter III, §3.4](https://johncremona.github.io/book/fulltext/chapter3.pdf).

`HeightPairing(P,Q)` = (h(P+Q)-h(P)-h(Q))/2, so a diagonal entry is h(P).
`HeightPairingMatrix(points)` does not assume independence.
`Regulator(points)` computes the determinant, certifying positive definiteness by
interval elimination. Dependent inputs, or insufficient precision to certify
independence, throw `ArithmeticException`. The empty determinant is 1. The result
is the curve regulator only for a full basis modulo torsion; a subgroup of index m
in the full Mordell–Weil lattice has regulator m² times the curve regulator.

## Numerical controls

All successful numerical results are `RealEnclosure` objects with exact rational
endpoints and width at most 10^(-`DecimalDigits`). The default is 12 absolute decimal
places. `Approximation` converts the midpoint to double for display and is not
used in any proof. Arithmetic rounds outward on a binary grid; logarithms have
geometric tail bounds, square roots use integer bounds, and pi uses alternating
arctangent series with explicit tails.

`RealComputationOptions` defaults to 256 fractional binary bits, four extra decimal
places for height truncation, 512 height/AGM iterations and 100000 root-isolation
steps. Logarithm and pi series lengths depend on `PrecisionBits`. `GuardDigits`
can tighten the height-series tail for large regulators; increase `PrecisionBits`
alongside it, since longer iterations also amplify rounding uncertainty.
For example, `DecimalDigits=30, PrecisionBits=512` and
`DecimalDigits=60, PrecisionBits=1024` are covered by reference tests.

An ambiguous chart, denominator or positive pivot, or an unmet width/work limit,
causes an explicit exception. Requested accuracy is never silently downgraded.
Large coefficients or nearly dependent generators may require more precision.
Minimalization and factorization costs are outside numerical work counters.

## Period conventions

`GetPeriods(options, token)` computes periods of dx/(2y+a1*x+a3) on the reduced
global minimal model. `PrimitiveRealPeriod` is the least positive real period w1.
The second basis period has positive imaginary part; its real part is 0 for
positive discriminant and w1/2 for negative discriminant. This convention does
not assert an SL(2,Z)-reduced basis. `MinimalModel` records the differential's model.

`RealPeriod` is the BSD real period, w1 times the number of real components;
`Area` is w1*Im(w2). These are the quantities compared with LMFDB. Scaling an
input equation does not scale these outputs, since they use its minimal model.
The distinction between the real period and w1 follows the
[LMFDB period-lattice convention](https://www.lmfdb.org/knowledge/show/ec.q.period_lattice).

The cubic roots are enclosed by exact Sturm isolation, then positive real AGM
formulas produce the period enclosures. The implementation handles both signs
of the discriminant. The AGM enclosure follows the bounds of its arithmetic and
geometric iterates; see [DLMF §19.8](https://dlmf.nist.gov/19.8).

## Exact saturation at requested primes

`Saturate(points, primes, options, token)` works with the subgroup generated by
the supplied independent points **and all rational torsion**. Full rank is not
required. Independence is first certified by interval heights. Dependent or
uncertified input produces `IndependenceCertified=false` and unresolved primes.

For each requested prime p, the algorithm examines projective representatives of
nonzero coefficient vectors in F_p^r and all rational torsion translates. Small
good reductions can disprove divisibility. Remaining candidates are tested using
the exact multiplication-by-p equation, rational-root isolation, rational square
tests and the exact group law. Every successful enlargement replaces a generator
and has proved index p; the search is restarted for the enlarged subgroup.

To see completeness at a fixed p, any missing p-division yields a nonzero relation
in the input lattice modulo p. Normalize its first nonzero coefficient to 1;
subtracting integral multiples of p changes the divided point by known generators.
It is therefore represented in the finite enumeration, including its torsion
translate. Rational roots are complete by the rational-root theorem: after clearing
denominators and content, the leading coefficient times every rational root is
an integer, which exact Sturm intervals enumerate.

The result records `Generators`, `IndexGain`, `CertifiedPrimes`, `UnresolvedPrimes`,
`Work` and `Reason`. All returned generators are on the input model. Previously
certified primes stay certified after enlargement at another prime, whose index
is coprime to them. A limit preserves verified enlargements and existing certificates.

`IsComplete` means all **requested** primes were certified, including the vacuous
empty-prime request when input independence is certified. It does not certify
unrequested primes, full rank, or a full Mordell–Weil basis. In particular, the
library does not yet prove a global bound on all primes dividing the saturation
index. LMFDB generator data can be used for comparison, but is not such a proof.

Defaults are 2000000 counted operations, division degree at most 1024, and 128
enlargements. The independence calculation defaults to 16 digits at 384 bits.
Preparation (minimalization, torsion, heights) is outside the work counter.
Polynomial coefficient sizes and exact root-isolation costs can grow rapidly;
the method is intended for modest rank and small saturation primes. Counters are
not wall-clock or memory limits. Cancellation is cooperative; the existing torsion
enumeration is not interruptible internally, and cancellation is observed when
that preparation completes.

## LMFDB adapter

`LmfdbEllipticCurve.FetchAsync(curve, httpClient, cancellationToken)` first queries
`ec_curvedata` by exact rational j, paginating candidates in groups of 100 up to
the server's offset limit. Each candidate is checked for rational isomorphism.
It then fetches the aggregate curve-data JSON snapshot. A caller-owned `HttpClient`
is never disposed. The synchronous constructor uses the same logic; all properties
read cached data and issue no requests. There are no automatic retries or native
fallbacks. HTTP failures, CAPTCHA, malformed/truncated used tables and mismatched
records are reported explicitly.

| Native computation | Stored LMFDB field |
|---|---|
| `GetLocalData` | `ec_localdata` valuations, Kodaira, splitting, c_p, w_p |
| `TamagawaProduct` | `ec_mwbsd.tamagawa_product` |
| `CanonicalHeight(P)` for each stored generator | `ec_mwbsd.heights` |
| `Regulator(storedGenerators)` | `ec_curvedata.regulator` |
| `GetPeriods().RealPeriod` | `ec_mwbsd.real_period` |
| `GetPeriods().Area` | `ec_mwbsd.area` |
| `FaltingsHeight()` | `ec_curvedata.faltings_height` |
| `StableFaltingsHeight()` | `ec_curvedata.stable_faltings_height` |
| Generators as saturation inputs/reference | `ec_mwbsd.gens`, `torsion_generators` |
| Proven rank bounds as a cross-check | `ec_mwbsd.rank_bounds` |

`MordellWeilGenerators` and `TorsionGenerators` are on `GlobalMinimalModel`.
Stored triples [a,b,c] use weighted coordinates x=a/c², y=b/c³; each parsed point
is checked against the model. `GetGeneratorsOnModel(curve)` performs the exact
coordinate change. Missing optional fields are null; present empty lists remain
empty, which matters for rank-zero curves. `FromStoredDataJson(json)` supports
offline aggregate snapshots such as
[the 37.a1 data page](https://www.lmfdb.org/EllipticCurve/Q/data/37.a1).

`LmfdbRealValue` retains decimal text and available precision metadata. Its
`AsRational()` exactly represents the supplied decimal, not the unknown true real
number. LMFDB does not supply guaranteed error bounds for these approximations;
see [its reliability statement](https://www.lmfdb.org/EllipticCurve/Q/Reliability).
The adapter supplies stored values, not remote methods for heights of arbitrary
points, arbitrary regulators, complex period bases or subgroup saturation.

Native Faltings-height formulas and their minimal-model normalization are described
in [the Faltings-height and finite-extension notes](faltings-and-finite-fields.md).

## Verification

Offline tests compare exact local invariants against PARI, and numerical enclosures
against 70-digit PARI reference values. Stored LMFDB snapshots cover ranks 0–3,
torsion, both period shapes and generators/regulators. Saturation tests deliberately
multiply and mix known generators, include torsion translates and nonminimal models,
then verify indices, subgroup regulators and honest partial results. Independent
exact group-law checks also exercise multiplication polynomials through degree 169.
See [fixture generation](../tests/Fixtures/README.md).

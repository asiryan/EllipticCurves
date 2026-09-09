# Unconditional native two-descent

`GetRankBounds()` selects descent by 2-isogeny when rational 2-torsion is present,
and general binary-quartic descent otherwise. Both run entirely in C#.
`GetRankBounds(RankComputationOptions, CancellationToken)` supplies work limits;
`PreferGeneralTwoDescent` allows the general construction on curves with 2-torsion.

The mathematical reference for the covering construction and reduction regions is
[Cremona, Algorithms for Modular Elliptic Curves, III.3.6, pp. 79–97](https://johncremona.github.io/book/fulltext/chapter3.pdf).
Equivalence uses the corrected criterion in
[Cremona–Fisher, On the equivalence of binary quartics, Lemma 10 and Theorem 12](https://johncremona.github.io/papers/quartequiv.pdf).
The implementation is written for this library from the mathematical formulas;
it does not include or translate eclib/PARI implementation code. PARI is used only
to produce independent test data.

## General descent: conditions needed for an upper bound

1. Reduce the input to its global minimal integral model. Set I=c4 and J=2*c6.
   If 16 divides I and 64 divides J, divide by these factors once. Search the
   resulting pair, and also (16I,64J) unless 4|I, 8|J and 16|(2I+J).
2. Isolate the real roots of z^3-3Iz+J with exact Sturm sequences. The three
   reduction types on p. 94 give finite bounds on a,b,c. Rational interval
   arithmetic encloses those bounds and rounds the integer endpoints outwards.
   Include overlapping regions and all boundary representatives. No floating-point
   root or comparison can exclude a covering. Insufficient root separation produces
   an incomplete result, never an asserted upper bound.
3. With p=3b^2-8ac, require p^3-48Ia^2p-64Ja^3=27r^2 exactly. Recover integer
   d and e and recheck both quartic invariants. Both signs of r are included.
4. Normalize the two invariant levels to the larger level by multiplying the smaller
   form by 4. A rational root (including a root at infinity) represents the identity
   class. Remove it and count that class exactly once.
5. Identify duplicate classes by the corrected Cremona–Fisher criterion. Choose a
   rational projective point where the first form's sextic covariant is nonzero.
   If its quartic and Hessian values there are a and p, equivalence holds precisely
   when a*Hessian(second)-p*second has a rational projective root. Seven affine
   integer points suffice to find a valid specialization when infinity is unsuitable.
   The older product-seminvariant test in the book can fail for reducible resolvents
   and is not used. Root-count signatures at good primes only reject duplicates;
   actual equivalence requires the exact rational-root test. Sturm sequences and
   the rational-root theorem give a finite check without factoring large constants.
6. Retain exactly the everywhere locally soluble classes. Their total number,
   including the identity, must be a power of two. Its base-2 logarithm is the
   dimension of Sel_2(E/Q). Subtract dim E(Q)[2] for the rank upper bound.

Stopping before all regions and both invariant levels have been searched makes
the upper bound unknown. A partial class count is never used as an upper bound.
The power-of-two check is a consistency check, not a substitute for completeness.

## Exact local solubility

For a nonsingular binary quartic F(u,v), two charts cover all Q_p-points:
v=1 with u in Z_p, and u=1 with v in p Z_p. Infinity is included in the latter.
An iterative depth-first search subdivides each chart into p-adic balls.

At an integer center x, a p-adic square F(x,1) proves solubility, including zero.
If v_p(F(x,1)) > 2 v_p(F'(x,1)), Hensel's lemma gives a root congruent to x
modulo p, also proving solubility in the original chart.

For a ball x+p^n Z_p, expand F(x+t,1) exactly. The minimum valuation of its
nonconstant Taylor terms at t=p^n gives a lower bound on every change in F.
When this fixes a nonsquare square class throughout the ball, reject it. Odd
valuations are nonsquares; unit squares are characterized modulo p for odd p,
and modulo 8 for p=2. Otherwise subdivide. No surviving finite congruence test
is accepted as a proof of solubility.

Nonsingularity ensures termination in unlimited arithmetic: around a nonzero
value the square class is locally constant, and around a root Hensel eventually
applies. Compactness of the two Z_p charts yields a finite cover by decided balls.
The work budget may stop earlier, in which case descent has no certified upper bound.

Over R, a positive leading coefficient or a real root proves solubility; a negative
definite form is insoluble. Odd primes not dividing the quartic discriminant need
no check: the smooth genus-one curve over F_p has a point, which lifts. General
descent uses 2 and the prime divisors of the minimal elliptic discriminant. Isogeny
descent uses 2 and the prime divisors of b(a^2-4b).

## Lower bounds from rational points

On the integral cubic obtained by X=4x, Y=8y+4a1*x+4a3, every simple root r
modulo a good odd prime p gives a Kummer character to F_2. Evaluate the Legendre
symbol of X-r; at (r,0) use the derivative of the cubic, and at infinity use 1.
These are homomorphisms. Combine their bits for all selected primes and use
exact binary elimination on the images of rational points.

If the image dimension is d and dim E(Q)[2]=t, then d<=rank(E(Q))+t,
so max(0,d-t) is an unconditional lower bound. This accounts for even torsion,
including torsion points of order 4 or 8. Repeated points and dependent multiples
cannot increase this bound beyond the actual rank. The finite list of primes can
miss independence; this only weakens the lower bound. Canonical-height arithmetic
is not required by this certificate.

An additional exact infinite-order test preserves a lower bound of one when a
point is invisible in the chosen quotients, for example when it is twice another
point. A rational torsion point over Q has order at most 12, so checking its first
12 multiples is sufficient.

General descent also counts inequivalent coverings with witnessed rational points.
With m such classes including the identity, ceil(log2(m))-t is a lower bound,
because the image E(Q)/2E(Q) has power-of-two order. The classical covariant map
transfers each witness to the minimal model and its curve equation is checked
exactly before accepting it. Isogeny descent instead spans witnessed square classes.

## Limits and result semantics

Defaults are SearchBound=32, MaxSquareClasses=65536, MaxDescentWork=5000000,
MaxPointSearchWork=1000000 and ReductionPrimeBound=101. SearchBound controls both
the numerator/denominator box for x on the minimal curve and primitive quartic
coordinates. Zero disables all point searches.

Descent work counts region enumeration, root-isolation/equivalence work and local
lifting steps. Point work is shared by the original curve and covering searches.
Neither counter measures wall time or bounds a single arbitrary-precision operation.
Initial minimization, factorization and rational 2-torsion detection are cancellable
but outside these counters. Very large invariants can require impractical work.

`TwoSelmerDimension` is populated only after complete general descent.
`UsedGeneralTwoDescent` and `UsedTwoIsogenyDescent` describe the selected method,
including incomplete attempts. `Reason` explains work exhaustion or a remaining
gap. `ExactRank` requires equal proved bounds. No numerical analytic result modifies
these bounds. No BSD, GRH, parity or Sha-finiteness assumption is used.

Full Selmer computation does not guarantee an exact rank: a nontrivial Sha[2] or
unfound rational points can leave a gap. Cassels-Tate pairings, higher descents,
height/regulator computation and a saturated Mordell-Weil basis are outside this API.

## Independent validation

`tests/Fixtures/generate-general-descent.gp` uses explicit coefficients and PARI's
`ellrank` and `ell2cover`; it performs no database lookup. The fixture records
algebraic rank bounds and the dimension of the 2-Selmer group, so testing does not
merely check agreement with a numerical analytic rank. See the
[fixture instructions](../tests/Fixtures/README.md) for regeneration.

The automated comparison covers 148 small reference curves, including full and
partial rational 2-torsion. Separate tests prove ranks 0 through 4, preserve a
rank gap on examples with nontrivial Sha[2], and verify nonminimal rational
coordinate changes. The Cremona–Fisher counterexample and zero-sextic specializations
are regression tests. Local tests include exhaustive modular obstructions, deep
p-adic coordinate changes and the everywhere locally soluble Reichardt–Lind covering.

Run the complete suite with `dotnet test EllipticCurves.sln -c Release`.

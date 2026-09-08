# Native arithmetic over Q

The public entry points are `EllipticCurveQ.GlobalMinimalModel`, `Conductor`,
`GetGlobalMinimalModel(CancellationToken)`, `GetConductor(CancellationToken)` and
`GetRankBounds(int searchBound, int maxSquareClasses, CancellationToken)`.
They perform no HTTP requests, start no processes and use no elliptic-curve database.
Singular input is rejected. Computations use `BigInteger` and `BigRational` throughout.

## Minimal model and conductor

Clear coefficient denominators using a rational change of variables. For each prime
dividing the integral discriminant, try dividing the invariants by p^4 and p^6.
A candidate is accepted only when an integral Weierstrass equation exists with those
invariants. This is checked by reconstructing the 12 possible reduced coefficient
patterns a1,a3 in {0,1}, a2 in {-1,0,1}. Repeat until no further division is possible.
The result is the reduced global minimal model.

On that model, the local conductor exponent is zero at good primes and one at
multiplicative primes. At additive primes p >= 5 it is two. At 2 and 3 we use Tate's
successive coordinate transformations, including the I_n* refinement loop and
the II*, III*, IV* branches. The conductor is the product of p raised to these exponents.

The mathematical reference is [Cremona, Algorithms for Modular Elliptic Curves,
Chapter III, sections 3.1–3.2](https://johncremona.github.io/book/fulltext/chapter3.pdf).
The implementation uses integral reconstruction for minimization and applies only
the required conductor portion of Tate's algorithm; it does not expose Tamagawa numbers.

## Rank bounds with rational 2-torsion

Completing the square and finding an integral root of the resulting monic cubic
gives an isomorphic integral equation

    E:  y² = x³ + a x² + b x.
    E': y² = x³ - 2a x² + (a² - 4b) x.

For each of these curves, consider the image of the descent homomorphism in
Q*/Q*²: infinity maps to 1, (0,0) maps to b, and other points map to the square
class of x. Its possible elements are the signed squarefree divisors d of b.
Such a class is in the image exactly when

    N² = d U⁴ + a U² V² + (b/d) V⁴

has a primitive integer solution (U,V), with U,V not both zero.
The two image dimensions satisfy

    rank(E(Q)) = dim(image on E) + dim(image on E') - 2.

The implementation obtains a lower bound for each image dimension by exact
quartic point searches and Gaussian elimination on square classes. It starts with
the known class of (0,0). The searches use 0 <= U,V <= searchBound and gcd(U,V)=1;
signs of U,V are immaterial because only even powers occur.

For an upper bound it eliminates classes with a real obstruction or no primitive
solution modulo one of 256, 81, 25, 49, 11, 13, 17, 19, 23, 29, 31. Both projective
charts are checked, including points for which a denominator is divisible by p.
If s classes survive, the actual image dimension is at most floor(log2(s)). This
uses containment and the power-of-two size of the actual image; the finite sieve
survivors themselves need not form a group.

Passing the sieve is **not** a claim of solubility over every Q_p, or over Q.
The sieve gives a valid but potentially weaker bound than full isogeny Selmer
groups. Equal final bounds prove the exact algebraic rank. The method uses neither
an analytic rank estimate nor BSD, GRH, or a parity conjecture.

For the descent construction and its distinction between local and global
solubility, see [Cremona, Chapter III, section 3.6](https://johncremona.github.io/book/fulltext/chapter3.pdf).

## Scope and computation limits

Without rational 2-torsion, the current fallback checks a bounded rational-point
box on the reduced minimal model. A found point is proved to have infinite order
if none of its first 12 multiples is infinity, using Mazur's torsion theorem.
This certifies a lower bound of one; otherwise the lower bound is zero. The upper
bound remains unknown (`null`). No general 2-descent or higher descent is implemented.
Found points are not presented as a Mordell-Weil basis.

`searchBound = 0` disables point searches. Increasing it can improve the lower
bound, but cannot guarantee equality of the bounds. `maxSquareClasses` bounds the
number of candidates per isogeny and throws `NotSupportedException` if exceeded.
It does not cap factorization time or the work of finding a cubic root.

All new factorizations certify their prime factors: deterministic Miller–Rabin
below 2^64, and recursive full n-1 primality proofs above it, with an exact trial
division fallback. Pollard rho supplies candidate factors. Unlike the older
torsion helper, no probable-prime result is accepted as a proof. Large integers
can still be prohibitively expensive to factor or prove prime. Cancellation
is checked in the factorization and search loops; a single BigInteger operation
cannot be interrupted midway.

The convenience properties recompute their results; callers doing repeated work
can retain the returned minimal model, conductor and rank bounds.

## Regression data

The tests include exact examples of ranks 0, 1 and 2, rational and non-minimal
coordinate changes, unknown results, cancellation and a strong pseudoprime which
passed the old helper's Miller–Rabin bases. Additional fixed fixtures were computed
independently using local PARI/GP, without curve labels or database lookup.
See [fixture generation instructions](../tests/Fixtures/README.md).
PARI/GP is not a build, test-run or runtime dependency; it is needed only to
regenerate those fixtures.

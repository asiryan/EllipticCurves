# Native arithmetic over Q

The public entry points are `EllipticCurveQ.GlobalMinimalModel`, `Conductor`,
`GetGlobalMinimalModel(CancellationToken)`, `GetConductor(CancellationToken)` and
`GetRankBounds(int searchBound, int maxSquareClasses, CancellationToken)`.
They perform no HTTP requests, start no processes and use no elliptic-curve database.
Singular input is rejected. The algebraic computations use `BigInteger` and
`BigRational` throughout. The additional analytic API uses numerical integration
and a separate rigorous interval calculation, as described below.

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

## Native analytic rank

`EstimateAnalyticRank(AnalyticRankOptions, CancellationToken)` is independent of the
descent API. `RootNumber` / `GetRootNumber(CancellationToken)` expose the exact sign
W of the functional equation, including the real factor -1.

Local signs use the reduction type for p >= 5 and invariant congruences for 2 and 3.
References: [Rizzo, Tables II and III](https://doi.org/10.1023/A:1022669121502) and
[Cowland Kellock–Dokchitser, Appendix A](https://arxiv.org/abs/2303.07883).
For normalized valuations (v(c4),v(c6),v(Delta)) = (2,5,0), the code retains
Rizzo's Table III condition for b >= 5; the 2023 table reverses that particular
row. The regression fixtures include examples distinguishing these conditions,
computed independently by PARI/GP. At 3 the c4 terms in Table II are interpreted
with their indicated normalized valuations, including the (>=4,6,9) branch.

On the minimal model, count points on its cubic modulo each prime p up to M and
set a_p = p+1-#E(F_p), including the point at infinity. At bad primes, counting
the singular cubic gives a_p = 0 or +/-1. A quadratic-residue sieve makes each
prime's point count O(p). Generate a_n multiplicatively with
a_(p^k) = a_p a_(p^(k-1)) - p a_(p^(k-2)) at good primes, and a_(p^k)=a_p^k
at bad primes. Coefficients are exact integers.

Put alpha=2*pi/sqrt(N), Lambda(s)=alpha^(-s) Gamma(s)L(E,s), and
F(z)=alpha Lambda(1+z). Splitting the Mellin integral at t=1 gives

    F^(k)(0) = alpha (1 + W*(-1)^k)
               sum_(n>=1) a_n integral_1^infinity exp(-alpha*n*t) log(t)^k dt.

After t=exp(u), integrate all moments on a finite interval with adaptive Simpson
quadrature. Repeat with different initial subdivisions and tighter tolerance.
The error diagnostics combine the quadrature estimates, disagreement between runs,
an allowance for rounding, and analytic truncation bounds using |a_n| <= 2n and
log(t)^k <= k! t. These diagnostics are not rigorous floating-point enclosures.
Convert moments to ordinary L derivatives using the Taylor expansion of
alpha^z/Gamma(1+z). `Derivatives[k]` is L^(k)(E,1), with no division by k!.

Rank detection uses F, whose parity is exact; ordinary L derivatives do not have
this parity. A derivative whose magnitude plus estimated error is below
`ZeroTolerance` is treated as a numerical zero. The first derivative of the
allowed parity exceeding 10 times that threshold, with error below a quarter
of the threshold, determines `EstimatedRank`. Ambiguity or failure to find a
nonzero derivative produces `Inconclusive`, not a claimed lower rank bound.

### Rigorous low-rank certificates

For candidate ranks 0 and 1, recompute respectively

    L(E,1)  = 2 sum_(n>=1) (a_n/n) exp(-alpha*n),          W=+1,
    L'(E,1) = 2 sum_(n>=1) (a_n/n) E1(alpha*n),          W=-1,
    E1(x)   = integral_x^infinity exp(-t)/t dt.

Every interval endpoint is an integer times 2^(-160); each operation rounds
outwards using integer division. Pi is enclosed by Machin's alternating arctangent
series, sqrt(N) by integer square roots, and exp(-x) by a positive exponential
Taylor series, inversion and repeated squaring. No `Math` transcendental call or
floating-point comparison enters certification.

To enclose E1(x), subdivide [x,96] into intervals [l,r] with r <= 2l.
For c=(l+r)/2, h=(r-l)/2, q=h/c <= 1/3 and S_j(c)=sum_(k=0)^j c^k/k!,

    integral_l^r exp(-t)/t dt
      = (2h exp(-c)/c) sum_(j even >=0) S_j(c) q^j/(j+1).

The terms are positive. Keep j=0,2,...,78 and bound the remainder above by
2q*q^80/[81(1-q^2)], since exp(-c)S_j(c) <= 1. The remaining integral from
96 to infinity is at most exp(-96)/96. Monotonicity of E1 encloses uncertain
lower integration limits. For an argument >=96, use [0,exp(-x)/x].

With M certificate terms and q=exp(-alpha), the omitted Fourier tails are bounded by

    4q^(M+1)/(1-q)                  for rank 0,
    4q^(M+1)/(alpha*(M+1)*(1-q))    for rank 1.

Only an interval excluding zero sets `Status=Certified` and `ProvenRank`.
For rank 1, the exact negative root number already proves L(E,1)=0.
The equality of algebraic and analytic rank in these two cases is unconditional
by modularity and the Gross–Zagier/Kolyvagin rank theorems; see
[Wiles's BSD survey, pp. 3–4](https://www.claymath.org/wp-content/uploads/2022/05/birchswin.pdf).
For higher ranks, small lower derivatives are not certified zeros. Even assuming
BSD does not eliminate this numerical uncertainty. The analytic method never
silently modifies the unconditional descent bounds.

### Analytic computation limits

Defaults: derivative order 4 (supported range 0–8), `ZeroTolerance=1e-9`,
`MaxTerms=20000`, `MaxPointCountingWork=20000000` (sum of primes),
`MaxIntegrationEvaluations=100000`, and `MaxCertificationTerms=1024`.
The numerical stage uses double precision. Certification uses at most the already
computed number of coefficients; increasing its limit cannot change a numerical
estimate into a proof unless the resulting interval excludes zero.

Coefficient demand grows approximately as sqrt(N) times the logarithm of the
requested accuracy. Direct point counting costs O(M^2/log M); the coefficient
and residue arrays use O(M) storage. Integration evaluates the finite Fourier
polynomial in O(M) time per sample. The method is intended for modest conductors,
not cryptographic-size curves. Work limits return `Inconclusive`; invalid options
and singular curves throw. Cancellation propagates as `OperationCanceledException`.
The earlier minimization and factorization steps can dominate the cost and are
not bounded by these numerical work limits. This is not a complete BSD leading-term
formula implementation: periods, regulators, Tamagawa factors and Sha are absent.

## Regression data

The tests include exact examples of ranks 0, 1 and 2, rational and non-minimal
coordinate changes, unknown results, cancellation and a strong pseudoprime which
passed the old helper's Miller–Rabin bases. Additional fixed fixtures were computed
independently using local PARI/GP, without curve labels or database lookup.
See [fixture generation instructions](../tests/Fixtures/README.md).
PARI/GP is not a build, test-run or runtime dependency; it is needed only to
regenerate those fixtures.

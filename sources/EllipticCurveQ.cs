using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace EllipticCurves
{
    /// <summary>
    /// Elliptic curve over ℚ in the general Weierstrass form
    ///     y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6 .
    /// All coefficients are exact rationals (BigRational). This type implements:
    ///  • exact invariants (b2,b4,b6,b8), c4, c6, Δ, j,
    ///  • the general Weierstrass group law (Add/Double/Multiply),
    ///  • conversion to the short model y^2 = x^3 + A x + B in characteristic 0,
    ///  • bounded search for rational/integral points,
    ///  • a self-contained computation of torsion points over ℚ
    ///    (via Lutz–Nagell + reductions), including mapping back to the original model.
    /// </summary>
    /// <remarks>Create an elliptic curve y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6.</remarks>
    public sealed partial class EllipticCurveQ(BigRational a1, BigRational a2, BigRational a3, BigRational a4, BigRational a6) : IEquatable<EllipticCurveQ>
    {
        #region Coefficients

        /// <summary>a1 coefficient in the general Weierstrass equation.</summary>
        public BigRational A1 { get; } = a1;
        /// <summary>a2 coefficient in the general Weierstrass equation.</summary>
        public BigRational A2 { get; } = a2;
        /// <summary>a3 coefficient in the general Weierstrass equation.</summary>
        public BigRational A3 { get; } = a3;
        /// <summary>a4 coefficient in the general Weierstrass equation.</summary>
        public BigRational A4 { get; } = a4;
        /// <summary>a6 coefficient in the general Weierstrass equation.</summary>
        public BigRational A6 { get; } = a6;

        #endregion

        #region Invariants (b2, b4, b6, b8), c4, c6, discriminant, j-invariant

        /// <summary>b2 = a1^2 + 4 a2.</summary>
        public BigRational B2 => A1 * A1 + BigRational.FromInt(4) * A2;

        /// <summary>b4 = 2 a4 + a1 a3.</summary>
        public BigRational B4 => BigRational.FromInt(2) * A4 + A1 * A3;

        /// <summary>b6 = a3^2 + 4 a6.</summary>
        public BigRational B6 => A3 * A3 + BigRational.FromInt(4) * A6;

        /// <summary>b8 = a1^2 a6 + 4 a2 a6 − a1 a3 a4 + a2 a3^2 − a4^2.</summary>
        public BigRational B8 => A1 * A1 * A6 + BigRational.FromInt(4) * A2 * A6 - A1 * A3 * A4 + A2 * A3 * A3 - A4 * A4;

        /// <summary>c4 = b2^2 − 24 b4.</summary>
        public BigRational C4 => B2 * B2 - BigRational.FromInt(24) * B4;

        /// <summary>c6 = −b2^3 + 36 b2 b4 − 216 b6.</summary>
        public BigRational C6 => -B2 * B2 * B2 + BigRational.FromInt(36) * B2 * B4 - BigRational.FromInt(216) * B6;

        /// <summary>
        /// Discriminant Δ = −b2^2 b8 − 8 b4^3 − 27 b6^2 + 9 b2 b4 b6.
        /// The curve is nonsingular iff Δ ≠ 0.
        /// </summary>
        public BigRational Discriminant
        {
            get
            {
                var t1 = -(B2 * B2) * B8;
                var t2 = -BigRational.FromInt(8) * B4 * B4 * B4;
                var t3 = -BigRational.FromInt(27) * B6 * B6;
                var t4 = BigRational.FromInt(9) * B2 * B4 * B6;
                return t1 + t2 + t3 + t4;
            }
        }

        /// <summary>j-invariant j = c4^3 / Δ (only defined if Δ ≠ 0).</summary>
        public BigRational JInvariant
        {
            get
            {
                var delta = Discriminant;
                if (delta.IsZero) throw new InvalidOperationException("Singular curve: discriminant = 0.");
                var c4val = C4;
                return c4val * c4val * c4val / delta;
            }
        }

        /// <summary>True iff the curve is singular (Δ = 0).</summary>
        public bool IsSingular => Discriminant.IsZero;

        #endregion

        #region Curve membership and group law (general Weierstrass formulas)

        /// <summary>Exact membership test for a point in affine coordinates.</summary>
        public bool IsOnCurve(EllipticCurvePoint P)
        {
            if (P.IsInfinity) return true;
            var x = P.X; var y = P.Y;
            var lhs = y * y + A1 * x * y + A3 * y;
            var rhs = x * x * x + A2 * x * x + A4 * x + A6;
            return lhs == rhs;
        }

        /// <summary>Group inverse: −(x,y) = (x, −y − a1 x − a3).</summary>
        public EllipticCurvePoint Negate(EllipticCurvePoint P)
        {
            if (P.IsInfinity) return P;
            return new EllipticCurvePoint(P.X, -(P.Y + A1 * P.X + A3));
        }

        /// <summary>
        /// Group law (general Weierstrass). Handles P, Q, doubling, and the vertical-tangent case.
        /// Slope:
        ///  • If x1 ≠ x2: λ = (y2 − y1)/(x2 − x1).
        ///  • If P = Q:   λ = (3x1^2 + 2a2 x1 + a4 − a1 y1) / (2y1 + a1 x1 + a3).
        /// Then:
        ///  x3 = λ^2 + a1 λ − a2 − x1 − x2,
        ///  ν = y1 − λ x1,
        ///  y3 = −(λ + a1) x3 − a3 − ν.
        /// </summary>
        public EllipticCurvePoint Add(EllipticCurvePoint P, EllipticCurvePoint Q)
        {
            if (P.IsInfinity) return Q;
            if (Q.IsInfinity) return P;

            // P == −Q ?
            if (P.X == Q.X && P.Y + Q.Y + A1 * Q.X + A3 == BigRational.Zero)
                return EllipticCurvePoint.Infinity;

            BigRational lambda;
            if (P.X != Q.X)
            {
                lambda = (Q.Y - P.Y) / (Q.X - P.X);
            }
            else
            {
                var num = BigRational.FromInt(3) * P.X * P.X
                        + BigRational.FromInt(2) * A2 * P.X
                        + A4 - A1 * P.Y;
                var den = BigRational.FromInt(2) * P.Y + A1 * P.X + A3;
                if (den.IsZero) return EllipticCurvePoint.Infinity; // vertical tangent
                lambda = num / den;
            }

            var x3 = lambda * lambda + A1 * lambda - A2 - P.X - Q.X;
            var nu = P.Y - lambda * P.X;
            var y3 = -(lambda + A1) * x3 - A3 - nu;
            return new EllipticCurvePoint(x3, y3);
        }

        /// <summary>Point doubling: 2P = P + P.</summary>
        public EllipticCurvePoint Double(EllipticCurvePoint P) => Add(P, P);

        /// <summary>Point subtraction: P − Q = P + (−Q).</summary>
        public EllipticCurvePoint Subtract(EllipticCurvePoint P, EllipticCurvePoint Q) => Add(P, Negate(Q));

        /// <summary>
        /// Scalar multiplication nP using left-to-right double-and-add.
        /// Supports n &lt; 0 (via negation) and n = 0 (returns O).
        /// </summary>
        public EllipticCurvePoint Multiply(EllipticCurvePoint P, BigInteger n)
        {
            if (n.Sign == 0) return EllipticCurvePoint.Infinity;
            if (n.Sign < 0) return Multiply(Negate(P), BigInteger.Negate(n));
            EllipticCurvePoint acc = EllipticCurvePoint.Infinity;
            EllipticCurvePoint basePt = P;
            var k = n;
            while (k > BigInteger.Zero)
            {
                if (!k.IsEven) acc = Add(acc, basePt);
                basePt = Double(basePt);
                k >>= 1;
            }
            return acc;
        }

        #endregion

        #region Short Weierstrass reduction: y^2 = x^3 + A x + B

        /// <summary>
        /// Convert to the short Weierstrass model in characteristic 0.
        /// Steps:
        ///  1) Complete the square: y = y' − (a1 x + a3)/2 ⇒ y'^2 = x^3 + a2' x^2 + a4' x + a6'
        ///  2) Remove x^2-term:    x = X − a2'/3          ⇒ y'^2 = X^3 + A X + B
        /// Returns a curve with (a1,a2,a3) = (0,0,0) and (a4,a6) = (A,B).
        /// </summary>
        public EllipticCurveQ ShortWeierstrass
        {
            get
            {
                var s1 = A1 / BigRational.FromInt(2);
                var s3 = A3 / BigRational.FromInt(2);
                var a2p = A2 + s1 * s1;                              // a2'
                var a4p = A4 + BigRational.FromInt(2) * s1 * s3;     // a4'
                var a6p = A6 + s3 * s3;                              // a6'

                var A = a4p - a2p * a2p / BigRational.FromInt(3);
                var B = a6p - a2p * a4p / BigRational.FromInt(3)
                             + BigRational.FromInt(2) * a2p * a2p * a2p / BigRational.FromInt(27);

                return new EllipticCurveQ(
                    BigRational.Zero,
                    BigRational.Zero,
                    BigRational.Zero,
                    A,
                    B
                );
            }
        }

        /// <summary>
        /// Enumerate rational points with bounded x = m/n (|m| ≤ numMax, 1 ≤ n ≤ denMax, gcd(m,n)=1).
        /// Uses the substitution y' = y + (a1 x + a3)/2 to test y'^2 being a rational square.
        /// Returns O first, then affine points found in the box.
        /// </summary>
        public IEnumerable<EllipticCurvePoint> RationalPoints(int numMax, int denMax)
        {
            yield return EllipticCurvePoint.Infinity;

            var two = BigRational.FromInt(2);

            for (BigInteger n = 1; n <= denMax; n++)
            {
                for (BigInteger m = -numMax; m <= numMax; m++)
                {
                    if (BigInteger.GreatestCommonDivisor(BigInteger.Abs(m), n) != BigInteger.One) continue;

                    var x = new BigRational(m, n);

                    // S(x) = x^3 + a2 x^2 + a4 x + a6 + ((a1 x + a3)^2)/4  must be a rational square
                    var t = A1 * x + A3;
                    var rhs = x * x * x + A2 * x * x + A4 * x + A6 + t * t / BigRational.FromInt(4);

                    if (BigRational.IsSquare(rhs, out var yp))
                    {
                        var y1 = yp - t / two;
                        var y2 = (-yp) - t / two;

                        var P1 = new EllipticCurvePoint(x, y1);
                        if (IsOnCurve(P1)) yield return P1;

                        if (y2 != y1)
                        {
                            var P2 = new EllipticCurvePoint(x, y2);
                            if (IsOnCurve(P2)) yield return P2;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Enumerate integral points with |x| ≤ xmax (a thin wrapper over RationalPoints with denMax=1).
        /// Returns O first, then affine integral points.
        /// </summary>
        public IEnumerable<EllipticCurvePoint> IntegralPoints(int xmax)
        {
            return RationalPoints(xmax, 1);
        }

        /// <summary>
        /// Returns the quadratic twist of this curve by a non-zero integer <paramref name="d"/>.
        /// The twisted curve E_d is isomorphic to E over the quadratic extension ℚ(√d),
        /// but generally not over ℚ (unless d is a square).
        /// </summary>
        public EllipticCurveQ QuadraticTwist(BigInteger d)
        {
            if (IsSingular) throw new InvalidOperationException("Cannot twist a singular curve.");
            if (d == 0) throw new ArgumentException("Twisting factor d cannot be zero.", nameof(d));

            var d_bi = d;               // d^1
            var d2 = d_bi * d_bi;       // d^2
            var d3 = d2 * d_bi;         // d^3

            // To ensure the twist is well-defined over ℚ without introducing square roots
            // into the a1/a3 terms, we first convert to the Short Weierstrass model.
            var shortCurve = this.ShortWeierstrass;

            // Apply the twisting formulas:
            // A' = A * d^2
            // B' = B * d^3
            // Note: Implicit conversion from BigInteger to BigRational is assumed for d2/d3.
            BigRational newA = shortCurve.A4 * new BigRational(d2);
            BigRational newB = shortCurve.A6 * new BigRational(d3);

            return new EllipticCurveQ(
                BigRational.Zero,
                BigRational.Zero,
                BigRational.Zero,
                newA,
                newB
            );
        }

        #endregion

        #region Torsion structure

        // Cached torsion set in ORIGINAL coordinates (including Infinity).
        private HashSet<EllipticCurvePoint> torsionPoints;

        /// <summary>
        /// Compute (and cache) the full set of rational torsion points of E(ℚ).
        /// Pipeline (no LMFDB):
        ///  0) Convert to short integral model Y^2 = X^3 + A'X + B'.
        ///  1) Use reductions at several good primes to restrict possible orders (Mazur admissible).
        ///  2) Apply Lutz–Nagell: Y^2 | |Δ'| and search via divisors to find integral torsion points.
        ///  3) Map the points back to the ORIGINAL model (inverse of the short/scale transform).
        /// Set contains Infinity and all affine torsion points; subsequent calls reuse the cache.
        /// </summary>
        public IEnumerable<EllipticCurvePoint> TorsionPoints
        {
            get
            {
                if (torsionPoints == null)
                {
                    // ---- 0) Short → integral short model ----
                    var Es = ShortWeierstrass; // y^2 = x^3 + A x + B
                    var A = Es.A4; var B = Es.A6;

                    var d = InternalMath.Lcm(A.Den, B.Den);
                    var Aint = A * BigRational.FromFraction(BigInteger.Pow(d, 4), BigInteger.One);
                    var Bint = B * BigRational.FromFraction(BigInteger.Pow(d, 6), BigInteger.One);
                    if (Aint.Den != BigInteger.One || Bint.Den != BigInteger.One)
                        throw new InvalidOperationException("Failed to obtain integral short model.");

                    var Ashort = Aint.Num;                  // integer A'
                    var Bshort = Bint.Num;                  // integer B'
                    var Eint = new EllipticCurveQ(
                        BigRational.Zero, BigRational.Zero, BigRational.Zero,
                        new BigRational(Ashort), new BigRational(Bshort));

                    // Δ' = −16(4A'^3 + 27B'^2). Nonsingular iff Δ' ≠ 0.
                    BigInteger Delta = -16 * (4 * BigInteger.Pow(Ashort, 3) + 27 * BigInteger.Pow(Bshort, 2));
                    if (Delta.IsZero) throw new InvalidOperationException("Singular curve.");

                    // ---- 1) Candidate orders via gcd of #E(F_p) for several good primes ----
                    // NOTE: if your C# doesn't support collection expressions, replace with new int[] { ... }.
                    int[] primes = [5, 7, 11, 13, 17, 19, 23, 29];
                    BigInteger gcdOrders = BigInteger.Zero;
                    for (int i = 0; i < primes.Length; i++)
                    {
                        int p = primes[i];
                        if (Delta % p == 0) continue; // only good reduction
                        var np = InternalMath.CountPointsFpShort(Ashort, Bshort, p);
                        gcdOrders = gcdOrders.IsZero ? new BigInteger(np) : BigInteger.GreatestCommonDivisor(gcdOrders, np);
                        if (gcdOrders.IsOne) break;
                    }
                    int[] mazur = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
                    var possibleOrders = new List<int>();
                    for (int i = 0; i < mazur.Length; i++)
                    {
                        int n = mazur[i];
                        if (gcdOrders.IsZero || (gcdOrders % n) == 0) possibleOrders.Add(n);
                    }

                    // ---- 2) Lutz–Nagell search on integral short model ----
                    var T = new HashSet<EllipticCurvePoint> { EllipticCurvePoint.Infinity };

                    // 2a) 2-torsion: X | B'  and X^3 + A'X + B' = 0
                    foreach (var x in InternalMath.EnumerateDivisorsAbs(Bshort))
                    {
                        if (InternalMath.EvalCubic(Ashort, Bshort, x) == 0)
                            T.Add(new EllipticCurvePoint(new BigRational(x), BigRational.Zero));
                        if (x != 0 && InternalMath.EvalCubic(Ashort, Bshort, -x) == 0)
                            T.Add(new EllipticCurvePoint(new BigRational(-x), BigRational.Zero));
                    }
                    // Edge: B' = 0 ⇒ X = ±sqrt(−A') if square
                    if (Bshort.IsZero)
                    {
                        var negA = BigInteger.Negate(Ashort);
                        if (negA.Sign >= 0)
                        {
                            var s = InternalMath.IntegerSqrt(negA);
                            if (s * s == negA && s != 0)
                            {
                                if (InternalMath.EvalCubic(Ashort, Bshort, s) == 0)
                                    T.Add(new EllipticCurvePoint(new BigRational(s), BigRational.Zero));
                                if (InternalMath.EvalCubic(Ashort, Bshort, -s) == 0)
                                    T.Add(new EllipticCurvePoint(new BigRational(-s), BigRational.Zero));
                            }
                        }
                    }

                    // 2b) odd torsion: y ≠ 0 and y^2 | |Δ'|
                    var factDelta = InternalMath.FactorAbs(Delta);
                    foreach (var y2 in InternalMath.EnumerateSquareDivisors(factDelta)) // y2 ≥ 1
                    {
                        if (y2.IsZero) continue;
                        var y = InternalMath.IntegerSqrt(y2); // exact sqrt

                        var C = Bshort - y2; // X^3 + A'X + (B' − y^2) = 0
                        if (C.IsZero)
                        {
                            var P1 = new EllipticCurvePoint(new BigRational(0), new BigRational(y));
                            var P2 = new EllipticCurvePoint(new BigRational(0), new BigRational(-y));
                            if (Eint.IsOnCurve(P1) && InternalMath.IsTorsionWithCandidates(Eint, P1, possibleOrders)) T.Add(P1);
                            if (Eint.IsOnCurve(P2) && InternalMath.IsTorsionWithCandidates(Eint, P2, possibleOrders)) T.Add(P2);

                            var negA = BigInteger.Negate(Ashort);
                            if (negA.Sign >= 0)
                            {
                                var s = InternalMath.IntegerSqrt(negA);
                                if (s * s == negA && s != 0)
                                {
                                    var xs = new[] { s, BigInteger.Negate(s) };
                                    for (int t = 0; t < xs.Length; t++)
                                    {
                                        var x = xs[t];
                                        if (InternalMath.EvalCubic(Ashort, Bshort, x) == y2)
                                        {
                                            var Q1 = new EllipticCurvePoint(new BigRational(x), new BigRational(y));
                                            var Q2 = new EllipticCurvePoint(new BigRational(x), new BigRational(-y));
                                            if (Eint.IsOnCurve(Q1) && InternalMath.IsTorsionWithCandidates(Eint, Q1, possibleOrders)) T.Add(Q1);
                                            if (Eint.IsOnCurve(Q2) && InternalMath.IsTorsionWithCandidates(Eint, Q2, possibleOrders)) T.Add(Q2);
                                        }
                                    }
                                }
                            }
                            continue;
                        }

                        var absC = C >= 0 ? C : -C;
                        foreach (var x in InternalMath.EnumerateDivisorsAbs(absC))
                        {
                            if (InternalMath.EvalCubic(Ashort, Bshort, x) == y2)
                            {
                                var R1 = new EllipticCurvePoint(new BigRational(x), new BigRational(y));
                                var R2 = new EllipticCurvePoint(new BigRational(x), new BigRational(-y));
                                if (Eint.IsOnCurve(R1) && InternalMath.IsTorsionWithCandidates(Eint, R1, possibleOrders)) T.Add(R1);
                                if (Eint.IsOnCurve(R2) && InternalMath.IsTorsionWithCandidates(Eint, R2, possibleOrders)) T.Add(R2);
                            }
                            if (x != 0)
                            {
                                var nx = -x;
                                if (InternalMath.EvalCubic(Ashort, Bshort, nx) == y2)
                                {
                                    var R1 = new EllipticCurvePoint(new BigRational(nx), new BigRational(y));
                                    var R2 = new EllipticCurvePoint(new BigRational(nx), new BigRational(-y));
                                    if (Eint.IsOnCurve(R1) && InternalMath.IsTorsionWithCandidates(Eint, R1, possibleOrders)) T.Add(R1);
                                    if (Eint.IsOnCurve(R2) && InternalMath.IsTorsionWithCandidates(Eint, R2, possibleOrders)) T.Add(R2);
                                }
                            }
                        }
                    }

                    // ---- 3) Map SHORT-INTEGRAL torsion back to ORIGINAL coordinates ----
                    var mapped = new HashSet<EllipticCurvePoint> { EllipticCurvePoint.Infinity };

                    var d2 = BigInteger.Pow(d, 2);
                    var d3 = BigInteger.Pow(d, 3);

                    // Inverse transform:
                    // short rational (X,Y) = (X_int/d^2, Y_int/d^3)
                    // a2' = A2 + (A1/2)^2
                    // original: x = X − a2'/3,   y = Y − (A1*x + A3)/2
                    var s1 = A1 / BigRational.FromInt(2);      // ORIGINAL A1
                    var a2p = A2 + s1 * s1;                     // a2'
                    var shift = a2p / BigRational.FromInt(3);   // a2'/3

                    foreach (var P in T)
                    {
                        if (P.IsInfinity) continue;

                        var X = P.X / new BigRational(d2);
                        var Y = P.Y / new BigRational(d3);

                        var x = X - shift;
                        var y = Y - (A1 * x + A3) / BigRational.FromInt(2);

                        var Porig = new EllipticCurvePoint(x, y);
                        if (IsOnCurve(Porig)) mapped.Add(Porig);
                    }

                    torsionPoints = mapped;
                }

                return torsionPoints;
            }
        }

        /// <summary>
        /// True iff the given point lies on the curve and is torsion (order &lt; ∞).
        /// </summary>
        public bool IsTorsionPoint(EllipticCurvePoint P)
        {
            if (!IsOnCurve(P)) throw new ArgumentException("Point is not on this curve.", nameof(P));
            return TorsionPoints.Contains(P);
        }

        /// <summary>
        /// Return the exact order of a torsion point or <c>null</c> if it has infinite order.
        /// The point must lie on the curve; Infinity is treated as order 1.
        /// </summary>
        public int? TorsionOrder(EllipticCurvePoint P)
        {
            if (!IsOnCurve(P)) throw new ArgumentException("Point is not on this curve.", nameof(P));
            if (P.IsInfinity) return 1;

            if (!TorsionPoints.Contains(P)) return null;

            var multiple = P;
            for (int n = 1; n <= TorsionPoints.Count(); n++)
            {
                if (multiple.IsInfinity) return n;
                multiple = Add(multiple, P);
            }

            throw new InvalidOperationException("Failed to determine torsion order for the given point.");
        }

        /// <summary>
        /// Group structure of the rational torsion subgroup E(ℚ)_tors inferred from TorsionPoints.
        /// Returns a label like "Z/1Z", "Z/4Z", or "Z/2Z x Z/4Z".
        /// </summary>
        public string TorsionStructure
        {
            get
            {
                // size includes Infinity and equals |E(ℚ)_tors|
                int size = TorsionPoints.Count();
                if (size == 1) return "Z/1Z";

                int twoTors = 0;
                var two = BigRational.FromInt(2);
                foreach (var P in TorsionPoints)
                {
                    if (P.IsInfinity) continue;

                    var doublingCondition = two * P.Y + A1 * P.X + A3;
                    if (doublingCondition.IsZero) twoTors++;
                }

                if (twoTors == 3)
                {
                    // Only non-cyclic case over ℚ: Z/2 × Z/2m (m=1..4), with |tors| = 4m.
                    int m = size / 4;
                    return $"Z/2Z x Z/{2 * m}Z";
                }
                else
                {
                    // Otherwise cyclic of order 'size'
                    return $"Z/{size}Z";
                }
            }
        }

        #endregion

        #region Isomorphism

        /// <summary>
        /// Return true iff this curve is Q–isomorphic to <paramref name="other"/>.
        /// Uses invariant scaling test (c4,c6,Δ) and outputs the scaling factor u (if requested).
        /// </summary>
        public bool IsIsomorphic(EllipticCurveQ other) => IsIsomorphic(other, out _);

        /// <summary>
        /// Return true iff this curve is Q–isomorphic to <paramref name="other"/>; also returns u s.t.
        /// c4_this = u^4 c4_other, c6_this = u^6 c6_other, Δ_this = u^12 Δ_other.
        /// </summary>
        public bool IsIsomorphic(EllipticCurveQ other, out BigRational u)
        {
            if (other is null) throw new ArgumentNullException(nameof(other));
            if (this.IsSingular || other.IsSingular)
                throw new InvalidOperationException("Isomorphism test requires nonsingular curves (Δ ≠ 0).");

            // “This” side (E): exact rationals
            var c4E = this.C4;
            var c6E = this.C6;
            var dE = this.Discriminant;

            // “Other” side (C): make an integral model by clearing denominators of a_i
            var (c4C, c6C, dC) = InternalMath.IntegralInvariants(other);

            return InternalMath.IsQIsomorphic(c4E, c6E, dE, c4C, c6C, dC, out u);
        }

        #endregion

        #region FromJInvariant

        /// <summary>
        /// Constructs an elliptic curve over ℚ with the specified j-invariant.
        /// </summary>
        public static EllipticCurveQ FromJInvariant(BigRational j)
        {
            // Case 1: CM curve with j = 0 (automorphism group order 6).
            // Standard minimal model: y^2 + y = x^3 (a3=1, others 0).
            // Discriminant = -27.
            if (j.IsZero)
                return new EllipticCurveQ(0, 0, 1, 0, 0);

            // Case 2: CM curve with j = 1728 (automorphism group order 4).
            // Standard minimal model: y^2 = x^3 - x (a4=-1, others 0).
            // Discriminant = 64.
            // Note: Assuming BigRational has implicit conversion or comparison with int.
            if (j == 1728)
                return new EllipticCurveQ(0, 0, 0, -1, 0);

            // Case 3: General case (j != 0 and j != 1728).
            // We construct a curve of the form y^2 = x^3 + Ax + A.
            // To satisfy the j-invariant equation, we set:
            // A = -27 * j / (4 * (j - 1728))
            var num = -BigRational.FromInt(27) * j;
            var den = BigRational.FromInt(4) * (j - 1728);

            var A = num / den;
            var B = A; // In this specific construction, we set B = A.

            // Returns y^2 = x^3 + Ax + A
            return new EllipticCurveQ(0, 0, 0, A, B);
        }

        #endregion

        #region Overrides

        /// <summary>Pretty printer that omits zero terms and formats signs compactly.</summary>
        public override string ToString()
        {
            var sb = new StringBuilder();

            // Left: y^2 + a1*x*y + a3*y
            sb.Append("y^2");
            AppendTerm(sb, A1, "x*y");
            AppendTerm(sb, A3, "y");

            // Separator
            sb.Append(" = ");

            // Right: x^3 + a2*x^2 + a4*x + a6
            sb.Append("x^3");
            AppendTerm(sb, A2, "x^2");
            AppendTerm(sb, A4, "x");
            AppendTerm(sb, A6, string.Empty);

            return sb.ToString();
        }

        /// <summary>Helper for ToString(): appends ± coeff*monomial, skipping zeros and eliding coeff=1 where appropriate.</summary>
        private static void AppendTerm(StringBuilder sb, BigRational coeff, string monomial)
        {
            if (coeff.IsZero) return;

            bool positive = coeff.Sign > 0;
            var abs = positive ? coeff : coeff.Negate();

            sb.Append(positive ? " + " : " - ");

            var showCoeff = string.IsNullOrEmpty(monomial) || abs != BigRational.One;

            if (showCoeff)
            {
                sb.Append(abs.ToString());
                if (!string.IsNullOrEmpty(monomial)) sb.Append('*');
            }

            if (!string.IsNullOrEmpty(monomial))
            {
                sb.Append(monomial);
            }
        }

        /// <summary>
        /// Value equality: finite points compare by coordinates; infinity compares only to infinity.
        /// </summary>
        public bool Equals(EllipticCurveQ other)
        {
            return A1 == other.A1
                && A2 == other.A2
                && A3 == other.A3
                && A4 == other.A4
                && A6 == other.A6;
        }

        /// <inheritdoc />
        public override bool Equals(object obj) => obj is EllipticCurveQ ec && Equals(ec);

        /// <inheritdoc />
        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 31 + A1.GetHashCode();
                hash = hash * 31 + A2.GetHashCode();
                hash = hash * 31 + A3.GetHashCode();
                hash = hash * 31 + A4.GetHashCode();
                hash = hash * 31 + A6.GetHashCode();
                return hash;
            }
        }

        #endregion
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        private EllipticCurveQ _reducedMinimalModel;
        private bool _conductorComputed;
        private BigInteger _conductor;
        private bool _minimalDiscriminantComputed;
        private BigInteger _minimalDiscriminant;

        /// <summary>
        /// The unique reduced global minimal model (integral a-invariants) isomorphic to this curve over Q.
        /// </summary>
        public EllipticCurveQ ReducedMinimalModel
        {
            get
            {
                if (_reducedMinimalModel is not null) return _reducedMinimalModel;

                // Compute integral invariants from an integral scaling of this curve, then recover
                // the reduced global minimal model via the Kraus–Laska–Connell algorithm.
                var inv = InternalMath.IntegralInvariants(this);
                BigInteger c4 = inv.c4;
                BigInteger c6 = inv.c6;

                var (ia1, ia2, ia3, ia4, ia6) = InternalMath.ReducedMinimalModelFromInvariants(c4, c6);
                _reducedMinimalModel = new EllipticCurveQ(new BigRational(ia1), new BigRational(ia2), new BigRational(ia3), new BigRational(ia4), new BigRational(ia6));
                return _reducedMinimalModel;
            }
        }

        /// <summary>
        /// Minimal discriminant Δ_min of the reduced minimal model (an integer).
        /// </summary>
        public BigInteger MinimalDiscriminant
        {
            get
            {
                if (_minimalDiscriminantComputed) return _minimalDiscriminant;
                var e = ReducedMinimalModel;
                var inv = InternalMath.ComputeInvariants(e.A1.Num, e.A2.Num, e.A3.Num, e.A4.Num, e.A6.Num);
                _minimalDiscriminant = inv.Delta;
                _minimalDiscriminantComputed = true;
                return _minimalDiscriminant;
            }
        }

        /// <summary>
        /// Conductor N of the curve (computed via Tate's algorithm on the reduced global minimal model).
        /// </summary>
        public BigInteger Conductor
        {
            get
            {
                if (_conductorComputed) return _conductor;
                _conductor = ComputeConductor();
                _conductorComputed = true;
                return _conductor;
            }
        }

        private BigInteger ComputeConductor()
        {
            var e = ReducedMinimalModel;
            var inv = InternalMath.ComputeInvariants(e.A1.Num, e.A2.Num, e.A3.Num, e.A4.Num, e.A6.Num);
            BigInteger delta = BigInteger.Abs(inv.Delta);
            if (delta.IsZero) throw new InvalidOperationException("Singular curve (Δ=0): conductor is undefined.");

            var factors = InternalMath.FactorAbs(delta);
            BigInteger N = BigInteger.One;
            foreach (var kv in factors)
            {
                BigInteger p = kv.Key;
                var (_, fp, _) = InternalMath.TateAlgorithm(e.A1.Num, e.A2.Num, e.A3.Num, e.A4.Num, e.A6.Num, p);
                if (fp > 0) N *= BigInteger.Pow(p, fp);
            }
            return N;
        }

        /// <summary>
        /// Compute (lower, upper) bounds for the Mordell–Weil rank of E(Q).
        /// 
        /// Lower bound is obtained by a rational point search with a numerical independence test
        /// using an approximate Néron–Tate height pairing.
        /// 
        /// Upper bound is a coarse fallback bound based on the size of S = primes dividing 2Δ_min
        /// (intended as a lightweight bound when full 2-descent is not implemented).
        /// </summary>
        /// <param name="searchNumMax">Maximum absolute numerator for x in the rational search.</param>
        /// <param name="searchDenMax">Maximum denominator for x in the rational search.</param>
        /// <param name="maxPoints">Maximum number of rational points to consider from the search output.</param>
        /// <param name="heightDoublings">Number of doublings used in the height approximation (>= 4 recommended).</param>
        public (int Lower, int Upper) RankBounds(int searchNumMax = 50, int searchDenMax = 50, int maxPoints = 512, int heightDoublings = 10)
        {
            if (searchNumMax < 1) throw new ArgumentOutOfRangeException(nameof(searchNumMax));
            if (searchDenMax < 1) throw new ArgumentOutOfRangeException(nameof(searchDenMax));
            if (maxPoints < 16) throw new ArgumentOutOfRangeException(nameof(maxPoints));
            if (heightDoublings < 0) throw new ArgumentOutOfRangeException(nameof(heightDoublings));

            int lower = RankLowerBoundFromSearch(searchNumMax, searchDenMax, maxPoints, heightDoublings);

            // Coarse bound (very safe but often loose): 2*ω(2*Δ_min) + 2.
            BigInteger absMinDisc = BigInteger.Abs(MinimalDiscriminant);
            int omega = InternalMath.FactorAbs(absMinDisc == 0 ? 0 : 2 * absMinDisc).Count;
            int upper = 2 * omega + 2;

            // Tighter bound when E has a rational 2-torsion point: 2-isogeny descent.
            // We take the minimum over all rational 2-torsion points (when there are several),
            // and then intersect with the coarse bound above.
            int? twoIsoUpper = TryRankUpperBoundByTwoIsogenyDescent();
            if (twoIsoUpper.HasValue) upper = Math.Min(upper, twoIsoUpper.Value);

            if (upper < lower) upper = lower;
            return (lower, upper);
        }

        private int? TryRankUpperBoundByTwoIsogenyDescent()
        {
            // Work on the reduced minimal model to keep integral coefficients small.
            var E = this.ReducedMinimalModel;

            if (E.A1.Den != 1 || E.A2.Den != 1 || E.A3.Den != 1 || E.A4.Den != 1 || E.A6.Den != 1)
                return null;

            BigInteger a1 = E.A1.Num, a2 = E.A2.Num, a3 = E.A3.Num, a4 = E.A4.Num, a6 = E.A6.Num;

            // Standard invariants
            BigInteger b2 = a1 * a1 + 4 * a2;
            BigInteger b4 = a1 * a3 + 2 * a4;
            BigInteger b6 = a3 * a3 + 4 * a6;

            // Rational 2-torsion x-coordinates are exactly the rational roots of:
            //   4 x^3 + b2 x^2 + 2 b4 x + b6 = 0.
            var roots = RationalRootsOfTwoTorsionCubic(b2, b4, b6);
            if (roots.Count == 0) return null;

            int best = int.MaxValue;

            // Completed-square cubic:
            //   x^3 + (b2/4) x^2 + (b4/2) x + (b6/4).
            var A = new BigRational(b2, 4);
            var B = new BigRational(b4, 2);

            foreach (var r in roots)
            {
                // Shift x = X + r so r becomes the root at X=0.
                // Then the model becomes: y^2 = X^3 + a X^2 + b X.
                BigRational aShift = A + 3 * r;
                BigRational bShift = B + 2 * A * r + 3 * r * r;

                // Choose u = 2^e so that aShift*u^2 and bShift*u^4 are integers.
                BigInteger u = MinimalPowerOfTwoScale(aShift.Den, bShift.Den);
                BigInteger u2 = u * u;
                BigInteger u4 = u2 * u2;

                BigInteger aIso = (aShift.Num * u2) / aShift.Den;
                BigInteger bIso = (bShift.Num * u4) / bShift.Den;

                if (bIso.IsZero) continue; // singular (shouldn't happen for nonsingular E)

                int? ub = RankUpperBoundFromTwoIsogenyDescent(aIso, bIso);
                if (ub.HasValue) best = Math.Min(best, ub.Value);
            }

            return best == int.MaxValue ? (int?)null : best;
        }

        private static List<BigRational> RationalRootsOfTwoTorsionCubic(BigInteger b2, BigInteger b4, BigInteger b6)
        {
            // P(x) = 4 x^3 + b2 x^2 + 2 b4 x + b6.
            // Return rational roots (in lowest terms). Any rational root has denominator dividing 4.
            var roots = new List<BigRational>();

            if (b6.IsZero)
                roots.Add(BigRational.Zero);

            if (!b6.IsZero)
            {
                foreach (var dAbs in InternalMath.EnumerateDivisorsAbs(b6))
                {
                    if (dAbs.IsZero) continue;

                    foreach (int q in new[] { 1, 2, 4 })
                    {
                        BigInteger qBI = q;

                        foreach (var p in new[] { dAbs, -dAbs })
                        {
                            if (IsTwoTorsionCubicRoot(p, qBI, b2, b4, b6))
                            {
                                var r = new BigRational(p, qBI);
                                if (!roots.Contains(r)) roots.Add(r);
                                if (roots.Count == 3) return roots;
                            }
                        }
                    }
                }
            }
            else
            {
                // b6 == 0: P(x) = x(4x^2 + b2 x + 2b4). Handle remaining roots explicitly.
                BigInteger disc = b2 * b2 - 32 * b4;
                if (disc.Sign >= 0)
                {
                    BigInteger s = InternalMath.IntegerSqrt(disc);
                    if (s * s == disc)
                    {
                        foreach (var sign in new[] { BigInteger.One, -BigInteger.One })
                        {
                            var r = new BigRational(-b2 + sign * s, 8);
                            if (!roots.Contains(r)) roots.Add(r);
                        }
                    }
                }
            }

            return roots;
        }

        private static bool IsTwoTorsionCubicRoot(BigInteger p, BigInteger q, BigInteger b2, BigInteger b4, BigInteger b6)
        {
            // Evaluate q^3 * P(p/q) to avoid rationals:
            //   4 p^3 + b2 p^2 q + 2 b4 p q^2 + b6 q^3
            BigInteger val =
                4 * p * p * p
                + b2 * p * p * q
                + 2 * b4 * p * q * q
                + b6 * q * q * q;

            return val.IsZero;
        }

        private static BigInteger MinimalPowerOfTwoScale(BigInteger denA, BigInteger denB)
        {
            // We need u = 2^e such that denA | u^2 and denB | u^4.
            // In our workflow (reduced minimal model + completing the square), denominators are powers of 2.
            int eA = TwoAdicExponent(denA);
            int eB = TwoAdicExponent(denB);

            // If a denominator is not a pure power of two, fall back to a conservative choice.
            if ((denA >> eA) != BigInteger.One || (denB >> eB) != BigInteger.One)
            {
                BigInteger u = denA * denB;
                return u.IsZero ? BigInteger.One : u;
            }

            int e = Math.Max((eA + 1) / 2, (eB + 3) / 4);
            return BigInteger.One << e;
        }

        private static int TwoAdicExponent(BigInteger n)
        {
            int e = 0;
            while (n > 1 && (n & BigInteger.One) == BigInteger.Zero)
            {
                n >>= 1;
                e++;
            }
            return e;
        }

        private static int? RankUpperBoundFromTwoIsogenyDescent(BigInteger a, BigInteger b)
        {
            // E: y^2 = x^3 + a x^2 + b x with kernel <(0,0)>.
            if (b.IsZero) return null;

            int? dim1 = TwoIsogenySelmerDimension(a, b);
            if (!dim1.HasValue) return null;

            BigInteger a2 = -2 * a;
            BigInteger b2 = a * a - 4 * b;
            if (b2.IsZero) return null;

            int? dim2 = TwoIsogenySelmerDimension(a2, b2);
            if (!dim2.HasValue) return null;

            int upper = dim1.Value + dim2.Value - 2;
            return upper < 0 ? 0 : upper;
        }

        private static int? TwoIsogenySelmerDimension(BigInteger a, BigInteger b)
        {
            if (b.IsZero) return null;

            var primes = InternalMath.FactorAbs(b).Keys.ToList();
            primes.Sort((x, y) => x.CompareTo(y));

            // Selmer elements live in the squareclass group generated by -1 and primes dividing b.
            int n = primes.Count + 1;

            // Prevent blow-ups for highly composite b.
            if (n > 30) return null;

            BigInteger discPart = 2 * b * (a * a - 4 * b);
            var checkPrimes = InternalMath.FactorAbs(discPart).Keys
                .Where(p => p >= 2 && p <= int.MaxValue)
                .Select(p => (int)p)
                .ToList();
            checkPrimes.Sort();

            int total = 1 << n;

            // Gaussian elimination basis (bit i corresponds to generator i).
            ulong[] basis = new ulong[n];
            int rank = 0;

            for (int mask = 0; mask < total; mask++)
            {
                BigInteger d = BigInteger.One;
                if ((mask & 1) != 0) d = -d;
                for (int i = 0; i < primes.Count; i++)
                {
                    if ((mask & (1 << (i + 1))) != 0) d *= primes[i];
                }

                if (!RealSolvableIsogenyQuartic(a, b, d)) continue;

                bool ok = true;
                foreach (int p in checkPrimes)
                {
                    if (!HasLocalPointIsogenyQuartic(a, b, d, p))
                    {
                        ok = false;
                        break;
                    }
                }
                if (!ok) continue;

                // Add vector to basis.
                ulong v = (ulong)mask;
                for (int bit = n - 1; bit >= 0; bit--)
                {
                    if (((v >> bit) & 1UL) == 0) continue;
                    if (basis[bit] != 0UL) v ^= basis[bit];
                    else
                    {
                        basis[bit] = v;
                        rank++;
                        break;
                    }
                }
            }

            return rank;
        }

        private static bool RealSolvableIsogenyQuartic(BigInteger a, BigInteger b, BigInteger d)
        {
            if (d.Sign > 0) return true;

            // d < 0, so the quartic is a downward-opening quadratic in t = z^2 >= 0:
            //   w^2 = d t^2 + a t + (b/d).
            if (a.Sign <= 0)
            {
                // Maximum on t >= 0 occurs at t=0.
                // Need b/d >= 0. Since d < 0, this is equivalent to b <= 0.
                return b.Sign <= 0;
            }

            // Vertex t0 = -a/(2d) is positive. Need max >= 0.
            // Condition simplifies to: 4*b - a^2 <= 0 (because d < 0).
            return (4 * b - a * a) <= 0;
        }

        private static bool HasLocalPointIsogenyQuartic(BigInteger a, BigInteger b, BigInteger d, int p)
        {
            if (p == 2)
            {
                // 2-adic obstructions are the most common; use a stronger modulus.
                const int k = 6; // 2^6 = 64
                return HasSolutionModulo(a, b, d, 1 << k, reciprocal: false) || HasSolutionModulo(a, b, d, 1 << k, reciprocal: true);
            }

            // For very large primes, a full search mod p is expensive; skipping is safe (may only loosen the bound).
            if (p > 50000) return true;

            int kOdd = (p <= 19) ? 2 : 1;

            long mLong = 1;
            for (int i = 0; i < kOdd; i++)
                mLong *= p;

            // Guardrail against big allocations/time.
            if (mLong > 2_000_000)
                mLong = p;

            int m = (int)mLong;

            return HasSolutionModulo(a, b, d, m, reciprocal: false) || HasSolutionModulo(a, b, d, m, reciprocal: true);
        }

        private static bool HasSolutionModulo(BigInteger a, BigInteger b, BigInteger d, int m, bool reciprocal)
        {
            if (m <= 1) return true;

            int dMod = ModToInt(d, m);
            int aMod = ModToInt(a, m);
            int bMod = ModToInt(b, m);

            int d2Mod = (int)((long)dMod * dMod % m);
            int adMod = (int)((long)aMod * dMod % m);

            // possibleRhs[r] == true iff r == d*w^2 (mod m) for some w.
            var possibleRhs = new bool[m];
            for (int w = 0; w < m; w++)
            {
                int w2 = (int)((long)w * w % m);
                int r = (int)((long)dMod * w2 % m);
                possibleRhs[r] = true;
            }

            for (int z = 0; z < m; z++)
            {
                int z2 = (int)((long)z * z % m);
                int z4 = (int)((long)z2 * z2 % m);

                long rhs = reciprocal
                    ? ((long)bMod * z4 + (long)adMod * z2 + d2Mod)
                    : ((long)d2Mod * z4 + (long)adMod * z2 + bMod);

                int r = (int)(rhs % m);
                if (r < 0) r += m;

                if (possibleRhs[r]) return true;
            }

            return false;
        }

        private static int ModToInt(BigInteger x, int m)
        {
            BigInteger r = x % m;
            if (r.Sign < 0) r += m;
            return (int)r;
        }

        private int RankLowerBoundFromSearch(int searchNumMax, int searchDenMax, int maxPoints, int heightDoublings)
        {
            // Enumerate rational points in a box. We only process up to maxPoints points to keep
            // the numerical linear algebra bounded.
            //
            // Torsion filter:
            // Over Q, the torsion subgroup is classified (Mazur). In particular, every torsion point
            // has order dividing lcm(1,2,3,4,5,6,7,8,9,10,12,16) = 5040, hence 5040*P = O  iff  P is torsion.
            //
            // This avoids computing the full torsion subgroup (your TorsionPoints cache), which can
            // dominate runtime when rank bounds are requested repeatedly.
            const int torsionKiller = 5040;

            var basis = new List<EllipticCurvePoint>();
            var heights = new List<double>();
            var gram = new List<double[]>(); // lower-triangular storage per row (row i has i+1 entries)

            int processed = 0;
            foreach (var p in RationalPoints(searchNumMax, searchDenMax))
            {
                if (processed++ >= maxPoints) break;
                if (p.IsInfinity) continue;
                if (Multiply(p, torsionKiller).IsInfinity) continue;

                // Approximate canonical height.
                double hp = CanonicalHeightApprox(p, heightDoublings);
                if (!(hp > 1e-10)) continue;

                if (basis.Count == 0)
                {
                    basis.Add(p);
                    heights.Add(hp);
                    gram.Add([hp]);
                    continue;
                }

                // Compute pairing vector with existing basis.
                int r = basis.Count;
                var v = new double[r];
                for (int i = 0; i < r; i++)
                {
                    double hbi = heights[i];
                    double hsum = CanonicalHeightApprox(Add(p, basis[i]), heightDoublings);
                    v[i] = 0.5 * (hsum - hp - hbi);
                }

                if (!NumericallyIndependent(hp, v, gram)) continue;

                // Extend Gram matrix.
                var newRow = new double[r + 1];
                for (int i = 0; i < r; i++) newRow[i] = v[i];
                newRow[r] = hp;

                basis.Add(p);
                heights.Add(hp);
                gram.Add(newRow);
            }

            return basis.Count;
        }

        private static bool NumericallyIndependent(double hp, double[] v, List<double[]> gram)
        {
            // Given Gram matrix G of current basis and vector v = <P, basis>,
            // compute residual = <P,P> - v^T G^{-1} v.
            // If residual is sufficiently positive, treat as independent.

            int r = gram.Count;
            if (r == 0) return hp > 1e-10;

            // Build dense matrix for small r.
            var G = new double[r, r];
            for (int i = 0; i < r; i++)
            {
                var row = gram[i];
                for (int j = 0; j <= i; j++)
                {
                    double gij = row[j];
                    G[i, j] = gij;
                    G[j, i] = gij;
                }
            }

            // Solve G * x = v using Gaussian elimination with partial pivoting.
            var x = SolveLinearSystem(G, v);
            if (x is null) return false;
            double proj = 0;
            for (int i = 0; i < r; i++) proj += v[i] * x[i];
            double residual = hp - proj;
            return residual > 1e-8;
        }

        private static double[]? SolveLinearSystem(double[,] A, double[] b)
        {
            // Dense Gaussian elimination for small matrices.
            int n = b.Length;
            var M = (double[,])A.Clone();
            var x = (double[])b.Clone();

            for (int k = 0; k < n; k++)
            {
                // Pivot.
                int piv = k;
                double best = Math.Abs(M[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    double v = Math.Abs(M[i, k]);
                    if (v > best)
                    {
                        best = v;
                        piv = i;
                    }
                }

                if (best < 1e-14) return null; // singular / ill-conditioned

                if (piv != k)
                {
                    // swap rows k and piv
                    for (int j = k; j < n; j++) (M[k, j], M[piv, j]) = (M[piv, j], M[k, j]);
                    (x[k], x[piv]) = (x[piv], x[k]);
                }

                double diag = M[k, k];
                for (int j = k; j < n; j++) M[k, j] /= diag;
                x[k] /= diag;

                for (int i = 0; i < n; i++)
                {
                    if (i == k) continue;
                    double f = M[i, k];
                    if (Math.Abs(f) < 1e-18) continue;
                    for (int j = k; j < n; j++) M[i, j] -= f * M[k, j];
                    x[i] -= f * x[k];
                }
            }

            return x;
        }

        private double CanonicalHeightApprox(EllipticCurvePoint p, int doublings)
        {
            // Very lightweight approximation:
            //   ˆh(P) ≈ average_{i=2..doublings} h(x(2^i P)) / 4^i
            // where h(x) is the logarithmic naive height of x.
            //
            // This is suitable for numerical independence tests; it is not intended as a
            // high-precision canonical height implementation.

            if (p.IsInfinity) return 0.0;
            if (doublings <= 0) return 0.0;

            EllipticCurvePoint q = p;
            double sum = 0.0;
            int cnt = 0;

            for (int i = 0; i <= doublings; i++)
            {
                if (q.IsInfinity) break;

                if (i >= 2)
                {
                    double hx = NaiveLogHeight(q.X);
                    sum += hx / Math.Pow(4.0, i);
                    cnt++;
                }

                q = Double(q);
            }

            if (cnt == 0) return NaiveLogHeight(p.X);
            return sum / cnt;
        }

        private static double NaiveLogHeight(BigRational x)
        {
            // h(x) = log(max(|num|, den))
            BigInteger n = BigInteger.Abs(x.Num);
            BigInteger d = x.Den;

            if (n.IsZero) return BigInteger.Log(d);
            double lnN = BigInteger.Log(n);
            double lnD = BigInteger.Log(d);
            return Math.Max(lnN, lnD);
        }
    }
}

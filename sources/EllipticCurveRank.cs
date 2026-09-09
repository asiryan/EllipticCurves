using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>
        /// Compute unconditional algebraic rank bounds locally. With rational 2-torsion,
        /// use descent by 2-isogeny, exact quartic point searches and local obstructions.
        /// Otherwise use general binary-quartic 2-descent. Good-reduction Kummer
        /// characters also certify lower bounds from several rational points.
        /// Equality of the bounds certifies the exact rank; no BSD or parity assumption is used.
        /// </summary>
        /// <param name="searchBound">Nonnegative bound on each primitive quartic coordinate.
        /// Also bounds both numerator and denominator of the searched x coordinates on the minimal model.
        /// Zero disables point search. Increasing this can improve the lower bound.</param>
        /// <param name="maxSquareClasses">Maximum number of signed square classes per isogeny.
        /// Exceeding this limit throws NotSupportedException in isogeny descent.
        /// In general descent this bounds covering classes and exhaustion gives a null upper bound.</param>
        /// <param name="cancellationToken">Cancels factorization, enumeration and point searches.</param>
        public RankBounds GetRankBounds(int searchBound = 32, int maxSquareClasses = 65536,
            CancellationToken cancellationToken = default)
        {
            if (searchBound < 0 || searchBound == int.MaxValue) throw new ArgumentOutOfRangeException(nameof(searchBound));
            if (maxSquareClasses < 2) throw new ArgumentOutOfRangeException(nameof(maxSquareClasses));
            return GetRankBounds(new RankComputationOptions { SearchBound = searchBound, MaxSquareClasses = maxSquareClasses }, cancellationToken);
        }

        /// <summary>
        /// Compute unconditional lower and upper rank bounds with explicit work limits.
        /// An incomplete descent returns a null upper bound; point-search exhaustion
        /// never invalidates already proved bounds. No BSD, GRH or parity assumption is used.
        /// </summary>
        public RankBounds GetRankBounds(RankComputationOptions options, CancellationToken cancellationToken = default)
        {
            if (options == null) throw new ArgumentNullException(nameof(options));
            options = options.Snapshot();
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("A singular curve has no Mordell-Weil rank.");
            var e = GetGlobalMinimalModel(cancellationToken);
            var a = e.B2.Num;
            var b = 8 * e.B4.Num;
            var c = 16 * e.B6.Num;
            bool hasTwoTorsion = TryCubicRoot(a, b, c, cancellationToken, out var root);
            int torsionDimension = 0;
            if (hasTwoTorsion)
            {
                var quadraticDiscriminant = (a + root) * (a + root) - 4 * (b + a * root + root * root);
                torsionDimension = quadraticDiscriminant > 0 && BigRational.IsSquare(new BigRational(quadraticDiscriminant), out _) ? 2 : 1;
            }
            var budget = new DescentBudget(options, cancellationToken);
            var pointRank = new RationalPointRank(e, torsionDimension, budget);
            pointRank.Search();
            bool general = !hasTwoTorsion || options.PreferGeneralTwoDescent;
            GeneralTwoDescent descent = null;
            int firstImageLower = 0, secondImageLower = 0;
            try
            {
                int lower, upper; int? selmerDimension = null;
                if (general)
                {
                    descent = new GeneralTwoDescent(e, torsionDimension, pointRank, budget);
                    selmerDimension = descent.ComputeSelmerDimension();
                    lower = descent.LowerBound;
                    upper = selmerDimension.Value - torsionDimension;
                }
                else
                {
                    b += 2 * a * root + 3 * root * root;
                    a += 3 * root;
                    foreach (var p in Factor(b, cancellationToken).Keys)
                    {
                        var p2 = p * p; var p4 = p2 * p2;
                        while (a % p2 == 0 && b % p4 == 0) { a /= p2; b /= p4; }
                    }
                    var dualB = a * a - 4 * b;
                    firstImageLower = BigRational.IsSquare(new BigRational(b), out _) ? 0 : 1;
                    secondImageLower = BigRational.IsSquare(new BigRational(dualB), out _) ? 0 : 1;
                    var first = IsogenyImageBounds(a, b, budget, dimension => firstImageLower = dimension);
                    var second = IsogenyImageBounds(-2 * a, dualB, budget, dimension => secondImageLower = dimension);
                    lower = Math.Max(pointRank.LowerBound, first.lower + second.lower - 2);
                    upper = first.upper + second.upper - 2;
                }
                if (upper < lower) throw new InvalidOperationException("Inconsistent certified rank bounds.");
                string reason = lower == upper ? "The proved lower and upper bounds agree."
                    : "The descent upper bound exceeds the proved lower bound; no claim of global solubility or of vanishing Sha is made.";
                if (budget.PointSearchExhausted) reason += " MaxPointSearchWork was reached.";
                return new RankBounds(lower, upper, !general, general, selmerDimension, reason, budget.Work, budget.PointWork);
            }
            catch (DescentLimitException ex)
            {
                int lower = descent?.LowerBound ?? Math.Max(pointRank.LowerBound, firstImageLower + secondImageLower - 2);
                string reason = ex.Message + (budget.PointSearchExhausted ? " MaxPointSearchWork was reached." : "");
                return new RankBounds(lower, null, !general, general,
                    reason: reason, descentWork: budget.Work, pointWork: budget.PointWork);
            }
        }

        private static bool TryCubicRoot(BigInteger a, BigInteger b, BigInteger c,
            CancellationToken token, out BigInteger root)
        {
            root = 0;
            if (c.IsZero) return true;
            // A rational root of a monic integral polynomial is an integral divisor of c.
            foreach (var d in Divisors(Factor(c, token)))
            {
                token.ThrowIfCancellationRequested();
                if (((d + a) * d + b) * d + c == 0) { root = d; return true; }
                if (((-d + a) * d - b) * d + c == 0) { root = -d; return true; }
            }
            return false;
        }

        private static (int lower, int upper) IsogenyImageBounds(BigInteger a, BigInteger b,
            DescentBudget budget, Action<int> recordLower)
        {
            var token = budget.Token;
            int maxClasses = budget.Options.MaxSquareClasses;
            var factors = Factor(b, token).OrderBy(pair => pair.Key).ToArray();
            long classCount = 2;
            foreach (var unused in factors)
            {
                classCount *= 2;
                if (classCount > maxClasses)
                    throw new NotSupportedException("The 2-isogeny descent exceeds maxSquareClasses. Increase the limit or use a smaller model.");
            }
            int count = (int)classCount;
            var values = new BigInteger[count];
            values[0] = 1;
            values[1] = -1;
            int used = 2, torsionClass = b.Sign < 0 ? 1 : 0;
            for (int i = 0; i < factors.Length; i++)
            {
                if ((factors[i].Value & 1) != 0) torsionClass |= 1 << (i + 1);
                for (int j = 0; j < used; j++) values[used + j] = values[j] * factors[i].Key;
                used *= 2;
            }
            var basis = new int[factors.Length + 1];
            AddSquareClass(basis, torsionClass);
            int survivors = 0;
            var badPrimes = Factor(2 * b * (a * a - 4 * b), token).Keys.OrderBy(p => p).ToArray();
            foreach (int mask in Enumerable.Range(0, count))
            {
                budget.Step();
                var d = values[mask];
                var other = Divide(b, d);
                var quartic = new BinaryQuartic(d, 0, a, 0, other);
                bool known = ReducesToZero(basis, mask);
                if (!known && !QuarticLocalSolubility.Everywhere(quartic, badPrimes, budget)) continue;
                survivors++;
                if (!known && quartic.TryPoint(budget, out _, out _, out _))
                {
                    AddSquareClass(basis, mask);
                    recordLower(basis.Count(x => x != 0));
                }
            }
            if (survivors == 0) throw new InvalidOperationException("The trivial descent class was eliminated.");
            if ((survivors & (survivors - 1)) != 0) throw new InvalidOperationException("The isogeny Selmer group does not have power-of-two order.");
            int upper = 0;
            for (int size = survivors; size > 1; size >>= 1) upper++;
            int lower = basis.Count(x => x != 0);
            if (lower > upper) throw new InvalidOperationException("Inconsistent isogeny image bounds.");
            return (lower, upper);
        }

        private static bool ReducesToZero(int[] basis, int value)
        {
            for (int i = basis.Length - 1; i >= 0; i--)
                if ((value & (1 << i)) != 0) value ^= basis[i];
            return value == 0;
        }

        private static void AddSquareClass(int[] basis, int value)
        {
            for (int i = basis.Length - 1; i >= 0; i--)
            {
                if ((value & (1 << i)) == 0) continue;
                if (basis[i] == 0) { basis[i] = value; return; }
                value ^= basis[i];
            }
        }
    }
}

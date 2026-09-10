using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>All rational Q satisfying [n]Q = point, on the input model. Requires n != 0.
        /// An empty list proves nondivisibility; a work or degree limit throws instead.</summary>
        public IReadOnlyList<EllipticCurvePoint> GetDivisionPoints(EllipticCurvePoint point, int n,
            PointDivisionOptions options = null, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("Division requires a nonsingular curve.");
            if (!IsOnCurve(point)) throw new ArgumentException("Point is not on this curve.", nameof(point));
            if (n == 0 || n == int.MinValue) throw new ArgumentOutOfRangeException(nameof(n));
            var limits = (options ?? new PointDivisionOptions()).Snapshot();
            int order = Math.Abs(n);
            if (order == 1) return Array.AsReadOnly(new[] { n > 0 ? point : Negate(point) });
            if (point.IsInfinity)
            {
                var torsion = TorsionPoints.ToArray();
                cancellationToken.ThrowIfCancellationRequested();
                return Array.AsReadOnly(torsion.Where(p => Multiply(p, order).IsInfinity).ToArray());
            }
            if ((long)order * order > limits.MaxDivisionDegree) throw new ArithmeticException("Division degree limit reached.");
            var map = GetShortModelIsomorphism(); var e = map.Target;
            var target = map.Map(n > 0 ? point : Negate(point));
            var budget = new DescentBudget(new RankComputationOptions { MaxDescentWork = limits.MaxWork }, cancellationToken);
            try
            {
                var arithmetic = new DivisionPolynomials(e.A4, e.A6, budget);
                var multiplication = arithmetic.Multiplication(order);
                var equation = arithmetic.Sub(multiplication.numerator, arithmetic.Scale(multiplication.denominator, target.X));
                var result = new HashSet<EllipticCurvePoint>();
                foreach (var x in arithmetic.RationalRoots(equation))
                {
                    budget.Step();
                    if (!BigRational.IsSquare(x * x * x + e.A4 * x + e.A6, out var y)) continue;
                    foreach (var candidate in new[] { new EllipticCurvePoint(x, y), new EllipticCurvePoint(x, -y) })
                        if (e.Multiply(candidate, order).Equals(target)) result.Add(map.MapBack(candidate));
                }
                return Array.AsReadOnly(result.OrderBy(p => p.X).ThenBy(p => p.Y).ToArray());
            }
            catch (DescentLimitException ex) { throw new ArithmeticException("Division work limit reached; divisibility is unresolved.", ex); }
        }

        /// <summary>Find a rational Q with [n]Q = point. False means proved nondivisibility; limits throw.
        /// When false, quotient is set to infinity and is not a solution.</summary>
        public bool TryDividePoint(EllipticCurvePoint point, int n, out EllipticCurvePoint quotient,
            PointDivisionOptions options = null, CancellationToken cancellationToken = default)
        {
            var points = GetDivisionPoints(point, n, options, cancellationToken);
            quotient = points.Count == 0 ? EllipticCurvePoint.Infinity : points[0];
            return points.Count != 0;
        }

        private WeierstrassIsomorphism GetShortModelIsomorphism()
        {
            var r = -B2 / 12;
            return ChangeModel(1, r, -A1 / 2, -(A3 + r * A1) / 2);
        }
    }
}

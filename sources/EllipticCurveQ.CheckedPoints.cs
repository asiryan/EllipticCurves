using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Create an affine rational point, rejecting singular models and coordinates outside this curve.</summary>
        /// <exception cref="InvalidOperationException">The curve is singular.</exception>
        /// <exception cref="ArgumentException">The coordinates do not satisfy this curve's equation.</exception>
        public EllipticCurvePoint CreatePoint(BigRational x, BigRational y)
        {
            ValidateCheckedCurve();
            var point = new EllipticCurvePoint(x, y);
            if (!IsOnCurve(point)) throw new ArgumentException("The coordinates do not lie on this curve.");
            return point;
        }

        /// <summary>Try to create an affine rational point on this nonsingular curve.
        /// Returns false for singular models or coordinates outside the curve, with an unused Infinity output.</summary>
        public bool TryCreatePoint(BigRational x, BigRational y, out EllipticCurvePoint point)
        {
            point = EllipticCurvePoint.Infinity;
            if (IsSingular) return false;
            var candidate = new EllipticCurvePoint(x, y);
            if (!IsOnCurve(candidate)) return false;
            point = candidate;
            return true;
        }

        /// <summary>Add points after checking nonsingularity and the membership of both inputs, including before identity shortcuts.</summary>
        public EllipticCurvePoint AddChecked(EllipticCurvePoint p, EllipticCurvePoint q)
        {
            ValidateCheckedCurve();
            ValidateCheckedPoint(p, nameof(p));
            ValidateCheckedPoint(q, nameof(q));
            return Add(p, q);
        }

        /// <summary>Double a point after checking nonsingularity and membership in this curve.</summary>
        public EllipticCurvePoint DoubleChecked(EllipticCurvePoint p)
        {
            ValidateCheckedCurve();
            ValidateCheckedPoint(p, nameof(p));
            return Double(p);
        }

        /// <summary>Negate a point after checking nonsingularity and membership in this curve.</summary>
        public EllipticCurvePoint NegateChecked(EllipticCurvePoint p)
        {
            ValidateCheckedCurve();
            ValidateCheckedPoint(p, nameof(p));
            return Negate(p);
        }

        /// <summary>Subtract points after checking nonsingularity and the membership of both inputs.</summary>
        public EllipticCurvePoint SubtractChecked(EllipticCurvePoint p, EllipticCurvePoint q)
        {
            ValidateCheckedCurve();
            ValidateCheckedPoint(p, nameof(p));
            ValidateCheckedPoint(q, nameof(q));
            return Subtract(p, q);
        }

        /// <summary>Multiply a point after checking nonsingularity and membership, even when n is zero.
        /// Checks the input once and observes cancellation between subsequent group operations.</summary>
        public EllipticCurvePoint MultiplyChecked(EllipticCurvePoint p, BigInteger n,
            CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            ValidateCheckedCurve();
            ValidateCheckedPoint(p, nameof(p));
            return Multiply(p, n, cancellationToken);
        }

        private void ValidateCheckedCurve()
        {
            if (IsSingular) throw new InvalidOperationException("Checked point operations require a nonsingular elliptic curve.");
        }

        private void ValidateCheckedPoint(EllipticCurvePoint point, string parameterName)
        {
            if (!IsOnCurve(point)) throw new ArgumentException("Point is not on this curve.", parameterName);
        }
    }
}

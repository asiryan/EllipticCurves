using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>A nonsingular general Weierstrass curve over a specified finite field F_(p^k).
    /// Supports characteristics two and three. Counting enumerates all q^2 affine coordinate pairs.</summary>
    public sealed class EllipticCurveFq : IEquatable<EllipticCurveFq>
    {
        /// <summary>Base field, including its defining polynomial.</summary>
        public FiniteField Field { get; }
        /// <summary>Coefficient of xy.</summary>
        public FiniteFieldElement A1 { get; }
        /// <summary>Coefficient of x^2.</summary>
        public FiniteFieldElement A2 { get; }
        /// <summary>Coefficient of y.</summary>
        public FiniteFieldElement A3 { get; }
        /// <summary>Coefficient of x.</summary>
        public FiniteFieldElement A4 { get; }
        /// <summary>Constant coefficient.</summary>
        public FiniteFieldElement A6 { get; }
        /// <summary>Nonzero discriminant in the base field.</summary>
        public FiniteFieldElement Discriminant { get; }
        /// <summary>Invariant c4.</summary>
        public FiniteFieldElement C4 { get; }
        /// <summary>Invariant c6.</summary>
        public FiniteFieldElement C6 { get; }
        /// <summary>j-invariant.</summary>
        public FiniteFieldElement JInvariant => C4 * C4 * C4 / Discriminant;

        /// <summary>Create a curve with coefficients in the specified field presentation; reject singular equations.</summary>
        public EllipticCurveFq(FiniteField field, FiniteFieldElement a1, FiniteFieldElement a2,
            FiniteFieldElement a3, FiniteFieldElement a4, FiniteFieldElement a6)
        {
            Field = field ?? throw new ArgumentNullException(nameof(field));
            foreach (var a in new[] { a1, a2, a3, a4, a6 }) field.Validate(a);
            A1 = a1; A2 = a2; A3 = a3; A4 = a4; A6 = a6;
            var b2 = a1 * a1 + 4 * a2; var b4 = 2 * a4 + a1 * a3; var b6 = a3 * a3 + 4 * a6;
            var b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
            C4 = b2 * b2 - 24 * b4; C6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
            Discriminant = -b2 * b2 * b8 - 8 * b4 * b4 * b4 - 27 * b6 * b6 + 9 * b2 * b4 * b6;
            if (Discriminant.IsZero) throw new ArgumentException("The equation is singular over this field.");
        }

        /// <summary>Create a curve by embedding integer coefficients into the prime subfield.</summary>
        public EllipticCurveFq(FiniteField field, BigInteger a1, BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6)
            : this(field, (field ?? throw new ArgumentNullException(nameof(field))).CreateElement(a1), field.CreateElement(a2),
                  field.CreateElement(a3), field.CreateElement(a4), field.CreateElement(a6)) { }

        private FiniteFieldElement Rhs(FiniteFieldElement x) => ((x + A2) * x + A4) * x + A6;
        private void Validate(EllipticCurvePointFq point)
        { if (!IsOnCurve(point)) throw new ArgumentException("Point is not on this curve over its field presentation.", nameof(point)); }

        /// <summary>Create an affine point, validating both the fields and the equation.</summary>
        public EllipticCurvePointFq CreatePoint(FiniteFieldElement x, FiniteFieldElement y)
        { var point = new EllipticCurvePointFq(x, y); Validate(point); return point; }
        /// <summary>Create a point with integer coordinates in the prime subfield.</summary>
        public EllipticCurvePointFq CreatePoint(BigInteger x, BigInteger y) => CreatePoint(Field.CreateElement(x), Field.CreateElement(y));
        /// <summary>Check the field presentations and curve equation. Default points are invalid.</summary>
        public bool IsOnCurve(EllipticCurvePointFq point) => point.IsInfinity ||
            (Field.Equals(point.X.Field) && Field.Equals(point.Y.Field) &&
             point.Y * point.Y + (A1 * point.X + A3) * point.Y == Rhs(point.X));
        /// <summary>Negate a point.</summary>
        public EllipticCurvePointFq Negate(EllipticCurvePointFq point)
        { Validate(point); return point.IsInfinity ? point : new EllipticCurvePointFq(point.X, -point.Y - A1 * point.X - A3); }
        /// <summary>Add two points using the general Weierstrass group law, including characteristics two and three.</summary>
        public EllipticCurvePointFq Add(EllipticCurvePointFq p, EllipticCurvePointFq q)
        {
            Validate(p); Validate(q);
            if (p.IsInfinity) return q;
            if (q.IsInfinity) return p;
            if (p.X == q.X && (p.Y + q.Y + A1 * p.X + A3).IsZero) return EllipticCurvePointFq.Infinity;
            var slope = p.X != q.X ? (q.Y - p.Y) / (q.X - p.X) :
                (3 * p.X * p.X + 2 * A2 * p.X + A4 - A1 * p.Y) / (2 * p.Y + A1 * p.X + A3);
            var x = slope * slope + A1 * slope - A2 - p.X - q.X;
            var y = -(slope + A1) * x - A3 - p.Y + slope * p.X;
            return new EllipticCurvePointFq(x, y);
        }
        /// <summary>Double a point.</summary>
        public EllipticCurvePointFq Double(EllipticCurvePointFq point) => Add(point, point);
        /// <summary>Subtract a point.</summary>
        public EllipticCurvePointFq Subtract(EllipticCurvePointFq p, EllipticCurvePointFq q) => Add(p, Negate(q));
        /// <summary>Multiply by a signed integer, checking cancellation between group operations.</summary>
        public EllipticCurvePointFq Multiply(EllipticCurvePointFq point, BigInteger n, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point);
            if (n.Sign < 0) { point = Negate(point); n = -n; }
            var result = EllipticCurvePointFq.Infinity;
            while (n > 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (!n.IsEven) result = Add(result, point);
                n >>= 1;
                if (n > 0) point = Double(point);
            }
            return result;
        }

        private void CheckSearchLimit(long maxWork, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (Field.Order * Field.Order > maxWork)
                throw new ArithmeticException("Direct point enumeration requires q^2 coordinate-pair checks; the work limit is too small.");
        }
        /// <summary>Count every point including infinity by direct search. maxWork bounds q^2 coordinate pairs,
        /// not field operations or elapsed time. A limit throws instead of returning a partial count.</summary>
        public BigInteger CountPoints(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            BigInteger count = 0;
            foreach (var point in Points(maxWork, cancellationToken)) count++;
            return count;
        }
        /// <summary>Enumerate all points, infinity first, with a prechecked q^2 coordinate-pair limit.
        /// Enumeration must finish to establish completeness.</summary>
        public IEnumerable<EllipticCurvePointFq> Points(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken);
            yield return EllipticCurvePointFq.Infinity;
            foreach (var x in Field.Elements((long)Field.Order, cancellationToken))
            {
                var rhs = Rhs(x); var b = A1 * x + A3;
                foreach (var y in Field.Elements((long)Field.Order, cancellationToken))
                    if (y * y + b * y == rhs) yield return new EllipticCurvePointFq(x, y);
            }
        }
        /// <inheritdoc/>
        public bool Equals(EllipticCurveFq other) => other != null && Field.Equals(other.Field) &&
            A1 == other.A1 && A2 == other.A2 && A3 == other.A3 && A4 == other.A4 && A6 == other.A6;
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is EllipticCurveFq curve && Equals(curve);
        /// <inheritdoc/>
        public override int GetHashCode() => unchecked(((((A1.GetHashCode() * 397 ^ A2.GetHashCode()) * 397 ^ A3.GetHashCode()) * 397 ^ A4.GetHashCode()) * 397 ^ A6.GetHashCode()));
    }
}

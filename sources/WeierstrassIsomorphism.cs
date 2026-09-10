using System;

namespace EllipticCurves
{
    /// <summary>An exact change x = u² x' + r, y = u³ y' + s u² x' + t from Source to Target.</summary>
    public sealed class WeierstrassIsomorphism
    {
        /// <summary>The original model.</summary>
        public EllipticCurveQ Source { get; }
        /// <summary>The transformed model.</summary>
        public EllipticCurveQ Target { get; }
        /// <summary>Nonzero scaling parameter.</summary>
        public BigRational U { get; }
        /// <summary>Horizontal translation.</summary>
        public BigRational R { get; }
        /// <summary>Shear parameter.</summary>
        public BigRational S { get; }
        /// <summary>Vertical translation.</summary>
        public BigRational T { get; }

        internal WeierstrassIsomorphism(EllipticCurveQ source, BigRational u, BigRational r, BigRational s, BigRational t)
        {
            Source = source ?? throw new ArgumentNullException(nameof(source));
            if (source.IsSingular) throw new InvalidOperationException("An elliptic isomorphism requires a nonsingular curve.");
            if (u.IsZero) throw new ArgumentOutOfRangeException(nameof(u));
            U = u; R = r; S = s; T = t;
            Target = new EllipticCurveQ((source.A1 + 2 * s) / u,
                (source.A2 - s * source.A1 + 3 * r - s * s) / BigRational.Pow(u, 2),
                (source.A3 + r * source.A1 + 2 * t) / BigRational.Pow(u, 3),
                (source.A4 - s * source.A3 + 2 * r * source.A2 - (t + r * s) * source.A1 + 3 * r * r - 2 * s * t) / BigRational.Pow(u, 4),
                (source.A6 + r * source.A4 + r * r * source.A2 + r * r * r - t * source.A3 - t * t - r * t * source.A1) / BigRational.Pow(u, 6));
        }

        /// <summary>Map a point on Source to Target.</summary>
        public EllipticCurvePoint Map(EllipticCurvePoint point)
        {
            if (!Source.IsOnCurve(point)) throw new ArgumentException("Point is not on the source curve.", nameof(point));
            if (point.IsInfinity) return point;
            return new EllipticCurvePoint((point.X - R) / BigRational.Pow(U, 2),
                (point.Y - S * (point.X - R) - T) / BigRational.Pow(U, 3));
        }

        /// <summary>Map a point on Target back to Source.</summary>
        public EllipticCurvePoint MapBack(EllipticCurvePoint point)
        {
            if (!Target.IsOnCurve(point)) throw new ArgumentException("Point is not on the target curve.", nameof(point));
            if (point.IsInfinity) return point;
            return new EllipticCurvePoint(U * U * point.X + R, BigRational.Pow(U, 3) * point.Y + S * U * U * point.X + T);
        }

        /// <summary>The inverse coordinate change.</summary>
        public WeierstrassIsomorphism Inverse() => new WeierstrassIsomorphism(Target, 1 / U, -R / (U * U), -S / U, (S * R - T) / BigRational.Pow(U, 3));
    }
}

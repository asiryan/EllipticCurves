using System;

namespace EllipticCurves
{
    /// <summary>A point on a curve over a finite field. Affine points are created by EllipticCurveFq.CreatePoint.</summary>
    public readonly struct EllipticCurvePointFq : IEquatable<EllipticCurvePointFq>
    {
        /// <summary>Affine x-coordinate, unused at infinity.</summary>
        public FiniteFieldElement X { get; }
        /// <summary>Affine y-coordinate, unused at infinity.</summary>
        public FiniteFieldElement Y { get; }
        /// <summary>Field presentation of an affine point; null at infinity or for an invalid default point.</summary>
        public FiniteField Field => IsInfinity ? null : X.Field;
        /// <summary>Whether this is the point at infinity.</summary>
        public bool IsInfinity { get; }
        /// <summary>Universal point at infinity, accepted on every EllipticCurveFq.</summary>
        public static EllipticCurvePointFq Infinity { get; } = new EllipticCurvePointFq(default, default, true);
        internal EllipticCurvePointFq(FiniteFieldElement x, FiniteFieldElement y, bool infinity = false)
        { X = x; Y = y; IsInfinity = infinity; }
        /// <inheritdoc/>
        public bool Equals(EllipticCurvePointFq other) => IsInfinity ? other.IsInfinity : !other.IsInfinity && X == other.X && Y == other.Y;
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is EllipticCurvePointFq point && Equals(point);
        /// <inheritdoc/>
        public override int GetHashCode() => IsInfinity ? 0 : unchecked(X.GetHashCode() * 397 ^ Y.GetHashCode());
        /// <summary>Exact equality, including the field presentation for affine points.</summary>
        public static bool operator ==(EllipticCurvePointFq a, EllipticCurvePointFq b) => a.Equals(b);
        /// <summary>Exact inequality.</summary>
        public static bool operator !=(EllipticCurvePointFq a, EllipticCurvePointFq b) => !a.Equals(b);
        /// <inheritdoc/>
        public override string ToString() => IsInfinity ? "O" : $"({X}, {Y})";
    }
}

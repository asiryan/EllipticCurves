using System;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>A point over a prime field. Create affine points with EllipticCurveFp.CreatePoint.</summary>
    public readonly struct EllipticCurvePointFp : IEquatable<EllipticCurvePointFp>
    {
        /// <summary>Canonical x-coordinate; unused at infinity.</summary>
        public BigInteger X { get; }
        /// <summary>Canonical y-coordinate; unused at infinity.</summary>
        public BigInteger Y { get; }
        /// <summary>Prime modulus of an affine point; zero for the universal point at infinity.</summary>
        public BigInteger Prime { get; }
        /// <summary>Whether this is the point at infinity.</summary>
        public bool IsInfinity { get; }
        /// <summary>The point at infinity, valid on every elliptic curve over a prime field.</summary>
        public static EllipticCurvePointFp Infinity { get; } = new EllipticCurvePointFp(0, 0, 0, true);
        internal EllipticCurvePointFp(BigInteger prime, BigInteger x, BigInteger y, bool infinity = false)
        { Prime = prime; X = x; Y = y; IsInfinity = infinity; }
        /// <inheritdoc/>
        public bool Equals(EllipticCurvePointFp other) => IsInfinity ? other.IsInfinity :
            !other.IsInfinity && Prime == other.Prime && X == other.X && Y == other.Y;
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is EllipticCurvePointFp point && Equals(point);
        /// <inheritdoc/>
        public override int GetHashCode() => IsInfinity ? 0 : unchecked((Prime.GetHashCode() * 397 ^ X.GetHashCode()) * 397 ^ Y.GetHashCode());
        /// <summary>Exact point equality, including the field modulus.</summary>
        public static bool operator ==(EllipticCurvePointFp a, EllipticCurvePointFp b) => a.Equals(b);
        /// <summary>Exact point inequality.</summary>
        public static bool operator !=(EllipticCurvePointFp a, EllipticCurvePointFp b) => !a.Equals(b);
        /// <inheritdoc/>
        public override string ToString() => IsInfinity ? "O" : $"({X}, {Y}) mod {Prime}";
    }
}

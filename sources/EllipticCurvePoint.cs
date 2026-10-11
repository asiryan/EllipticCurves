using System;

namespace EllipticCurves
{
    /// <summary>
    /// Immutable rational point coordinates, or the point at infinity.
    /// 
    /// <para>
    /// Representation:
    /// <list type="bullet">
    ///   <item><description>Finite points are stored as exact rationals <see cref="X"/>, <see cref="Y"/> with <see cref="IsInfinity"/> = <c>false</c>.</description></item>
    ///   <item><description>The neutral element (point at infinity) is represented canonically by <see cref="Infinity"/> with <see cref="IsInfinity"/> = <c>true</c> and <see cref="X"/>=<see cref="Y"/>=<c>0</c>.</description></item>
    /// </list>
    /// </para>
    /// 
    /// <para>
    /// Equality and hashing follow value semantics:
    /// two finite points are equal iff both coordinates are equal as rationals; the
    /// infinity point is only equal to itself.
    /// </para>
    /// <para>Coordinates carry no curve identity and need not satisfy a curve equation.
    /// default(EllipticCurvePoint) is the finite point (0,0); its membership depends on the receiving curve.</para>
    /// </summary>
    public readonly partial struct EllipticCurvePoint : IEquatable<EllipticCurvePoint>
    {
        /// <summary>
        /// Affine <c>x</c>-coordinate (exact rational). For <see cref="Infinity"/>, this is <c>0</c>.
        /// </summary>
        public readonly BigRational X;

        /// <summary>
        /// Affine <c>y</c>-coordinate (exact rational). For <see cref="Infinity"/>, this is <c>0</c>.
        /// </summary>
        public readonly BigRational Y;

        /// <summary>
        /// Indicates whether the point is the neutral element (the point at infinity).
        /// </summary>
        public readonly bool IsInfinity;

        /// <summary>
        /// Create a finite point with given affine coordinates <paramref name="x"/>, <paramref name="y"/>.
        /// This constructor stores coordinates without checking any curve equation.
        /// Use <see cref="EllipticCurveQ.CreatePoint"/> or <see cref="EllipticCurveQ.TryCreatePoint"/>
        /// to validate coordinates against a nonsingular rational curve.
        /// </summary>
        public EllipticCurvePoint(BigRational x, BigRational y) { X = x; Y = y; IsInfinity = false; }

        /// <summary>
        /// Private constructor used to create the canonical infinity element.
        /// </summary>
        private EllipticCurvePoint(bool inf) { X = BigRational.Zero; Y = BigRational.Zero; IsInfinity = inf; }

        /// <summary>
        /// The neutral element of the elliptic curve group (point at infinity).
        /// </summary>
        public static EllipticCurvePoint Infinity => new EllipticCurvePoint(true);

        /// <summary>
        /// Value equality: finite points compare by coordinates; infinity compares only to infinity.
        /// </summary>
        public bool Equals(EllipticCurvePoint other)
        {
            if (IsInfinity || other.IsInfinity) return IsInfinity == other.IsInfinity;
            return X == other.X && Y == other.Y;
        }

        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is EllipticCurvePoint p && Equals(p);

        /// <inheritdoc/>
        public override int GetHashCode()
        {
            unchecked
            {
                int h = 17;
                h = h * 31 + X.GetHashCode();
                h = h * 31 + Y.GetHashCode();
                h = h * 31 + (IsInfinity ? 1 : 0);
                return h;
            }
        }

        /// <summary>
        /// Human-readable representation: <c>"O"</c> for infinity; <c>"(x, y)"</c> for finite points.
        /// Uses <see cref="BigRational.ToString()"/> for coordinates (culture-invariant).
        /// </summary>
        public override string ToString() => IsInfinity ? "O" : $"({X}, {Y})";
    }
}

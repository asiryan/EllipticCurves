using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>A nonsingular general Weierstrass curve over F_p, including characteristics two and three.
    /// Counting and enumeration use direct search, not SEA.</summary>
    /// <remarks>Point operations validate membership and the prime modulus, including before identity shortcuts and a zero multiplier.
    /// Infinity is accepted; a default affine point is invalid. Group operations throw ArgumentException for invalid points.
    /// Token-taking operations observe cooperative cancellation between arithmetic steps; individual BigInteger operations are not interruptible.</remarks>
    public sealed class EllipticCurveFp : IEquatable<EllipticCurveFp>
    {
        /// <summary>Prime defining the base field.</summary>
        public BigInteger Prime { get; }
        /// <summary>Coefficient of xy.</summary>
        public BigInteger A1 { get; }
        /// <summary>Coefficient of x^2.</summary>
        public BigInteger A2 { get; }
        /// <summary>Coefficient of y.</summary>
        public BigInteger A3 { get; }
        /// <summary>Coefficient of x.</summary>
        public BigInteger A4 { get; }
        /// <summary>Constant coefficient.</summary>
        public BigInteger A6 { get; }
        /// <summary>Nonzero discriminant in F_p.</summary>
        public BigInteger Discriminant { get; }
        /// <summary>Invariant c4 in F_p.</summary>
        public BigInteger C4 { get; }
        /// <summary>Invariant c6 in F_p.</summary>
        public BigInteger C6 { get; }
        /// <summary>j-invariant in F_p.</summary>
        public BigInteger JInvariant => Mod(C4 * C4 * C4 * Inverse(Discriminant));

        /// <summary>Create a curve, prove primality of p and reject singular equations.
        /// Primality verification can be expensive for primes beyond 64 bits.</summary>
        /// <exception cref="ArgumentOutOfRangeException">prime is not prime.</exception>
        /// <exception cref="ArgumentException">The equation is singular modulo prime.</exception>
        public EllipticCurveFp(BigInteger prime, BigInteger a1, BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6,
            CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            Prime = prime; A1 = Mod(a1); A2 = Mod(a2); A3 = Mod(a3); A4 = Mod(a4); A6 = Mod(a6);
            var b2 = Mod(A1 * A1 + 4 * A2); var b4 = Mod(2 * A4 + A1 * A3); var b6 = Mod(A3 * A3 + 4 * A6);
            var b8 = Mod(A1 * A1 * A6 + 4 * A2 * A6 - A1 * A3 * A4 + A2 * A3 * A3 - A4 * A4);
            C4 = Mod(b2 * b2 - 24 * b4); C6 = Mod(-b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6);
            Discriminant = Mod(-b2 * b2 * b8 - 8 * b4 * b4 * b4 - 27 * b6 * b6 + 9 * b2 * b4 * b6);
            cancellationToken.ThrowIfCancellationRequested();
            if (Discriminant.IsZero) throw new ArgumentException("The equation is singular over the requested prime field.");
        }

        private BigInteger Mod(BigInteger value) => NativeNumberTheory.Mod(value, Prime);
        private BigInteger ModPow(BigInteger value, BigInteger exponent, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            var result = BigInteger.ModPow(value, exponent, Prime);
            token.ThrowIfCancellationRequested();
            return result;
        }
        private BigInteger Inverse(BigInteger value, CancellationToken token = default)
        {
            token.ThrowIfCancellationRequested();
            value = Mod(value);
            if (value.IsZero) throw new DivideByZeroException();
            return ModPow(value, Prime - 2, token);
        }
        private BigInteger Rhs(BigInteger x, CancellationToken token = default)
        {
            token.ThrowIfCancellationRequested();
            var result = Mod(((x + A2) * x + A4) * x + A6);
            token.ThrowIfCancellationRequested();
            return result;
        }
        private void Validate(EllipticCurvePointFp point, CancellationToken token = default)
        { if (!IsOnCurveCore(point, token)) throw new ArgumentException("Point is not on this curve over its base field.", nameof(point)); }

        /// <summary>Reduce coordinates modulo p and create a point, rejecting points outside the curve.</summary>
        public EllipticCurvePointFp CreatePoint(BigInteger x, BigInteger y)
        { var point = new EllipticCurvePointFp(Prime, Mod(x), Mod(y)); Validate(point); return point; }
        /// <summary>Accept infinity, or check both the field modulus and the curve equation. Default points return false.</summary>
        public bool IsOnCurve(EllipticCurvePointFp point) => IsOnCurveCore(point, default);
        private bool IsOnCurveCore(EllipticCurvePointFp point, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            bool result = point.IsInfinity || (point.Prime == Prime &&
                Mod(point.Y * point.Y + (A1 * point.X + A3) * point.Y - Rhs(point.X, token)).IsZero);
            token.ThrowIfCancellationRequested();
            return result;
        }
        /// <summary>Negate a point.</summary>
        public EllipticCurvePointFp Negate(EllipticCurvePointFp point) => NegateCore(point, default);
        private EllipticCurvePointFp NegateCore(EllipticCurvePointFp point, CancellationToken token)
        {
            Validate(point, token);
            var result = point.IsInfinity ? point : new EllipticCurvePointFp(Prime, point.X, Mod(-point.Y - A1 * point.X - A3));
            token.ThrowIfCancellationRequested();
            return result;
        }
        /// <summary>Add two points using the general Weierstrass group law.</summary>
        public EllipticCurvePointFp Add(EllipticCurvePointFp p, EllipticCurvePointFp q) => AddCore(p, q, default);
        private EllipticCurvePointFp AddCore(EllipticCurvePointFp p, EllipticCurvePointFp q, CancellationToken token)
        {
            Validate(p, token); Validate(q, token);
            if (p.IsInfinity) return q; if (q.IsInfinity) return p;
            bool inversePoints = p.X == q.X && Mod(p.Y + q.Y + A1 * p.X + A3).IsZero;
            token.ThrowIfCancellationRequested();
            if (inversePoints) return EllipticCurvePointFp.Infinity;
            var slope = p.X != q.X ? Mod((q.Y - p.Y) * Inverse(q.X - p.X, token)) :
                Mod((3 * p.X * p.X + 2 * A2 * p.X + A4 - A1 * p.Y) * Inverse(2 * p.Y + A1 * p.X + A3, token));
            var x = Mod(slope * slope + A1 * slope - A2 - p.X - q.X);
            var y = Mod(-(slope + A1) * x - A3 - p.Y + slope * p.X);
            token.ThrowIfCancellationRequested();
            return new EllipticCurvePointFp(Prime, x, y);
        }
        /// <summary>Double a point.</summary>
        public EllipticCurvePointFp Double(EllipticCurvePointFp point) => Add(point, point);
        /// <summary>Subtract a point.</summary>
        public EllipticCurvePointFp Subtract(EllipticCurvePointFp p, EllipticCurvePointFp q) => Add(p, Negate(q));
        /// <summary>Multiply by an arbitrary signed integer.</summary>
        public EllipticCurvePointFp Multiply(EllipticCurvePointFp point, BigInteger n, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point, cancellationToken);
            if (n.Sign < 0) { point = NegateCore(point, cancellationToken); n = -n; }
            var result = EllipticCurvePointFp.Infinity;
            while (n > 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (!n.IsEven) result = AddCore(result, point, cancellationToken);
                n >>= 1; if (n > 0) point = AddCore(point, point, cancellationToken);
            }
            cancellationToken.ThrowIfCancellationRequested();
            return result;
        }

        private void CheckSearchLimit(long maxWork, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (Prime > maxWork) throw new ArithmeticException("The direct search needs at least p x-coordinate steps.");
        }
        /// <summary>Count all points, including infinity. maxWork must be nonnegative and at least p;
        /// it counts x-coordinates, not modular operations or elapsed time. A limit throws before counting starts.</summary>
        /// <exception cref="ArgumentOutOfRangeException">maxWork is negative.</exception>
        /// <exception cref="ArithmeticException">p exceeds maxWork.</exception>
        public BigInteger CountPoints(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken); BigInteger count = 1;
            for (BigInteger x = 0; x < Prime; x++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var rhs = Rhs(x, cancellationToken); var b = Mod(A1 * x + A3);
                if (Prime == 2)
                { for (int y = 0; y < 2; y++) if (Mod(y * y + b * y - rhs).IsZero) count++; }
                else
                {
                    var d = Mod(b * b + 4 * rhs);
                    if (d.IsZero) count++;
                    else if (ModPow(d, (Prime - 1) / 2, cancellationToken).IsOne) count += 2;
                }
            }
            cancellationToken.ThrowIfCancellationRequested();
            return count;
        }
        /// <summary>Lazily enumerate all points, with infinity first. maxWork counts p x-coordinates.
        /// Limits and cancellation are checked on the first MoveNext, before yielding infinity.
        /// The search is complete when enumeration finishes.</summary>
        /// <exception cref="ArgumentOutOfRangeException">maxWork is negative.</exception>
        /// <exception cref="ArithmeticException">p exceeds maxWork.</exception>
        public IEnumerable<EllipticCurvePointFp> Points(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken); yield return EllipticCurvePointFp.Infinity;
            var half = Prime == 2 ? BigInteger.Zero : (Prime + 1) / 2;
            for (BigInteger x = 0; x < Prime; x++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var rhs = Rhs(x, cancellationToken); var b = Mod(A1 * x + A3);
                if (Prime == 2)
                {
                    for (int y = 0; y < 2; y++)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                        if (Mod(y * y + b * y - rhs).IsZero) yield return new EllipticCurvePointFp(Prime, x, y);
                    }
                }
                else if (TrySquareRoot(Mod(b * b + 4 * rhs), cancellationToken, out var root))
                {
                    yield return new EllipticCurvePointFp(Prime, x, Mod((root - b) * half));
                    cancellationToken.ThrowIfCancellationRequested();
                    if (!root.IsZero) yield return new EllipticCurvePointFp(Prime, x, Mod((-root - b) * half));
                }
            }
            cancellationToken.ThrowIfCancellationRequested();
        }
        private bool TrySquareRoot(BigInteger value, CancellationToken token, out BigInteger root)
        {
            token.ThrowIfCancellationRequested();
            root = 0; if (value.IsZero) return true;
            if (ModPow(value, (Prime - 1) / 2, token) != 1) return false;
            if (Prime % 4 == 3) { root = ModPow(value, (Prime + 1) / 4, token); return true; }
            BigInteger odd = Prime - 1; int s = 0;
            while (odd.IsEven) { token.ThrowIfCancellationRequested(); odd >>= 1; s++; }
            BigInteger z = 2;
            while (ModPow(z, (Prime - 1) / 2, token) != Prime - 1) z++;
            var c = ModPow(z, odd, token); root = ModPow(value, (odd + 1) / 2, token);
            var t = ModPow(value, odd, token);
            while (t != 1)
            {
                token.ThrowIfCancellationRequested(); int i = 0; var square = t;
                while (square != 1 && i < s) { token.ThrowIfCancellationRequested(); square = Mod(square * square); i++; }
                if (i == s) throw new ArithmeticException("Modular square-root iteration failed.");
                var b = ModPow(c, BigInteger.One << (s - i - 1), token);
                root = Mod(root * b); c = Mod(b * b); t = Mod(t * c); s = i;
            }
            token.ThrowIfCancellationRequested();
            return true;
        }
        /// <summary>Exact point order. Infinity returns one without counting; maxWork must still be nonnegative.
        /// For affine points, uses a complete direct point count and certified factorization of the group order,
        /// so maxWork must be at least p and has the same meaning as in CountPoints.</summary>
        public BigInteger GetPointOrder(EllipticCurvePointFp point, long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point, cancellationToken);
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (point.IsInfinity) return 1;
            var order = CountPoints(maxWork, cancellationToken);
            foreach (var prime in NativeNumberTheory.Factor(order, cancellationToken).Keys)
                while (order % prime == 0 && Multiply(point, order / prime, cancellationToken).IsInfinity) order /= prime;
            cancellationToken.ThrowIfCancellationRequested();
            return order;
        }
        /// <inheritdoc/>
        public bool Equals(EllipticCurveFp other) => other != null && Prime == other.Prime && A1 == other.A1 && A2 == other.A2 && A3 == other.A3 && A4 == other.A4 && A6 == other.A6;
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is EllipticCurveFp curve && Equals(curve);
        /// <inheritdoc/>
        public override int GetHashCode() => unchecked((((((Prime.GetHashCode() * 397 ^ A1.GetHashCode()) * 397 ^ A2.GetHashCode()) * 397 ^ A3.GetHashCode()) * 397 ^ A4.GetHashCode()) * 397 ^ A6.GetHashCode()));
    }
}

using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>A nonsingular general Weierstrass curve over F_p, including characteristics two and three.
    /// Counting and enumeration use direct search, not SEA.</summary>
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
            if (Discriminant.IsZero) throw new ArgumentException("The equation is singular over the requested prime field.");
        }

        private BigInteger Mod(BigInteger value) => NativeNumberTheory.Mod(value, Prime);
        private BigInteger Inverse(BigInteger value)
        {
            value = Mod(value);
            if (value.IsZero) throw new DivideByZeroException();
            return BigInteger.ModPow(value, Prime - 2, Prime);
        }
        private BigInteger Rhs(BigInteger x) => Mod(((x + A2) * x + A4) * x + A6);
        private void Validate(EllipticCurvePointFp point)
        { if (!IsOnCurve(point)) throw new ArgumentException("Point is not on this curve over its base field.", nameof(point)); }

        /// <summary>Reduce coordinates modulo p and create a point, rejecting points outside the curve.</summary>
        public EllipticCurvePointFp CreatePoint(BigInteger x, BigInteger y)
        { var point = new EllipticCurvePointFp(Prime, Mod(x), Mod(y)); Validate(point); return point; }
        /// <summary>Check both the field modulus and the curve equation.</summary>
        public bool IsOnCurve(EllipticCurvePointFp point) => point.IsInfinity ||
            (point.Prime == Prime && Mod(point.Y * point.Y + (A1 * point.X + A3) * point.Y - Rhs(point.X)).IsZero);
        /// <summary>Negate a point.</summary>
        public EllipticCurvePointFp Negate(EllipticCurvePointFp point)
        { Validate(point); return point.IsInfinity ? point : new EllipticCurvePointFp(Prime, point.X, Mod(-point.Y - A1 * point.X - A3)); }
        /// <summary>Add two points using the general Weierstrass group law.</summary>
        public EllipticCurvePointFp Add(EllipticCurvePointFp p, EllipticCurvePointFp q)
        {
            Validate(p); Validate(q);
            if (p.IsInfinity) return q; if (q.IsInfinity) return p;
            if (p.X == q.X && Mod(p.Y + q.Y + A1 * p.X + A3).IsZero) return EllipticCurvePointFp.Infinity;
            var slope = p.X != q.X ? Mod((q.Y - p.Y) * Inverse(q.X - p.X)) :
                Mod((3 * p.X * p.X + 2 * A2 * p.X + A4 - A1 * p.Y) * Inverse(2 * p.Y + A1 * p.X + A3));
            var x = Mod(slope * slope + A1 * slope - A2 - p.X - q.X);
            var y = Mod(-(slope + A1) * x - A3 - p.Y + slope * p.X);
            return new EllipticCurvePointFp(Prime, x, y);
        }
        /// <summary>Double a point.</summary>
        public EllipticCurvePointFp Double(EllipticCurvePointFp point) => Add(point, point);
        /// <summary>Subtract a point.</summary>
        public EllipticCurvePointFp Subtract(EllipticCurvePointFp p, EllipticCurvePointFp q) => Add(p, Negate(q));
        /// <summary>Multiply by an arbitrary signed integer.</summary>
        public EllipticCurvePointFp Multiply(EllipticCurvePointFp point, BigInteger n, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point);
            if (n.Sign < 0) { point = Negate(point); n = -n; }
            var result = EllipticCurvePointFp.Infinity;
            while (n > 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (!n.IsEven) result = Add(result, point);
                n >>= 1; if (n > 0) point = Double(point);
            }
            return result;
        }

        private void CheckSearchLimit(long maxWork, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (Prime > maxWork) throw new ArithmeticException("The direct search needs at least p x-coordinate steps.");
        }
        /// <summary>Count all points, including infinity. The limit counts x-coordinates, not modular operations.</summary>
        public BigInteger CountPoints(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken); BigInteger count = 1;
            for (BigInteger x = 0; x < Prime; x++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var rhs = Rhs(x); var b = Mod(A1 * x + A3);
                if (Prime == 2)
                { for (int y = 0; y < 2; y++) if (Mod(y * y + b * y - rhs).IsZero) count++; }
                else
                {
                    var d = Mod(b * b + 4 * rhs);
                    if (d.IsZero) count++;
                    else if (BigInteger.ModPow(d, (Prime - 1) / 2, Prime).IsOne) count += 2;
                }
            }
            return count;
        }
        /// <summary>Enumerate all points, with infinity first. The search is complete when enumeration finishes.</summary>
        public IEnumerable<EllipticCurvePointFp> Points(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken); yield return EllipticCurvePointFp.Infinity;
            var half = Prime == 2 ? BigInteger.Zero : (Prime + 1) / 2;
            for (BigInteger x = 0; x < Prime; x++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var rhs = Rhs(x); var b = Mod(A1 * x + A3);
                if (Prime == 2)
                { for (int y = 0; y < 2; y++) if (Mod(y * y + b * y - rhs).IsZero) yield return new EllipticCurvePointFp(Prime, x, y); }
                else if (TrySquareRoot(Mod(b * b + 4 * rhs), cancellationToken, out var root))
                {
                    yield return new EllipticCurvePointFp(Prime, x, Mod((root - b) * half));
                    if (!root.IsZero) yield return new EllipticCurvePointFp(Prime, x, Mod((-root - b) * half));
                }
            }
        }
        private bool TrySquareRoot(BigInteger value, CancellationToken token, out BigInteger root)
        {
            root = 0; if (value.IsZero) return true;
            if (BigInteger.ModPow(value, (Prime - 1) / 2, Prime) != 1) return false;
            if (Prime % 4 == 3) { root = BigInteger.ModPow(value, (Prime + 1) / 4, Prime); return true; }
            BigInteger odd = Prime - 1; int s = 0;
            while (odd.IsEven) { odd >>= 1; s++; }
            BigInteger z = 2;
            while (BigInteger.ModPow(z, (Prime - 1) / 2, Prime) != Prime - 1) { token.ThrowIfCancellationRequested(); z++; }
            var c = BigInteger.ModPow(z, odd, Prime); root = BigInteger.ModPow(value, (odd + 1) / 2, Prime);
            var t = BigInteger.ModPow(value, odd, Prime);
            while (t != 1)
            {
                token.ThrowIfCancellationRequested(); int i = 0; var square = t;
                while (square != 1 && i < s) { square = Mod(square * square); i++; }
                if (i == s) throw new ArithmeticException("Modular square-root iteration failed.");
                var b = BigInteger.ModPow(c, BigInteger.One << (s - i - 1), Prime);
                root = Mod(root * b); c = Mod(b * b); t = Mod(t * c); s = i;
            }
            return true;
        }
        /// <summary>Exact point order, using a complete direct point count and factorization of the group order.</summary>
        public BigInteger GetPointOrder(EllipticCurvePointFp point, long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point);
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (point.IsInfinity) return 1;
            var order = CountPoints(maxWork, cancellationToken);
            foreach (var prime in NativeNumberTheory.Factor(order, cancellationToken).Keys)
                while (order % prime == 0 && Multiply(point, order / prime, cancellationToken).IsInfinity) order /= prime;
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

using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>A nonsingular general Weierstrass curve over a specified finite field F_(p^k).
    /// Supports characteristics two and three. Counting and point enumeration have a q x-coordinate work limit.</summary>
    /// <remarks>Point operations validate membership and field presentations before identity shortcuts or a zero multiplier.
    /// Infinity is accepted; a default affine point is invalid. Group operations throw ArgumentException for invalid points.
    /// Construction and basic group operations have no cancellation token. Multiply, CountPoints and Points check cooperative cancellation;
    /// individual BigInteger operations and field addition are not interruptible.</remarks>
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
        /// <exception cref="ArgumentNullException">field is null.</exception>
        /// <exception cref="ArgumentException">A coefficient is invalid or belongs to an incompatible field, or the equation is singular.</exception>
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

        private FiniteFieldElement Rhs(FiniteFieldElement x, CancellationToken token = default)
        {
            token.ThrowIfCancellationRequested();
            var result = Field.Multiply(Field.Multiply(x + A2, x, token) + A4, x, token) + A6;
            token.ThrowIfCancellationRequested();
            return result;
        }
        private void Validate(EllipticCurvePointFq point, CancellationToken token = default)
        { if (!IsOnCurveCore(point, token)) throw new ArgumentException("Point is not on this curve over its field presentation.", nameof(point)); }

        /// <summary>Create an affine point, validating both the fields and the equation.</summary>
        public EllipticCurvePointFq CreatePoint(FiniteFieldElement x, FiniteFieldElement y)
        { var point = new EllipticCurvePointFq(x, y); Validate(point); return point; }
        /// <summary>Create a point with integer coordinates in the prime subfield.</summary>
        public EllipticCurvePointFq CreatePoint(BigInteger x, BigInteger y) => CreatePoint(Field.CreateElement(x), Field.CreateElement(y));
        /// <summary>Accept infinity, or check the field presentations and curve equation. Default points return false.</summary>
        public bool IsOnCurve(EllipticCurvePointFq point) => IsOnCurveCore(point, default);
        private bool IsOnCurveCore(EllipticCurvePointFq point, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            bool result = point.IsInfinity || (Field.Equals(point.X.Field) && Field.Equals(point.Y.Field) &&
                Field.Multiply(point.Y, point.Y, token) +
                Field.Multiply(Field.Multiply(A1, point.X, token) + A3, point.Y, token) == Rhs(point.X, token));
            token.ThrowIfCancellationRequested();
            return result;
        }
        /// <summary>Negate a point.</summary>
        public EllipticCurvePointFq Negate(EllipticCurvePointFq point) => NegateCore(point, default);
        private EllipticCurvePointFq NegateCore(EllipticCurvePointFq point, CancellationToken token)
        {
            Validate(point, token);
            var result = point.IsInfinity ? point : new EllipticCurvePointFq(point.X, -point.Y - Field.Multiply(A1, point.X, token) - A3);
            token.ThrowIfCancellationRequested();
            return result;
        }
        /// <summary>Add two points using the general Weierstrass group law, including characteristics two and three.</summary>
        public EllipticCurvePointFq Add(EllipticCurvePointFq p, EllipticCurvePointFq q) => AddCore(p, q, default);
        private EllipticCurvePointFq AddCore(EllipticCurvePointFq p, EllipticCurvePointFq q, CancellationToken token)
        {
            Validate(p, token); Validate(q, token);
            if (p.IsInfinity) return q;
            if (q.IsInfinity) return p;
            var a1x = p.X == q.X ? Field.Multiply(A1, p.X, token) : Field.Zero;
            bool inversePoints = p.X == q.X && (p.Y + q.Y + a1x + A3).IsZero;
            token.ThrowIfCancellationRequested();
            if (inversePoints) return EllipticCurvePointFq.Infinity;
            FiniteFieldElement slope;
            if (p.X != q.X) slope = Field.Divide(q.Y - p.Y, q.X - p.X, token);
            else
            {
                var two = Field.CreateElement(2); var three = Field.CreateElement(3);
                var numerator = Field.Multiply(three, Field.Multiply(p.X, p.X, token), token) +
                    Field.Multiply(two, Field.Multiply(A2, p.X, token), token) + A4 - Field.Multiply(A1, p.Y, token);
                var denominator = Field.Multiply(two, p.Y, token) + a1x + A3;
                slope = Field.Divide(numerator, denominator, token);
            }
            var x = Field.Multiply(slope, slope, token) + Field.Multiply(A1, slope, token) - A2 - p.X - q.X;
            var y = -Field.Multiply(slope + A1, x, token) - A3 - p.Y + Field.Multiply(slope, p.X, token);
            token.ThrowIfCancellationRequested();
            return new EllipticCurvePointFq(x, y);
        }
        /// <summary>Double a point.</summary>
        public EllipticCurvePointFq Double(EllipticCurvePointFq point) => Add(point, point);
        /// <summary>Subtract a point.</summary>
        public EllipticCurvePointFq Subtract(EllipticCurvePointFq p, EllipticCurvePointFq q) => Add(p, Negate(q));
        /// <summary>Multiply by a signed integer, checking cancellation within field arithmetic and between group operations.</summary>
        public EllipticCurvePointFq Multiply(EllipticCurvePointFq point, BigInteger n, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(point, cancellationToken);
            if (n.Sign < 0) { point = NegateCore(point, cancellationToken); n = -n; }
            var result = EllipticCurvePointFq.Infinity;
            while (n > 0)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (!n.IsEven) result = AddCore(result, point, cancellationToken);
                n >>= 1;
                if (n > 0) point = AddCore(point, point, cancellationToken);
            }
            cancellationToken.ThrowIfCancellationRequested();
            return result;
        }

        private readonly CachedComputation<PointCountValue> pointCountCache = new CachedComputation<PointCountValue>();
        private sealed class PointCountValue
        {
            internal BigInteger Value { get; }
            internal PointCountValue(BigInteger value) { Value = value; }
        }

        private void CheckSearchLimit(long maxWork, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (maxWork < 0) throw new ArgumentOutOfRangeException(nameof(maxWork));
            if (Field.Order > maxWork)
                throw new ArithmeticException("Point search requires q x-coordinates; the work limit is too small.");
        }
        /// <summary>Count every point including infinity. maxWork must be at least q and bounds x-coordinates,
        /// not field operations or elapsed time. Curves defined over the prime subfield use the exact Frobenius trace recurrence.
        /// Only completed counts are cached; limits and cancellation are checked even for cached or prime-subfield results.</summary>
        /// <exception cref="ArgumentOutOfRangeException">maxWork is negative.</exception>
        /// <exception cref="ArithmeticException">q exceeds maxWork.</exception>
        public BigInteger CountPoints(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken);
            return pointCountCache.Get(() => new PointCountValue(ComputePointCount(cancellationToken)), cancellationToken).Value;
        }

        private BigInteger ComputePointCount(CancellationToken token)
        {
            if (A1.Coordinates.Length <= 1 && A2.Coordinates.Length <= 1 && A3.Coordinates.Length <= 1 &&
                A4.Coordinates.Length <= 1 && A6.Coordinates.Length <= 1)
            {
                BigInteger Constant(FiniteFieldElement value) => value.IsZero ? BigInteger.Zero : value.Coordinates[0];
                var primeCurve = new EllipticCurveFp(Field.Characteristic, Constant(A1), Constant(A2), Constant(A3),
                    Constant(A4), Constant(A6), token);
                var primeTrace = Field.Characteristic + 1 - primeCurve.CountPoints((long)Field.Order, token);
                BigInteger previous = 2, trace = primeTrace;
                for (int i = 2; i <= Field.Degree; i++)
                {
                    token.ThrowIfCancellationRequested();
                    var next = primeTrace * trace - Field.Characteristic * previous;
                    previous = trace; trace = next;
                }
                token.ThrowIfCancellationRequested();
                return Field.Order + 1 - trace;
            }
            var solver = new FiniteFieldQuadraticSolver(Field, token);
            bool characteristicTwo = Field.Characteristic == 2;
            var four = Field.CreateElement(4);
            // Constant b needs only one inversion in characteristic two.
            var constantScale = characteristicTwo && A1.IsZero && !A3.IsZero
                ? solver.InverseSquare(A3, token) : Field.Zero;
            BigInteger count = 1;
            foreach (var x in Field.Elements((long)Field.Order, token))
            {
                var rhs = Rhs(x, token); var b = Field.Multiply(A1, x, token) + A3;
                if (!characteristicTwo)
                {
                    var discriminant = Field.Multiply(b, b, token) + Field.Multiply(four, rhs, token);
                    if (discriminant.IsZero) count++;
                    else if (solver.IsSquare(discriminant, token)) count += 2;
                }
                else if (b.IsZero) count++;
                else
                {
                    var scale = A1.IsZero ? constantScale : solver.InverseSquare(b, token);
                    if (solver.HasArtinSchreierRoot(Field.Multiply(rhs, scale, token), token)) count += 2;
                }
            }
            token.ThrowIfCancellationRequested();
            return count;
        }

        /// <summary>Lazily enumerate all points, infinity first, with a q x-coordinate limit.
        /// The quadratic in y is solved exactly at each x; no q^2 coordinate-pair search is performed.
        /// Limits and cancellation are checked on the first MoveNext, before yielding infinity.
        /// Enumeration must finish to establish completeness. The order of the two ordinates is unspecified.</summary>
        /// <exception cref="ArgumentOutOfRangeException">maxWork is negative.</exception>
        /// <exception cref="ArithmeticException">q exceeds maxWork.</exception>
        public IEnumerable<EllipticCurvePointFq> Points(long maxWork = 1000000, CancellationToken cancellationToken = default)
        {
            CheckSearchLimit(maxWork, cancellationToken);
            yield return EllipticCurvePointFq.Infinity;
            var solver = new FiniteFieldQuadraticSolver(Field, cancellationToken);
            bool characteristicTwo = Field.Characteristic == 2;
            var half = characteristicTwo ? Field.Zero : Field.CreateElement((Field.Characteristic + 1) / 2);
            var four = Field.CreateElement(4);
            var constantScale = characteristicTwo && A1.IsZero && !A3.IsZero
                ? solver.InverseSquare(A3, cancellationToken) : Field.Zero;
            foreach (var x in Field.Elements((long)Field.Order, cancellationToken))
            {
                var rhs = Rhs(x, cancellationToken); var b = Field.Multiply(A1, x, cancellationToken) + A3;
                FiniteFieldElement y;
                if (!characteristicTwo)
                {
                    var discriminant = Field.Multiply(b, b, cancellationToken) + Field.Multiply(four, rhs, cancellationToken);
                    if (!solver.TrySquareRoot(discriminant, out var root, cancellationToken)) continue;
                    y = Field.Multiply(root - b, half, cancellationToken);
                    yield return new EllipticCurvePointFq(x, y);
                    cancellationToken.ThrowIfCancellationRequested();
                    if (!root.IsZero) yield return new EllipticCurvePointFq(x, -y - b);
                }
                else if (b.IsZero)
                {
                    if (!solver.TrySquareRoot(rhs, out y, cancellationToken))
                        throw new ArithmeticException("Squaring must be invertible in characteristic two.");
                    yield return new EllipticCurvePointFq(x, y);
                }
                else
                {
                    var scale = A1.IsZero ? constantScale : solver.InverseSquare(b, cancellationToken);
                    if (!solver.TryArtinSchreierRoot(Field.Multiply(rhs, scale, cancellationToken), out var root, cancellationToken)) continue;
                    y = Field.Multiply(b, root, cancellationToken);
                    yield return new EllipticCurvePointFq(x, y);
                    cancellationToken.ThrowIfCancellationRequested();
                    yield return new EllipticCurvePointFq(x, y + b);
                }
                cancellationToken.ThrowIfCancellationRequested();
            }
            cancellationToken.ThrowIfCancellationRequested();
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

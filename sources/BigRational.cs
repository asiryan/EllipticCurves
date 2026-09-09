using System;
using System.Globalization;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>
    /// A minimal, allocation-friendly rational type for exact arithmetic over Q.
    /// 
    /// Invariants:
    ///  • Num and Den expose canonical form: gcd(|Num|, Den) = 1
    ///  • Denominator is strictly positive (Den > 0)
    ///  • Zero, including default(BigRational), is represented as 0/1
    /// 
    /// Notes:
    ///  • This is an immutable value type with value semantics (Equals/GetHashCode implemented).
    ///  • Designed for exact math in number theory (elliptic curves, invariants, etc.),
    ///    not for floating-point approximations.
    /// </summary>
    public readonly partial struct BigRational : IEquatable<BigRational>, IComparable<BigRational>
    {
        /// <summary>Numerator (can be negative).</summary>
        public BigInteger Num { get; }

        /// <summary>Denominator (strictly positive by invariant).</summary>
        public BigInteger Den => denominator.IsZero ? BigInteger.One : denominator;
        private readonly BigInteger denominator;

        /// <summary>The rational number 0 (stored as 0/1).</summary>
        public static readonly BigRational Zero = new BigRational(BigInteger.Zero, BigInteger.One);

        /// <summary>The rational number 1 (stored as 1/1).</summary>
        public static readonly BigRational One = new BigRational(BigInteger.One, BigInteger.One);

        /// <summary>The rational number 2 (stored as 2/1).</summary>
        public static readonly BigRational Two = new BigRational(2);

        /// <summary>Create from a 64-bit integer (n/1).</summary>
        public BigRational(long n) : this(new BigInteger(n), BigInteger.One) { }

        /// <summary>Create from a BigInteger (n/1).</summary>
        public BigRational(BigInteger n) : this(n, BigInteger.One) { }

        /// <summary>
        /// Create from a numerator/denominator pair and normalize:
        ///  • throw if den == 0
        ///  • move the sign into the numerator
        ///  • divide by gcd(|num|, den)
        ///  • normalize zero to 0/1
        /// </summary>
        public BigRational(BigInteger num, BigInteger den)
        {
            if (den.IsZero) throw new DivideByZeroException("Rational with zero denominator.");

            // Normalize zero as 0/1 immediately.
            if (num.IsZero)
            {
                Num = BigInteger.Zero; denominator = BigInteger.One; return;
            }

            // Keep denominator positive; move sign to numerator.
            if (den.Sign < 0) { num = BigInteger.Negate(num); den = BigInteger.Negate(den); }

            // Reduce by gcd for a canonical representation.
            var g = BigInteger.GreatestCommonDivisor(BigInteger.Abs(num), den);
            Num = num / g; denominator = den / g;
        }

        /// <summary>Create from an integer (n/1).</summary>
        public static BigRational FromInt(long n) => new BigRational(n);

        /// <summary>Create from an unreduced fraction (num/den); will be normalized.</summary>
        public static BigRational FromFraction(BigInteger num, BigInteger den) => new BigRational(num, den);

        /// <summary>True if this is exactly zero.</summary>
        public bool IsZero => Num.IsZero;

        /// <summary>Sign of the rational: −1, 0, or +1.</summary>
        public int Sign => Num.Sign;

        /// <summary>Absolute value (|Num|/Den).</summary>
        public BigRational Abs() => new BigRational(BigInteger.Abs(Num), Den);

        /// <summary>Negation (−Num/Den).</summary>
        public BigRational Negate() => new BigRational(BigInteger.Negate(Num), Den);

        /// <summary>Multiplicative inverse (Den/Num). Throws if zero.</summary>
        public BigRational Reciprocal()
        {
            if (IsZero) throw new DivideByZeroException();
            return new BigRational(Den, Num);
        }

        /// <summary>Exact addition with normalization.</summary>
        public static BigRational operator +(BigRational a, BigRational b)
            => new BigRational(a.Num * b.Den + b.Num * a.Den, a.Den * b.Den);

        /// <summary>Exact subtraction with normalization.</summary>
        public static BigRational operator -(BigRational a, BigRational b)
            => new BigRational(a.Num * b.Den - b.Num * a.Den, a.Den * b.Den);

        /// <summary>Exact multiplication with normalization.</summary>
        public static BigRational operator *(BigRational a, BigRational b)
            => new BigRational(a.Num * b.Num, a.Den * b.Den);

        /// <summary>Exact division with normalization. Throws if divisor is zero.</summary>
        public static BigRational operator /(BigRational a, BigRational b)
        {
            if (b.IsZero) throw new DivideByZeroException();
            return new BigRational(a.Num * b.Den, a.Den * b.Num);
        }

        /// <summary>Unary negation.</summary>
        public static BigRational operator -(BigRational a) => a.Negate();

        /// <summary>
        /// r^k for integer k ≥ 0. For k == 0 returns 1 by convention.
        /// Note: if you need negative powers, use Reciprocal with a positive exponent.
        /// </summary>
        public static BigRational Pow(BigRational r, int k)
        {
            if (k == 0) return BigRational.One; // also handles r == 0
            // BigInteger.Pow requires k >= 0; if k < 0 it will throw — by design here.
            var num = BigInteger.Pow(r.Num, k);
            var den = BigInteger.Pow(r.Den, k);
            return new BigRational(num, den);
        }

        /// <summary>
        /// Exact test for a rational square: returns true iff r = s^2 for some s ∈ Q
        /// and outputs s in canonical form. Works by requiring both numerator and
        /// denominator to be perfect squares in Z.
        /// </summary>
        public static bool IsSquare(BigRational r, out BigRational s)
        {
            s = default;
            if (r.Sign < 0) return false; // no real square roots for negative rationals

            // Check that both numerator and denominator are perfect squares.
            var an = BigInteger.Abs(r.Num);
            var ad = r.Den;

            if (!TryIntegerSquareRoot(an, out var sn)) return false;
            if (!TryIntegerSquareRoot(ad, out var sd)) return false;

            s = new BigRational(sn, sd);
            // Safety: verify s^2 == r exactly (guards against any corner cases).
            return (s * s) == r;
        }

        /// <summary>
        /// Try exact integer square root: returns true iff n is a perfect square and
        /// outputs rt = sqrt(n). Otherwise returns false and sets rt = floor(sqrt(n)).
        /// </summary>
        private static bool TryIntegerSquareRoot(BigInteger n, out BigInteger rt)
        {
            rt = IntegerSquareRoot(n);
            return rt * rt == n;
        }

        /// <summary>
        /// floor(sqrt(n)) for n ≥ 0 using a monotone Newton iteration.
        /// For n ∈ {0,1} returns n; otherwise converges from above.
        /// </summary>
        private static BigInteger IntegerSquareRoot(BigInteger n)
        {
            if (n <= 1) return n;
            BigInteger x0 = n, x1 = (n >> 1) + 1;
            while (x1 < x0) { x0 = x1; x1 = (x1 + n / x1) >> 1; }
            return x0;
        }

        /// <summary>Value equality (uses canonical representation, so structural equality works).</summary>
        public static bool operator ==(BigRational a, BigRational b) => a.Num == b.Num && a.Den == b.Den;

        /// <summary>Value inequality.</summary>
        public static bool operator !=(BigRational a, BigRational b) => !(a == b);

        /// <summary>Strict less-than (exact ℚ order via cross-products).</summary>
        public static bool operator <(BigRational a, BigRational b) => a.CompareTo(b) < 0;

        /// <summary>Strict greater-than (exact order on ℚ; cross-products).</summary>
        public static bool operator >(BigRational a, BigRational b) => a.CompareTo(b) > 0;

        /// <summary>Less-than or equal (exact order on ℚ; cross-products).</summary>
        public static bool operator <=(BigRational a, BigRational b) => a.CompareTo(b) <= 0;

        /// <summary>Greater-than or equal (exact order on ℚ; cross-products).</summary>
        public static bool operator >=(BigRational a, BigRational b) => a.CompareTo(b) >= 0;

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(long value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(ulong value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(short value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(ushort value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(int value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(uint value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(byte value) => new BigRational(value, 1);

        /// <summary>Defines an implicit conversion of a number to BigRational.</summary>
        public static implicit operator BigRational(sbyte value) => new BigRational(value, 1);

        /// <summary>
        /// Lexicographic-free comparison via cross-multiplication:
        /// compare Num/Den and other.Num/other.Den exactly.
        /// </summary>
        public int CompareTo(BigRational other)
        {
            var lhs = Num * other.Den;
            var rhs = other.Num * Den;
            return lhs.CompareTo(rhs);
        }

        /// <summary>Typed equality (same as operator ==).</summary>
        public bool Equals(BigRational other) => this == other;

        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is BigRational r && Equals(r);

        /// <inheritdoc/>
        public override int GetHashCode()
        {
            unchecked
            {
                int h = 17;
                h = h * 31 + Num.GetHashCode();
                h = h * 31 + Den.GetHashCode();
                return h;
            }
        }

        /// <summary>
        /// Culture-invariant string: "n" for integers, "n/d" for nonintegral rationals.
        /// Intended for logs/debugging and round-trippable parsing in simple cases.
        /// </summary>
        public override string ToString()
        {
            if (Den.IsOne) return Num.ToString(CultureInfo.InvariantCulture);
            return $"{Num.ToString(CultureInfo.InvariantCulture)}/{Den.ToString(CultureInfo.InvariantCulture)}";
        }
    }
}

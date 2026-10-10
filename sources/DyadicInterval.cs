using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Outward-rounded dyadic intervals. No floating-point operation enters a proof.
    internal readonly struct DyadicInterval
    {
        internal static readonly BigInteger Scale = BigInteger.One << 160;
        internal readonly BigInteger Lower, Upper;
        internal BigRational LowerRational => new BigRational(Lower, Scale);
        internal BigRational UpperRational => new BigRational(Upper, Scale);
        internal DyadicInterval(BigInteger lower, BigInteger upper)
        {
            if (lower > upper) throw new ArgumentException("Reversed interval.");
            Lower = lower; Upper = upper;
        }
        internal static DyadicInterval Integer(BigInteger n) => new DyadicInterval(n * Scale, n * Scale);
        internal static DyadicInterval Fraction(BigInteger n, BigInteger d) => Integer(n) / Integer(d);
        private static BigInteger Floor(BigInteger n, BigInteger d)
        {
            var q = BigInteger.DivRem(n, d, out var remainder);
            return remainder.Sign < 0 ? q - 1 : q; // d > 0
        }
        private static BigInteger Ceiling(BigInteger n, BigInteger d) => -Floor(-n, d);
        public static DyadicInterval operator +(DyadicInterval a, DyadicInterval b)
            => new DyadicInterval(a.Lower + b.Lower, a.Upper + b.Upper);
        public static DyadicInterval operator -(DyadicInterval a, DyadicInterval b)
            => new DyadicInterval(a.Lower - b.Upper, a.Upper - b.Lower);
        public static DyadicInterval operator *(DyadicInterval a, DyadicInterval b)
        {
            var ll = a.Lower * b.Lower; var lu = a.Lower * b.Upper;
            var ul = a.Upper * b.Lower; var uu = a.Upper * b.Upper;
            return new DyadicInterval(Floor(BigInteger.Min(BigInteger.Min(ll, lu), BigInteger.Min(ul, uu)), Scale),
                Ceiling(BigInteger.Max(BigInteger.Max(ll, lu), BigInteger.Max(ul, uu)), Scale));
        }
        public static DyadicInterval operator /(DyadicInterval a, DyadicInterval b)
        {
            if (b.Lower <= 0 && b.Upper >= 0) throw new DivideByZeroException("Interval contains zero.");
            if (b.Upper < 0) return (Integer(0) - a) / (Integer(0) - b);
            return a * new DyadicInterval(Floor(Scale * Scale, b.Upper), Ceiling(Scale * Scale, b.Lower));
        }
        internal DyadicInterval Widen(BigInteger units) => new DyadicInterval(Lower - units, Upper + units);

        internal static DyadicInterval SqrtInteger(BigInteger n, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            var square = n * Scale * Scale;
            var root = InternalMath.IntegerSqrt(square, cancellationToken);
            return new DyadicInterval(root, root * root == square ? root : root + 1);
        }

        internal static DyadicInterval ExpNegative(DyadicInterval x, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (x.Lower.Sign < 0) throw new ArgumentOutOfRangeException(nameof(x));
            int shifts = 0;
            while (x.Upper > Scale)
            {
                cancellationToken.ThrowIfCancellationRequested();
                x /= Integer(2); shifts++;
            }
            var sum = Integer(1); var term = sum;
            for (int k = 1; ; k++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                term = term * x / Integer(k); sum += term;
                if (term.Upper <= 1)
                {
                    // The remaining positive Taylor terms are bounded by twice the last term.
                    sum = new DyadicInterval(sum.Lower, sum.Upper + 2);
                    break;
                }
            }
            var result = Integer(1) / sum;
            for (int k = 0; k < shifts; k++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                result *= result;
            }
            cancellationToken.ThrowIfCancellationRequested();
            return result;
        }

        internal static DyadicInterval Pi(CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            // Machin's formula and the alternating-series remainder.
            DyadicInterval Atan(int q)
            {
                var power = new BigInteger(q); var sum = Integer(0);
                for (int k = 0; ; k++)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    var term = Fraction(1, (2 * k + 1) * power);
                    sum = k % 2 == 0 ? sum + term : sum - term;
                    if (term.Upper <= 1)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                        return sum.Widen(1);
                    }
                    power *= q * q;
                }
            }
            var result = Integer(16) * Atan(5) - Integer(4) * Atan(239);
            cancellationToken.ThrowIfCancellationRequested();
            return result;
        }
    }
}

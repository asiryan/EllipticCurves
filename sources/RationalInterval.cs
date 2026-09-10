using System;
using System.Linq;
using System.Numerics;

namespace EllipticCurves
{
    // Rational interval endpoints enclose the algebraic numbers used in reduction
    // bounds. Integer search endpoints are always rounded towards the outside.
    internal readonly struct RationalInterval
    {
        internal readonly BigRational Lower, Upper;
        internal RationalInterval(BigRational lower, BigRational upper)
        {
            if (lower > upper) throw new ArgumentException("Reversed interval.");
            Lower = lower; Upper = upper;
        }
        internal static RationalInterval Exact(BigRational n) => new RationalInterval(n, n);
        public static implicit operator RationalInterval(int n) => Exact(n);
        public static implicit operator RationalInterval(BigInteger n) => Exact(new BigRational(n));
        public static RationalInterval operator +(RationalInterval x, RationalInterval y) => new RationalInterval(x.Lower + y.Lower, x.Upper + y.Upper);
        public static RationalInterval operator -(RationalInterval x, RationalInterval y) => new RationalInterval(x.Lower - y.Upper, x.Upper - y.Lower);
        public static RationalInterval operator -(RationalInterval x) => new RationalInterval(-x.Upper, -x.Lower);
        public static RationalInterval operator *(RationalInterval x, RationalInterval y)
        {
            var products = new[] { x.Lower * y.Lower, x.Lower * y.Upper, x.Upper * y.Lower, x.Upper * y.Upper };
            return new RationalInterval(products.Min(), products.Max());
        }
        public static RationalInterval operator /(RationalInterval x, RationalInterval y)
        {
            if (y.Lower <= 0 && y.Upper >= 0) throw new DescentLimitException("Root intervals do not separate a denominator in the quartic bounds.");
            return x * new RationalInterval(1 / y.Upper, 1 / y.Lower);
        }
        internal RationalInterval Square()
        {
            var l = Lower * Lower; var u = Upper * Upper;
            return new RationalInterval(Lower <= 0 && Upper >= 0 ? BigRational.Zero : (l < u ? l : u), l > u ? l : u);
        }
        internal RationalInterval Sqrt()
        {
            if (Upper < 0) throw new InvalidOperationException("Negative radical in quartic reduction bounds.");
            var scale = BigInteger.One << 96;
            BigRational lower = Lower > 0 ? Lower : BigRational.Zero;
            var a = InternalMath.IntegerSqrt(lower.Num * scale * scale / lower.Den);
            var b = InternalMath.IntegerSqrt(Upper.Num * scale * scale / Upper.Den);
            // Widening even an exact root is harmless for enumeration completeness.
            return new RationalInterval(new BigRational(a, scale), new BigRational(b + 1, scale));
        }
    }
}

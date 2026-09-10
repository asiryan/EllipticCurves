using System;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>A real value enclosed by exact rational endpoints.</summary>
    public sealed class RealEnclosure
    {
        /// <summary>Proved lower bound.</summary>
        public BigRational LowerBound { get; }
        /// <summary>Proved upper bound.</summary>
        public BigRational UpperBound { get; }
        /// <summary>Exact width of the enclosure.</summary>
        public BigRational Width => UpperBound - LowerBound;
        /// <summary>Floating-point approximation to the midpoint; use the rational endpoints for proofs.</summary>
        public double Approximation => ToDouble((LowerBound + UpperBound) / 2);
        internal RealEnclosure(BigRational lower, BigRational upper)
        { if (lower > upper) throw new ArgumentException("Reversed enclosure."); LowerBound = lower; UpperBound = upper; }
        /// <summary>Test membership of an exact rational.</summary>
        public bool Contains(BigRational value) => LowerBound <= value && value <= UpperBound;
        /// <summary>Display the approximate midpoint.</summary>
        public override string ToString() => Approximation.ToString("G17", System.Globalization.CultureInfo.InvariantCulture);
        internal static double ToDouble(BigRational r)
        {
            if (r.IsZero) return 0;
            int sign = r.Num.Sign;
            var numerator = BigInteger.Abs(r.Num);
            var denominator = r.Den;
            int exponent = RealArithmetic.BitLength(numerator) - RealArithmetic.BitLength(denominator);
            if (exponent > 1024) return sign * double.PositiveInfinity;
            if (exponent < -1075) return sign * 0.0;
            if ((exponent >= 0 ? numerator.CompareTo(denominator << exponent)
                : (numerator << -exponent).CompareTo(denominator)) < 0) exponent--;
            if (exponent > 1023) return sign * double.PositiveInfinity;

            // Round the exact quotient once, to nearest with ties to even.
            // Subnormals use the fixed 2^-1074 grid; avoid underflowing a scale
            // factor before the significand has been taken into account.
            int power = Math.Max(exponent - 52, -1074);
            if (power >= 0) denominator <<= power;
            else numerator <<= -power;
            var significand = BigInteger.DivRem(numerator, denominator, out var remainder);
            int comparison = (remainder << 1).CompareTo(denominator);
            if (comparison > 0 || (comparison == 0 && !significand.IsEven)) significand++;
            return sign * (double)significand * Math.Pow(2, power);
        }
    }
}

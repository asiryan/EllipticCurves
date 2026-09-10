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
            int n = RealArithmetic.BitLength(r.Num), d = RealArithmetic.BitLength(r.Den);
            int ns = Math.Max(0, n - 54), ds = Math.Max(0, d - 54);
            return (double)(r.Num >> ns) / (double)(r.Den >> ds) * Math.Pow(2, ns - ds);
        }
    }
}

using System;

namespace EllipticCurves
{
    /// <summary>Stored lower and upper bounds for the Mordell–Weil rank.</summary>
    public sealed class LmfdbRankBounds
    {
        /// <summary>The lower bound.</summary>
        public int LowerBound { get; }
        /// <summary>The upper bound.</summary>
        public int UpperBound { get; }
        internal LmfdbRankBounds(int lower, int upper)
        {
            if (lower < 0 || upper < lower) throw new FormatException("LMFDB: invalid rank bounds.");
            LowerBound = lower; UpperBound = upper;
        }
    }
}

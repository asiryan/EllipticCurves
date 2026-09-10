using System;

namespace EllipticCurves
{
    /// <summary>Limits for exact rational division points. Exhaustion throws ArithmeticException.</summary>
    public sealed class PointDivisionOptions
    {
        /// <summary>Maximum counted polynomial and root-isolation steps; excludes torsion preparation.</summary>
        public long MaxWork { get; set; } = 2000000;
        /// <summary>Maximum degree n^2 of a multiplication equation, between 1 and 100000.</summary>
        public int MaxDivisionDegree { get; set; } = 1024;
        internal PointDivisionOptions Snapshot()
        {
            var result = (PointDivisionOptions)MemberwiseClone();
            if (result.MaxWork < 0) throw new ArgumentOutOfRangeException(nameof(MaxWork));
            if (result.MaxDivisionDegree < 1 || result.MaxDivisionDegree > 100000) throw new ArgumentOutOfRangeException(nameof(MaxDivisionDegree));
            return result;
        }
    }
}

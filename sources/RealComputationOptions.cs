using System;

namespace EllipticCurves
{
    /// <summary>Limits for certified real computations. Exceeding a limit throws ArithmeticException.</summary>
    public sealed class RealComputationOptions
    {
        /// <summary>Requested absolute accuracy: enclosure width at most 10^(-DecimalDigits).</summary>
        public int DecimalDigits { get; set; } = 12;
        /// <summary>Bits after the binary point used for outward-rounded arithmetic.</summary>
        public int PrecisionBits { get; set; } = 256;
        /// <summary>Extra decimal places in the height-series truncation bound. Increase with PrecisionBits for large regulators.</summary>
        public int GuardDigits { get; set; } = 4;
        /// <summary>Maximum Tate height-series or AGM iterations. Logarithm and pi series lengths follow PrecisionBits.</summary>
        public int MaxIterations { get; set; } = 512;
        /// <summary>Maximum counted root-isolation steps.</summary>
        public long MaxRootWork { get; set; } = 100000;
        internal RealComputationOptions Snapshot()
        {
            var o = (RealComputationOptions)MemberwiseClone();
            if (o.DecimalDigits < 1 || o.DecimalDigits > 100) throw new ArgumentOutOfRangeException(nameof(DecimalDigits));
            if (o.PrecisionBits < 64 || o.PrecisionBits > 4096) throw new ArgumentOutOfRangeException(nameof(PrecisionBits));
            if (o.GuardDigits < 4 || o.GuardDigits > 100) throw new ArgumentOutOfRangeException(nameof(GuardDigits));
            if (o.MaxIterations < 1 || o.MaxIterations > 100000) throw new ArgumentOutOfRangeException(nameof(MaxIterations));
            if (o.MaxRootWork < 0) throw new ArgumentOutOfRangeException(nameof(MaxRootWork));
            return o;
        }
    }
}

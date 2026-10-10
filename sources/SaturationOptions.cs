namespace EllipticCurves
{
    /// <summary>Limits for exact saturation at explicitly requested primes.</summary>
    public sealed class SaturationOptions
    {
        /// <summary>Maximum counted polynomial operations, root-isolation steps and candidate tests. Default: 2000000.
        /// Excludes height, torsion and minimal-model preparation. Exhaustion returns a partial saturation result.</summary>
        public long MaxWork { get; set; } = 2000000;
        /// <summary>Maximum degree of a division equation (p² for p-saturation), from 4 through 100000. Default: 1024.</summary>
        public int MaxDivisionDegree { get; set; } = 1024;
        /// <summary>Maximum successful subgroup enlargements. Default: 128; zero allows certification without enlargement.</summary>
        public int MaxEnlargements { get; set; } = 128;
        /// <summary>Precision used to certify initial point independence. Defaults to 16 decimal digits and 384 fractional bits.
        /// Null selects the general RealComputationOptions defaults.</summary>
        public RealComputationOptions HeightOptions { get; set; } = new RealComputationOptions { DecimalDigits = 16, PrecisionBits = 384 };
    }
}

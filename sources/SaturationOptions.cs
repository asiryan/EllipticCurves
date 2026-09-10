namespace EllipticCurves
{
    /// <summary>Limits for exact saturation at explicitly requested primes.</summary>
    public sealed class SaturationOptions
    {
        /// <summary>Maximum counted polynomial operations, root-isolation steps and candidate tests.</summary>
        public long MaxWork { get; set; } = 2000000;
        /// <summary>Maximum degree of a division equation (p² for p-saturation).</summary>
        public int MaxDivisionDegree { get; set; } = 1024;
        /// <summary>Maximum successful subgroup enlargements.</summary>
        public int MaxEnlargements { get; set; } = 128;
        /// <summary>Precision used to certify initial point independence.</summary>
        public RealComputationOptions HeightOptions { get; set; } = new RealComputationOptions { DecimalDigits = 16, PrecisionBits = 384 };
    }
}

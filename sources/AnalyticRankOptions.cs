using System;

namespace EllipticCurves
{
    /// <summary>Resource limits and numerical thresholds for native analytic rank estimation.</summary>
    public sealed class AnalyticRankOptions
    {
        /// <summary>Largest derivative order to examine (0 through 8). Default: 4.</summary>
        public int MaxDerivativeOrder { get; set; } = 4;
        /// <summary>Absolute threshold for recognizing numerical zeros (1e-12 through 1e-3). Default: 1e-9.</summary>
        public double ZeroTolerance { get; set; } = 1e-9;
        /// <summary>Maximum number of Fourier coefficients (1 through 1000000). Default: 20000.</summary>
        public int MaxTerms { get; set; } = 20000;
        /// <summary>Maximum sum of primes for direct finite-field point counting. Default: 20000000.</summary>
        public long MaxPointCountingWork { get; set; } = 20000000;
        /// <summary>Maximum integrand evaluations, including the second integration. Default: 100000.</summary>
        public int MaxIntegrationEvaluations { get; set; } = 100000;
        /// <summary>Attempt rigorous certification when the estimated rank is zero or one. Default: true.</summary>
        public bool CertifyLowRanks { get; set; } = true;
        /// <summary>Maximum terms used for the separate interval certificate. Default: 1024; zero disables it.</summary>
        public int MaxCertificationTerms { get; set; } = 1024;

        internal void Validate()
        {
            if (MaxDerivativeOrder < 0 || MaxDerivativeOrder > 8) throw new ArgumentOutOfRangeException(nameof(MaxDerivativeOrder));
            if (!(ZeroTolerance >= 1e-12 && ZeroTolerance <= 1e-3)) throw new ArgumentOutOfRangeException(nameof(ZeroTolerance));
            if (MaxTerms < 1 || MaxTerms > 1000000) throw new ArgumentOutOfRangeException(nameof(MaxTerms));
            if (MaxPointCountingWork < 1) throw new ArgumentOutOfRangeException(nameof(MaxPointCountingWork));
            if (MaxIntegrationEvaluations < 1) throw new ArgumentOutOfRangeException(nameof(MaxIntegrationEvaluations));
            if (MaxCertificationTerms < 0 || MaxCertificationTerms > 1000000) throw new ArgumentOutOfRangeException(nameof(MaxCertificationTerms));
        }
    }
}

using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>Whether an analytic rank was proved, estimated numerically, or left unresolved.</summary>
    public enum AnalyticRankStatus
    {
        /// <summary>No rank is returned: a computation limit or numerical ambiguity was encountered.</summary>
        Inconclusive,
        /// <summary>A numerical estimate; small derivatives have not been proved to vanish.</summary>
        NumericalEstimate,
        /// <summary>Rank zero or one proved by an outward-rounded interval excluding zero.</summary>
        Certified
    }

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

    /// <summary>A native analytic-rank estimate, with its numerical diagnostics and optional rank 0/1 proof.</summary>
    public sealed class AnalyticRankResult
    {
        /// <summary>Distinguishes proofs, numerical estimates, and unresolved computations.</summary>
        public AnalyticRankStatus Status { get; }
        /// <summary>Estimated order of vanishing of L(E,s) at s=1; null if inconclusive.</summary>
        public int? EstimatedRank { get; }
        /// <summary>Proved algebraic and analytic rank (zero or one), or null. Does not assume BSD.</summary>
        public int? ProvenRank => Status == AnalyticRankStatus.Certified ? EstimatedRank : null;
        /// <summary>Whether the rank has been proved.</summary>
        public bool IsCertified => Status == AnalyticRankStatus.Certified;
        /// <summary>The exact conductor used in the functional equation.</summary>
        public BigInteger Conductor { get; }
        /// <summary>The exact global root number, either +1 or -1.</summary>
        public int RootNumber { get; }
        /// <summary>Number of computed Fourier coefficients.</summary>
        public int TermsUsed { get; }
        /// <summary>Estimates of L^(k)(E,1), indexed by k. These are derivatives, not Taylor coefficients.</summary>
        public IReadOnlyList<double> Derivatives { get; }
        /// <summary>Estimated absolute numerical errors, not rigorous interval bounds.</summary>
        public IReadOnlyList<double> EstimatedErrors { get; }
        /// <summary>Rigorous lower endpoint for L(1) or L'(1), when certification was attempted.</summary>
        public BigRational? CertifiedLeadingLowerBound { get; }
        /// <summary>Rigorous upper endpoint for L(1) or L'(1), when certification was attempted.</summary>
        public BigRational? CertifiedLeadingUpperBound { get; }
        /// <summary>Explanation of the result, including any resource limit or unproved assumption.</summary>
        public string Reason { get; }

        internal AnalyticRankResult(AnalyticRankStatus status, int? rank, BigInteger conductor, int rootNumber,
            int terms, double[] derivatives, double[] errors, string reason, BigRational? lower = null, BigRational? upper = null)
        {
            Status = status; EstimatedRank = rank; Conductor = conductor; RootNumber = rootNumber; TermsUsed = terms;
            Derivatives = Array.AsReadOnly(derivatives); EstimatedErrors = Array.AsReadOnly(errors);
            Reason = reason; CertifiedLeadingLowerBound = lower; CertifiedLeadingUpperBound = upper;
        }

        /// <summary>A short description that always includes the evidence status.</summary>
        public override string ToString() => EstimatedRank.HasValue ? $"{EstimatedRank.Value} ({Status})" : $"Inconclusive: {Reason}";
    }
}

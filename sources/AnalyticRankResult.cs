using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
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

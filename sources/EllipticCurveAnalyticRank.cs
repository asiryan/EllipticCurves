using System;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>
        /// Estimate ord_(s=1) L(E,s) natively. BSD predicts that it equals the algebraic rank.
        /// Only an interval certificate for rank 0 or 1 produces ProvenRank; higher estimates are numerical.
        /// Limits apply to Fourier coefficients and integration, not to integer factorization.
        /// </summary>
        public AnalyticRankResult EstimateAnalyticRank(AnalyticRankOptions options = null, CancellationToken cancellationToken = default)
        {
            options = options ?? new AnalyticRankOptions();
            options.Validate();
            cancellationToken.ThrowIfCancellationRequested();
            // Snapshot caller-owned settings before any long computation.
            int order = options.MaxDerivativeOrder, maxTerms = options.MaxTerms;
            int maxEvaluations = options.MaxIntegrationEvaluations, certificateTerms = options.MaxCertificationTerms;
            long maxWork = options.MaxPointCountingWork;
            double tolerance = options.ZeroTolerance;
            bool certify = options.CertifyLowRanks;
            var e = GetGlobalMinimalModel(cancellationToken);
            var conductor = e.GetConductor(cancellationToken);
            int sign = e.GetRootNumber(cancellationToken);
            AnalyticRankResult Unknown(string reason, int terms = 0, double[] values = null, double[] errors = null)
                => new AnalyticRankResult(AnalyticRankStatus.Inconclusive, null, conductor, sign, terms,
                    values ?? Array.Empty<double>(), errors ?? Array.Empty<double>(), reason);
            double alpha = 2 * Math.PI / Math.Sqrt((double)conductor);
            if (!(alpha > 0)) return Unknown("The conductor exceeds the numerical range.");
            double cutoff = 40 + 2 * order;
            while (AnalyticIntegration.TailBound(alpha, cutoff, order) > tolerance * 1e-6) cutoff += 8;
            double needed = Math.Ceiling(cutoff / alpha);
            if (needed > maxTerms) return Unknown("MaxTerms is too small for this conductor and tolerance.");
            int count = Math.Max(1, (int)needed);
            var coefficients = AnalyticCoefficients.Compute(e, count, maxWork, cancellationToken);
            if (coefficients == null) return Unknown("MaxPointCountingWork was exceeded.");
            var integration = new AnalyticIntegration(coefficients, alpha, sign, order, maxEvaluations, cancellationToken);
            double[] completed, completedErrors;
            try
            {
                integration.Compute(cutoff, tolerance, out completed, out completedErrors);
            }
            catch (AnalyticIntegrationLimitException ex)
            {
                return Unknown(ex.Message, count);
            }
            AnalyticIntegration.ToLDerivatives(completed, completedErrors, alpha, out var derivatives, out var errors);
            // Rank detection uses the completed function, whose parity is exact. Ordinary L derivatives
            // do not have this parity, because the gamma factor also has derivatives.
            int? rank = null;
            for (int k = sign == 1 ? 0 : 1; k <= order; k += 2)
            {
                double magnitude = Math.Abs(completed[k]), error = completedErrors[k];
                if (double.IsNaN(magnitude) || double.IsInfinity(magnitude) || error >= tolerance / 4)
                    return Unknown("Numerical errors are too large for ZeroTolerance.", count, derivatives, errors);
                if (magnitude + error < tolerance) continue;
                if (magnitude - error <= 10 * tolerance)
                    return Unknown("A derivative is too close to the numerical zero threshold.", count, derivatives, errors);
                rank = k;
                break;
            }
            if (!rank.HasValue) return Unknown("No nonzero derivative was resolved up to MaxDerivativeOrder.", count, derivatives, errors);

            BigRational? lower = null, upper = null;
            if (certify && rank.Value <= 1 && certificateTerms > 0)
            {
                var interval = AnalyticCertificate.Leading(conductor, coefficients, rank.Value,
                    Math.Min(count, certificateTerms), cancellationToken);
                lower = interval.LowerRational; upper = interval.UpperRational;
                if (interval.Lower.Sign > 0 || interval.Upper.Sign < 0)
                    return new AnalyticRankResult(AnalyticRankStatus.Certified, rank, conductor, sign, count, derivatives, errors,
                        "A rigorous interval excludes zero. The algebraic rank agrees by the rank 0/1 theorems; BSD is not assumed.", lower, upper);
            }
            return new AnalyticRankResult(AnalyticRankStatus.NumericalEstimate, rank, conductor, sign, count, derivatives, errors,
                rank.Value <= 1
                    ? "The rank is a numerical estimate; a rigorous interval certificate was not obtained."
                    : "Small derivatives are numerical zeros only. Identifying the estimate with algebraic rank also uses BSD.", lower, upper);
        }
    }
}

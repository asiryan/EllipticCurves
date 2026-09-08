using System;
using System.Numerics;
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
            catch (AnalyticIntegration.LimitException ex)
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

    internal sealed class AnalyticIntegration
    {
        internal sealed class LimitException : Exception { internal LimitException(string message) : base(message) { } }
        private readonly long[] coefficients;
        private readonly double alpha;
        private readonly int sign, order, maxEvaluations;
        private readonly CancellationToken token;
        private int evaluations;
        private const double Epsilon = 2.2204460492503131e-16;
        internal AnalyticIntegration(long[] coefficients, double alpha, int sign, int order, int maxEvaluations, CancellationToken token)
        {
            this.coefficients = coefficients; this.alpha = alpha; this.sign = sign;
            this.order = order; this.maxEvaluations = maxEvaluations; this.token = token;
        }

        internal static double Factorial(int n) { double f = 1; for (int i = 2; i <= n; i++) f *= i; return f; }

        // |a_n| <= d(n)sqrt(n) <= 2n and log(t)^k <= k! t for t>=1.
        internal static double TailBound(double alpha, double cutoff, int k)
        {
            double q = Math.Exp(-alpha), z = Math.Exp(-cutoff);
            return 4 * Factorial(k) * z * ((1 + 1 / cutoff) / (1 - q) + (cutoff + 1) / (alpha * (1 - z) * (1 - z)));
        }

        internal void Compute(double cutoff, double tolerance, out double[] values, out double[] errors)
        {
            double end = Math.Log(cutoff / alpha);
            var first = Integrate(end, tolerance * 1e-2, 16, out var firstError);
            values = Integrate(end, tolerance * 1e-3, 24, out errors);
            double absoluteSum = 0;
            for (int n = 1; n < coefficients.Length; n++)
                absoluteSum += Math.Abs(coefficients[n]) * Math.Exp(-alpha * n) / n * (1 + 1 / (alpha * n));
            for (int k = 0; k <= order; k++)
            {
                if ((k % 2 == 0 ? 1 : -1) != sign) { values[k] = 0; errors[k] = 0; continue; }
                errors[k] += firstError[k] + Math.Abs(first[k] - values[k]) + TailBound(alpha, cutoff, k)
                    + 64 * Epsilon * Factorial(k) * absoluteSum;
            }
        }

        // F(z)=alpha Lambda(1+z). Its derivatives are moments of the Fourier series on t>=1.
        private double[] Integrand(double u)
        {
            token.ThrowIfCancellationRequested();
            if (evaluations >= maxEvaluations) throw new LimitException("MaxIntegrationEvaluations was exceeded.");
            evaluations++;
            double x = alpha * Math.Exp(u), q = Math.Exp(-x), f = 0;
            // Horner's rule evaluates the finite q-series without separately exponentiating each term.
            for (int n = coefficients.Length - 1; n >= 1; n--)
            {
                if ((n & 4095) == 0) token.ThrowIfCancellationRequested();
                f = (f + coefficients[n]) * q;
            }
            var result = new double[order + 1];
            double moment = 2 * x * f;
            for (int k = 0; k <= order; k++, moment *= u)
                if ((k % 2 == 0 ? 1 : -1) == sign) result[k] = moment;
            return result;
        }

        private double[] Integrate(double end, double accuracy, int pieces, out double[] errors)
        {
            var sum = new double[order + 1]; errors = new double[order + 1];
            for (int i = 0; i < pieces; i++)
            {
                double l = end * i / pieces, r = end * (i + 1) / pieces;
                var fl = Integrand(l); var fm = Integrand((l + r) / 2); var fr = Integrand(r);
                Refine(l, r, fl, fm, fr, Simpson(r - l, fl, fm, fr), accuracy / pieces, 0, sum, errors);
            }
            return sum;
        }

        private double[] Simpson(double width, double[] a, double[] b, double[] c)
        {
            var result = new double[order + 1];
            for (int k = 0; k <= order; k++) result[k] = width * (a[k] + 4 * b[k] + c[k]) / 6;
            return result;
        }
        private void Refine(double l, double r, double[] fl, double[] fm, double[] fr, double[] whole,
            double accuracy, int depth, double[] sum, double[] errors)
        {
            double mid = (l + r) / 2;
            var leftMid = Integrand((l + mid) / 2); var rightMid = Integrand((mid + r) / 2);
            var left = Simpson(mid - l, fl, leftMid, fm); var right = Simpson(r - mid, fm, rightMid, fr);
            bool converged = true;
            for (int k = 0; k <= order; k++)
                if (Math.Abs(left[k] + right[k] - whole[k]) > Math.Max(15 * accuracy,
                    64 * Epsilon * (Math.Abs(left[k]) + Math.Abs(right[k]) + Math.Abs(whole[k])))) converged = false;
            if (converged)
            {
                for (int k = 0; k <= order; k++)
                {
                    double correction = (left[k] + right[k] - whole[k]) / 15;
                    sum[k] += left[k] + right[k] + correction;
                    errors[k] += Math.Abs(correction);
                }
                return;
            }
            if (depth == 24) throw new LimitException("The integration depth limit was exceeded.");
            Refine(l, mid, fl, leftMid, fm, left, accuracy / 2, depth + 1, sum, errors);
            Refine(mid, r, fm, rightMid, fr, right, accuracy / 2, depth + 1, sum, errors);
        }

        internal static void ToLDerivatives(double[] f, double[] fErrors, double alpha, out double[] values, out double[] errors)
        {
            // L(1+z) = F(z) alpha^z / Gamma(1+z). Expand the logarithm, then exponentiate.
            double[] zeta = { 0, 0, 1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
                1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268, 1.0040773561979443394 };
            var log = new double[f.Length]; var gamma = new double[f.Length]; gamma[0] = 1;
            if (f.Length > 1) log[1] = Math.Log(alpha) + 0.57721566490153286061;
            for (int k = 2; k < log.Length; k++) log[k] = (k % 2 == 0 ? -1 : 1) * zeta[k] / k;
            for (int n = 1; n < gamma.Length; n++)
                for (int k = 1; k <= n; k++) gamma[n] += k * log[k] * gamma[n - k] / n;
            values = new double[f.Length]; errors = new double[f.Length];
            for (int n = 0; n < f.Length; n++) for (int k = 0; k <= n; k++)
            {
                double factor = Factorial(n) / Factorial(k) * gamma[n - k];
                values[n] += factor * f[k]; errors[n] += Math.Abs(factor) * fErrors[k];
            }
        }
    }
}

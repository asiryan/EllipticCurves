using System;
using System.Threading;

namespace EllipticCurves
{
    internal sealed class AnalyticIntegration
    {
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
            if (evaluations >= maxEvaluations) throw new AnalyticIntegrationLimitException("MaxIntegrationEvaluations was exceeded.");
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
            if (depth == 24) throw new AnalyticIntegrationLimitException("The integration depth limit was exceeded.");
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

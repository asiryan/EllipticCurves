using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    internal static class AnalyticCertificate
    {
        internal static DyadicInterval Leading(BigInteger conductor, long[] a, int rank, int count, CancellationToken token)
        {
            var alpha = DyadicInterval.Integer(2) * DyadicInterval.Pi() / DyadicInterval.SqrtInteger(conductor);
            var q = DyadicInterval.ExpNegative(alpha); var power = DyadicInterval.Integer(1);
            var sum = DyadicInterval.Integer(0);
            for (int n = 1; n <= count; n++)
            {
                token.ThrowIfCancellationRequested();
                power *= q;
                if (a[n] == 0) continue;
                var weight = rank == 0 ? power : ExponentialIntegral(alpha * DyadicInterval.Integer(n), token);
                sum += DyadicInterval.Fraction(2 * a[n], n) * weight;
            }
            // |a_n| <= 2n. E1(x) <= exp(-x)/x. This bounds every omitted coefficient.
            var tail = DyadicInterval.Integer(4) * power * q / (DyadicInterval.Integer(1) - q);
            if (rank == 1) tail /= alpha * DyadicInterval.Integer(count + 1);
            return sum.Widen(tail.Upper);
        }

        internal static DyadicInterval ExponentialIntegral(DyadicInterval x, CancellationToken token)
        {
            if (x.Lower <= 0) throw new ArgumentOutOfRangeException(nameof(x));
            // E1 is decreasing: evaluate each rational endpoint independently.
            var lower = IntegrateFrom(x.Upper, token);
            var upper = IntegrateFrom(x.Lower, token);
            return new DyadicInterval(lower.Lower, upper.Upper);
        }

        private static DyadicInterval IntegrateFrom(BigInteger start, CancellationToken token)
        {
            var l = new DyadicInterval(start, start); var end = DyadicInterval.Integer(96);
            if (start >= end.Lower)
                return new DyadicInterval(0, (DyadicInterval.ExpNegative(l) / l).Upper);
            var sum = DyadicInterval.Integer(0);
            while (l.Lower < end.Lower)
            {
                token.ThrowIfCancellationRequested();
                var right = BigInteger.Min(2 * l.Lower, end.Lower);
                var r = new DyadicInterval(right, right);
                var c = (l + r) / DyadicInterval.Integer(2);
                var h = (r - l) / DyadicInterval.Integer(2);
                var q = h / c;
                var qPower = DyadicInterval.Integer(1);
                var exponentialTerm = DyadicInterval.Integer(1);
                var partialExponential = exponentialTerm;
                var series = DyadicInterval.Integer(0);
                // exp(-t)/t about c, integrated on [c-h,c+h]. Odd powers integrate to zero.
                // The even coefficients are exp(-c)/c * sum_(k=0)^j c^k/k! / c^j.
                for (int j = 0; j < 80; j++)
                {
                    if (j > 0)
                    {
                        exponentialTerm = exponentialTerm * c / DyadicInterval.Integer(j);
                        partialExponential += exponentialTerm;
                    }
                    if (j % 2 == 0) series += partialExponential * qPower / DyadicInterval.Integer(j + 1);
                    qPower *= q;
                }
                var integral = DyadicInterval.Integer(2) * h * DyadicInterval.ExpNegative(c) / c * series;
                // q<=1/3; S_j(c)<=exp(c) gives a geometric bound for j>=80 (positive terms).
                var remainder = DyadicInterval.Integer(2) * q * qPower
                    / (DyadicInterval.Integer(81) * (DyadicInterval.Integer(1) - q * q));
                sum += new DyadicInterval(integral.Lower, integral.Upper + remainder.Upper);
                l = r;
            }
            return new DyadicInterval(sum.Lower, sum.Upper + (DyadicInterval.ExpNegative(end) / end).Upper);
        }
    }
}

using System.Threading;

namespace EllipticCurves
{
    internal static class AnalyticCoefficients
    {
        // Input is a minimal integral model. a[0] is unused and a[1]=1.
        internal static long[] Compute(EllipticCurveQ e, int count, long maxWork, CancellationToken token)
        {
            var primes = new int[count + 1];
            long work = 0;
            for (int p = 2; p <= count; p++)
            {
                token.ThrowIfCancellationRequested();
                if (primes[p] != 0) continue;
                work += p;
                if (work > maxWork) return null;
                for (int n = p; n <= count; n += p)
                {
                    if ((n & 1023) == 0) token.ThrowIfCancellationRequested();
                    if (primes[n] == 0) primes[n] = p;
                }
            }
            var a = new long[count + 1]; a[1] = 1;
            var delta = e.Discriminant.Num;
            for (int n = 2; n <= count; n++)
            {
                token.ThrowIfCancellationRequested();
                int p = primes[n], m = n / p;
                if (p == n) a[n] = Trace(e, p, token);
                else if (m % p != 0) a[n] = checked(a[p] * a[m]);
                else a[n] = checked(a[p] * a[m] - (delta % p == 0 ? 0 : p * a[m / p]));
            }
            return a;
        }

        internal static long Trace(EllipticCurveQ e, int p, CancellationToken token)
        {
            // Counting the singular cubic also gives the correct bad-prime coefficient: 0 or +/-1.
            return (long)p + 1 - (long)PrimeFieldPointSearch.Count(p, e.A1.Num, e.A2.Num,
                e.A3.Num, e.A4.Num, e.A6.Num, token);
        }
    }
}

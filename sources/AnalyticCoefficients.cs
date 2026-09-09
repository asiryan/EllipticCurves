using System;
using System.Numerics;
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
                for (int n = p; n <= count; n += p) if (primes[n] == 0) primes[n] = p;
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

        private static long Trace(EllipticCurveQ e, int p, CancellationToken token)
        {
            long a1 = (long)NativeNumberTheory.Mod(e.A1.Num, p), a2 = (long)NativeNumberTheory.Mod(e.A2.Num, p);
            long a3 = (long)NativeNumberTheory.Mod(e.A3.Num, p), a4 = (long)NativeNumberTheory.Mod(e.A4.Num, p);
            long a6 = (long)NativeNumberTheory.Mod(e.A6.Num, p);
            if (p == 2)
            {
                int points = 1;
                for (int x = 0; x < 2; x++) for (int y = 0; y < 2; y++)
                    if ((y * y + a1 * x * y + a3 * y - x * x * x - a2 * x * x - a4 * x - a6) % 2 == 0) points++;
                return 3 - points;
            }
            // A quadratic-residue sieve avoids a modular exponentiation for each x.
            var squares = new bool[p];
            for (long y = 1; y <= p / 2; y++) squares[y * y % p] = true;
            long sum = 0;
            for (long x = 0; x < p; x++)
            {
                if ((x & 1023) == 0) token.ThrowIfCancellationRequested();
                long b = (a1 * x + a3) % p;
                long rhs = (((x + a2) * x + a4) % p * x + a6) % p;
                long d = (b * b + 4 * rhs) % p;
                if (d != 0) sum += squares[d] ? 1 : -1;
            }
            // Counting the singular cubic also gives the correct bad-prime coefficient: 0 or +/-1.
            return -sum;
        }
    }
}

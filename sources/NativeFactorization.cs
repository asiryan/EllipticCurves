using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    internal static partial class NativeNumberTheory
    {
        private static BigInteger FindDivisor(BigInteger n, CancellationToken token, int maxWorkers)
        {
            int bits = RealArithmetic.BitLength(n);
            int rhoWork = n <= ulong.MaxValue || bits > 200 ? 131072 : bits < 160 ? 8192 : 32768;
            // Bound each rho attempt: a large least prime factor must not trap the
            // factorizer in rho indefinitely. SIQS handles the remaining cofactors.
            for (int attempt = 1; ; attempt++)
            {
                var divisor = BrentDivisor(n, attempt, rhoWork, token);
                if (divisor > 1) return divisor;
                if (n > ulong.MaxValue && attempt == 2)
                {
                    // For small composites a sieve is already cheaper than an
                    // ECM campaign. Scale preliminary work with the input size.
                    if (bits >= 160)
                    {
                        divisor = NativeEcmFactorization.FindDivisor(n, bits > 200 ? 16 : 8, 500, 5000, token);
                        if (divisor > 1) return divisor;
                    }
                    // Extra ECM work pays off before sending a much larger
                    // cofactor to the sieve. Medium composites are cheaper to sieve.
                    if (bits > 200)
                    {
                        divisor = NativeEcmFactorization.FindDivisor(n, 16, 2000, 20000, token);
                        if (divisor > 1) return divisor;
                    }
                    return NativeQuadraticSieve.FindDivisor(n, token, maxWorkers);
                }
            }
        }

        // Brent's power-of-two cycle search with batched gcds. A batch that
        // contains more than one factor is replayed one difference at a time.
        private static BigInteger BrentDivisor(BigInteger n, int attempt, int limit, CancellationToken token)
        {
            BigInteger y = 2 + attempt, c = attempt, x = 0, saved = 0, divisor = 1;
            int length = 1, steps = 0;
            while (divisor.IsOne && steps < limit)
            {
                x = y;
                for (int i = 0; i < length && steps < limit; i++, steps++)
                {
                    if ((i & 255) == 0) token.ThrowIfCancellationRequested();
                    y = (y * y + c) % n;
                }
                for (int offset = 0; offset < length && divisor.IsOne && steps < limit; offset += 64)
                {
                    token.ThrowIfCancellationRequested();
                    saved = y;
                    BigInteger product = 1;
                    int batch = Math.Min(64, Math.Min(length - offset, limit - steps));
                    for (int i = 0; i < batch; i++, steps++)
                    {
                        y = (y * y + c) % n;
                        product = product * (x - y) % n;
                    }
                    divisor = BigInteger.GreatestCommonDivisor(product, n);
                }
                length *= 2;
            }
            if (divisor == n)
            {
                do
                {
                    token.ThrowIfCancellationRequested();
                    saved = (saved * saved + c) % n;
                    divisor = BigInteger.GreatestCommonDivisor(x - saved, n);
                } while (divisor.IsOne);
            }
            return divisor > 1 && divisor < n ? divisor : BigInteger.One;
        }

        private static bool TryPerfectPower(BigInteger n, CancellationToken token,
            out BigInteger root, out int exponent)
        {
            root = InternalMath.IntegerSqrt(n);
            exponent = 2;
            if (root * root == n) return true;
            int bits = RealArithmetic.BitLength(n);
            // All primes through 37 have already been removed, so a remaining
            // base is at least 41. Composite exponents have a prime divisor.
            foreach (int power in SievePrimes(bits / 5))
            {
                if (power == 2) continue;
                token.ThrowIfCancellationRequested();
                var value = BigInteger.One << ((bits + power - 1) / power);
                while (true)
                {
                    token.ThrowIfCancellationRequested();
                    var next = ((power - 1) * value + n / BigInteger.Pow(value, power - 1)) / power;
                    if (next >= value) break;
                    value = next;
                }
                if (BigInteger.Pow(value, power) == n)
                {
                    root = value;
                    exponent = power;
                    return true;
                }
            }
            return false;
        }

        internal static int[] SievePrimes(int bound)
        {
            var composite = new bool[bound + 1];
            var primes = new List<int>();
            for (int p = 2; p <= bound; p++)
            {
                if (composite[p]) continue;
                primes.Add(p);
                if (p <= bound / p)
                    for (int multiple = p * p; multiple <= bound; multiple += p)
                        composite[multiple] = true;
            }
            return primes.ToArray();
        }
    }
}

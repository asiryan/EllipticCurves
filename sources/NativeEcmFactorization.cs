using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Lenstra ECM on Montgomery curves, with Suyama's parametrization and a
    // baby-step/giant-step second stage. It supplies divisors, never prime claims.
    internal static class NativeEcmFactorization
    {
        private readonly struct Point
        {
            internal readonly BigInteger X, Z;
            internal Point(BigInteger x, BigInteger z) { X = x; Z = z; }
        }

        internal static BigInteger FindDivisor(BigInteger n, int curves, int firstBound,
            int secondBound, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (n <= 1 || curves < 0 || firstBound < 2 || secondBound < firstBound)
                throw new ArgumentOutOfRangeException();
            if (n.IsEven) return n > 2 ? 2 : BigInteger.One;
            var primes = NativeNumberTheory.SievePrimes(secondBound);
            for (int curve = 0; curve < curves; curve++)
            {
                token.ThrowIfCancellationRequested();
                BigInteger sigma = 6 + curve;
                var u = NativeNumberTheory.Mod(sigma * sigma - 5, n);
                var v = 4 * sigma % n;
                var u3 = u * u % n * u % n;
                var denominator = 16 * u3 * v % n;
                var divisor = BigInteger.GreatestCommonDivisor(denominator, n);
                if (divisor > 1 && divisor < n) return divisor;
                if (divisor == n) continue;
                var a24 = NativeNumberTheory.Mod(BigInteger.Pow(v - u, 3) * (3 * u + v)
                    * Inverse(denominator, n), n);
                // Exclude curves singular modulo a factor; a nontrivial gcd is
                // already a factorization success.
                divisor = BigInteger.GreatestCommonDivisor(a24 * (a24 - 1), n);
                if (divisor > 1 && divisor < n) return divisor;
                if (divisor == n) continue;
                var point = new Point(u3, v * v % n * v % n);
                bool degenerate = false;
                foreach (int prime in primes)
                {
                    if (prime > firstBound) break;
                    token.ThrowIfCancellationRequested();
                    int power = prime;
                    while (power <= firstBound / prime) power *= prime;
                    point = Multiply(point, power, a24, n);
                    divisor = BigInteger.GreatestCommonDivisor(point.Z, n);
                    if (divisor > 1 && divisor < n) return divisor;
                    if (divisor == n) { degenerate = true; break; }
                }
                if (degenerate || secondBound == firstBound) continue;
                divisor = SecondStage(point, a24, n, primes, firstBound, token);
                if (divisor > 1) return divisor;
            }
            return BigInteger.One;
        }

        private static BigInteger SecondStage(Point point, BigInteger a24, BigInteger n,
            int[] primes, int firstBound, CancellationToken token)
        {
            // p = m*210 +/- r. Equality of x-coordinates for [m*210]Q and
            // [r]Q detects either sign, using projective cross-products only.
            const int stride = 210;
            var babies = new Point[stride / 2 + 1];
            var twice = Double(point, a24, n);
            babies[1] = point;
            babies[3] = Add(twice, point, point, n);
            for (int r = 5; r < babies.Length; r += 2)
                babies[r] = Add(babies[r - 2], twice, babies[r - 4], n);
            var step = Multiply(point, stride, a24, n);
            int giantIndex = -1;
            Point giant = default, previous = default;
            var terms = new BigInteger[64];
            BigInteger product = 1;
            int count = 0;
            foreach (int prime in primes)
            {
                if (prime <= firstBound) continue;
                token.ThrowIfCancellationRequested();
                // The small primes dividing the wheel use direct multiplication.
                if (prime <= 7)
                {
                    var smallDivisor = BigInteger.GreatestCommonDivisor(Multiply(point, prime, a24, n).Z, n);
                    if (smallDivisor > 1 && smallDivisor < n) return smallDivisor;
                    continue;
                }
                int m = (prime + stride / 2) / stride;
                int r = Math.Abs(prime - m * stride);
                if (giantIndex < 0)
                {
                    giantIndex = m;
                    giant = Multiply(step, m, a24, n);
                    previous = Multiply(step, Math.Max(0, m - 1), a24, n);
                }
                while (giantIndex < m)
                {
                    var next = giantIndex == 0 ? step : giantIndex == 1
                        ? Double(step, a24, n) : Add(giant, step, previous, n);
                    previous = giant; giant = next; giantIndex++;
                }
                var baby = babies[r];
                terms[count] = (giant.X * baby.Z - baby.X * giant.Z) % n;
                product = product * terms[count++] % n;
                if (count == terms.Length)
                {
                    var divisor = BatchDivisor(product, terms, count, n);
                    if (divisor > 1 && divisor < n) return divisor;
                    if (divisor == n) return BigInteger.One;
                    product = 1; count = 0;
                }
            }
            var final = BatchDivisor(product, terms, count, n);
            return final > 1 && final < n ? final : BigInteger.One;
        }

        private static BigInteger BatchDivisor(BigInteger product, BigInteger[] terms, int count, BigInteger n)
        {
            var divisor = BigInteger.GreatestCommonDivisor(product, n);
            if (divisor != n) return divisor;
            for (int i = 0; i < count; i++)
            {
                divisor = BigInteger.GreatestCommonDivisor(terms[i], n);
                if (divisor > 1 && divisor < n) return divisor;
            }
            return n;
        }

        private static Point Multiply(Point point, int scalar, BigInteger a24, BigInteger n)
        {
            if (scalar == 0) return new Point(1, 0);
            if (scalar == 1) return point;
            var lo = point;
            var hi = Double(point, a24, n);
            int bit = 1;
            while (bit <= scalar / 2) bit <<= 1;
            for (bit >>= 1; bit != 0; bit >>= 1)
            {
                var sum = Add(lo, hi, point, n);
                if ((scalar & bit) == 0) { lo = Double(lo, a24, n); hi = sum; }
                else { hi = Double(hi, a24, n); lo = sum; }
            }
            return lo;
        }

        private static Point Double(Point point, BigInteger a24, BigInteger n)
        {
            var sum = point.X + point.Z;
            var difference = point.X - point.Z;
            var aa = sum * sum % n;
            var bb = difference * difference % n;
            var delta = aa - bb;
            return new Point(aa * bb % n, delta * (bb + a24 * delta % n) % n);
        }

        private static Point Add(Point left, Point right, Point difference, BigInteger n)
        {
            var first = (left.X + left.Z) * (right.X - right.Z) % n;
            var second = (left.X - left.Z) * (right.X + right.Z) % n;
            var sum = first + second;
            var delta = first - second;
            return new Point(difference.Z * (sum * sum % n) % n, difference.X * (delta * delta % n) % n);
        }

        private static BigInteger Inverse(BigInteger value, BigInteger modulus)
        {
            BigInteger r = modulus, next = value, t = 0, coefficient = 1;
            while (!next.IsZero)
            {
                var quotient = BigInteger.DivRem(r, next, out var remainder);
                r = next; next = remainder;
                var previous = t;
                t = coefficient; coefficient = previous - quotient * coefficient;
            }
            return NativeNumberTheory.Mod(t, modulus);
        }
    }
}

using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Exact direct searches. The discriminant of the quadratic in y is
    // 4*x^3 + b2*x^2 + 2*b4*x + b6; finite differences evaluate it using
    // three modular additions per x. Tables have fixed caps, independent of p.
    // Algorithmic reference: PARI Fl_elltrace_naive (FpE.c), character tables
    // and forward differences: https://pari.math.u-bordeaux.fr/lcov-report/basemath/FpE.c.gcov.html
    // This is an independent implementation for the general Weierstrass equation.
    internal static class PrimeFieldPointSearch
    {
        internal const int MaxResidueTablePrime = 1 << 24; // At most 2 MiB of bits.
        internal const int MaxRootTablePrime = 1 << 20;    // At most 4 MiB of roots.

        internal static BigInteger Count(BigInteger prime, BigInteger a1, BigInteger a2,
            BigInteger a3, BigInteger a4, BigInteger a6, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (prime == 2)
            {
                int count = 1;
                a1 = Mod(a1, prime); a2 = Mod(a2, prime); a3 = Mod(a3, prime);
                a4 = Mod(a4, prime); a6 = Mod(a6, prime);
                for (int x = 0; x < 2; x++) for (int y = 0; y < 2; y++)
                    if (((y * y + a1 * x * y + a3 * y - x * x * x - a2 * x * x - a4 * x - a6) % 2).IsZero) count++;
                token.ThrowIfCancellationRequested();
                return count;
            }

            if (prime <= MaxResidueTablePrime)
            {
                int p = (int)prime;
                var residues = new uint[(p + 31) / 32];
                long square = 0, increment = 1;
                for (int y = 0; y <= p / 2; y++)
                {
                    if ((y & 1023) == 0) token.ThrowIfCancellationRequested();
                    residues[(int)square >> 5] |= 1u << ((int)square & 31);
                    square = Add(square, increment, p);
                    increment = Add(increment, 2, p);
                }
                var values = new CubicDifferences(p, a1, a2, a3, a4, a6);
                long count = 1;
                for (int x = 0; x < p; x++)
                {
                    if ((x & 1023) == 0) token.ThrowIfCancellationRequested();
                    int d = values.Value;
                    if (d == 0) count++;
                    else if ((residues[d >> 5] & (1u << (d & 31))) != 0) count += 2;
                    values.Advance();
                }
                token.ThrowIfCancellationRequested();
                return count;
            }

            if (prime <= int.MaxValue)
                return CountWordWithoutTable((int)prime, a1, a2, a3, a4, a6, token);
            return CountWithoutTable(prime, a1, a2, a3, a4, a6, token);
        }

        // Proven odd primes fitting Int32 need no BigInteger operations in the
        // search loop. Values and differences are reduced; their sums fit Int64.
        // The count can exceed Int32 even though the modulus fits it.
        internal static BigInteger CountWordWithoutTable(int prime, BigInteger a1, BigInteger a2,
            BigInteger a3, BigInteger a4, BigInteger a6, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            var values = new CubicDifferences(prime, a1, a2, a3, a4, a6);
            long total = (long)prime + 1;
            for (int x = 0; x < prime; x++)
            {
                if ((x & 1023) == 0) token.ThrowIfCancellationRequested();
                total += JacobiWord(values.Value, prime);
                values.Advance();
            }
            token.ThrowIfCancellationRequested();
            return total;
        }

        // Keep the large-modulus path bounded in memory. Jacobi equals
        // Legendre because callers have already proved the odd modulus prime.
        internal static BigInteger CountWithoutTable(BigInteger prime, BigInteger a1, BigInteger a2,
            BigInteger a3, BigInteger a4, BigInteger a6, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            a1 = Mod(a1, prime); a2 = Mod(a2, prime); a3 = Mod(a3, prime);
            a4 = Mod(a4, prime); a6 = Mod(a6, prime);
            var c2 = Mod(a1 * a1 + 4 * a2, prime);
            var c1 = Mod(2 * a1 * a3 + 4 * a4, prime);
            var value = Mod(a3 * a3 + 4 * a6, prime);
            var first = Mod(4 + c2 + c1, prime);
            var second = Mod(24 + 2 * c2, prime);
            var third = Mod(24, prime);
            BigInteger total = prime + 1;
            for (BigInteger x = 0; x < prime; x++)
            {
                token.ThrowIfCancellationRequested();
                total += Jacobi(value, prime, token);
                value = Add(value, first, prime);
                first = Add(first, second, prime);
                second = Add(second, third, prime);
            }
            token.ThrowIfCancellationRequested();
            return total;
        }

        // Encoded root r+1; zero denotes a nonsquare. Only odd primes are used.
        internal static int[] RootTable(int prime, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            var roots = new int[prime];
            long square = 0, increment = 1;
            for (int y = 0; y <= prime / 2; y++)
            {
                if ((y & 1023) == 0) token.ThrowIfCancellationRequested();
                roots[(int)square] = y + 1;
                square = Add(square, increment, prime);
                increment = Add(increment, 2, prime);
            }
            token.ThrowIfCancellationRequested();
            return roots;
        }

        internal struct CubicDifferences
        {
            private readonly long prime;
            private long value, first, second;
            private readonly long third;
            internal int Value => (int)value;

            internal CubicDifferences(int prime, BigInteger a1, BigInteger a2,
                BigInteger a3, BigInteger a4, BigInteger a6)
            {
                this.prime = prime;
                a1 = Mod(a1, prime); a2 = Mod(a2, prime); a3 = Mod(a3, prime);
                a4 = Mod(a4, prime); a6 = Mod(a6, prime);
                var c2 = Mod(a1 * a1 + 4 * a2, prime);
                var c1 = Mod(2 * a1 * a3 + 4 * a4, prime);
                value = (long)Mod(a3 * a3 + 4 * a6, prime);
                first = (long)Mod(4 + c2 + c1, prime);
                second = (long)Mod(24 + 2 * c2, prime);
                third = 24 % prime;
            }

            internal void Advance()
            {
                value = Add(value, first, prime);
                first = Add(first, second, prime);
                second = Add(second, third, prime);
            }
        }

        // Arguments remain nonnegative Int32 values; shifts and remainders
        // decrease them, so neither the reciprocity steps nor stripping powers
        // of two can overflow. Each call has a bounded number of word steps.
        private static int JacobiWord(int a, int n)
        {
            int sign = 1;
            while (a != 0)
            {
                while ((a & 1) == 0)
                {
                    a >>= 1;
                    int residue = n & 7;
                    if (residue == 3 || residue == 5) sign = -sign;
                }
                int swap = a; a = n; n = swap;
                if ((a & 3) == 3 && (n & 3) == 3) sign = -sign;
                a %= n;
            }
            return n == 1 ? sign : 0;
        }

        private static int Jacobi(BigInteger a, BigInteger n, CancellationToken token)
        {
            int sign = 1;
            while (!a.IsZero)
            {
                token.ThrowIfCancellationRequested();
                while (a.IsEven)
                {
                    token.ThrowIfCancellationRequested();
                    a >>= 1;
                    int residue = (int)(n & 7);
                    if (residue == 3 || residue == 5) sign = -sign;
                }
                var swap = a; a = n; n = swap;
                if ((a & 3) == 3 && (n & 3) == 3) sign = -sign;
                a %= n;
            }
            return n.IsOne ? sign : 0;
        }

        private static long Add(long a, long b, long prime)
        { long sum = a + b; return sum >= prime ? sum - prime : sum; }
        private static BigInteger Add(BigInteger a, BigInteger b, BigInteger prime)
        { var sum = a + b; return sum >= prime ? sum - prime : sum; }
        private static BigInteger Mod(BigInteger a, BigInteger prime) => NativeNumberTheory.Mod(a, prime);
    }
}

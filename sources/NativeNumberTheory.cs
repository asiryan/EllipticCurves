using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Exact arithmetic for certified results. In particular, the older FactorAbs
    // probable-prime shortcut must not be used to certify rank or conductor.
    internal static class NativeNumberTheory
    {
        internal static BigInteger Mod(BigInteger a, BigInteger m) => (a % m + m) % m;

        internal static BigInteger Divide(BigInteger a, BigInteger b)
        {
            var q = BigInteger.DivRem(a, b, out var r);
            if (!r.IsZero) throw new InvalidOperationException("Non-integral arithmetic in native curve computation.");
            return q;
        }

        internal static Dictionary<BigInteger, int> Factor(BigInteger n, CancellationToken token)
        {
            var result = new Dictionary<BigInteger, int>();
            n = BigInteger.Abs(n);
            if (n.IsZero) throw new ArgumentException("Cannot factor zero.", nameof(n));
            foreach (int p in new[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 })
            {
                token.ThrowIfCancellationRequested();
                while (n % p == 0)
                {
                    AddFactor(result, p);
                    n /= p;
                }
            }
            FactorRecursive(n, result, token);
            return result;
        }

        private static void AddFactor(Dictionary<BigInteger, int> result, BigInteger p)
        {
            result.TryGetValue(p, out int count);
            result[p] = count + 1;
        }

        private static void FactorRecursive(BigInteger n, Dictionary<BigInteger, int> result, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (n.IsOne) return;
            if (IsPrime(n, token)) { AddFactor(result, n); return; }
            for (BigInteger c = 1; ; c++)
            {
                BigInteger x = 2, y = 2, d = 1;
                while (d.IsOne)
                {
                    token.ThrowIfCancellationRequested();
                    x = (x * x + c) % n;
                    y = (y * y + c) % n;
                    y = (y * y + c) % n;
                    d = BigInteger.GreatestCommonDivisor(BigInteger.Abs(x - y), n);
                }
                if (d == n) continue;
                FactorRecursive(d, result, token);
                FactorRecursive(n / d, result, token);
                return;
            }
        }

        private static bool IsPrime(BigInteger n, CancellationToken token)
        {
            if (n < 2) return false;
            foreach (int p in new[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 })
                if (n % p == 0) return n == p;

            // Deterministic Miller-Rabin on the entire unsigned 64-bit range.
            BigInteger oddPart = n - 1;
            int twos = 0;
            while (oddPart.IsEven) { oddPart >>= 1; twos++; }
            foreach (int basis in new[] { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 })
            {
                token.ThrowIfCancellationRequested();
                var a = basis % n;
                if (a.IsZero) continue;
                var x = BigInteger.ModPow(a, oddPart, n);
                if (x.IsOne || x == n - 1) continue;
                bool passed = false;
                for (int j = 1; j < twos; j++)
                {
                    x = x * x % n;
                    if (x == n - 1) { passed = true; break; }
                }
                if (!passed) return false;
            }
            if (n <= ulong.MaxValue) return true;

            // Above 64 bits, prove primality by the full n-1 (Lucas) criterion.
            // The factors are recursively proved; probable primality is never enough.
            foreach (var q in Factor(n - 1, token).Keys)
            {
                bool proved = false;
                for (int a = 2; a <= 128; a++)
                {
                    token.ThrowIfCancellationRequested();
                    if (BigInteger.ModPow(a, n - 1, n) != 1) return false;
                    if (BigInteger.GreatestCommonDivisor(BigInteger.ModPow(a, (n - 1) / q, n) - 1, n).IsOne)
                    { proved = true; break; }
                }
                if (!proved)
                {
                    // Rare deterministic fallback; cancellable, but potentially slow.
                    for (BigInteger d = 41; d * d <= n; d += 2)
                    {
                        token.ThrowIfCancellationRequested();
                        if (n % d == 0) return false;
                    }
                    return true;
                }
            }
            return true;
        }

        internal static IEnumerable<BigInteger> Divisors(Dictionary<BigInteger, int> factors)
        {
            var values = new List<BigInteger> { BigInteger.One };
            foreach (var pair in factors)
            {
                int count = values.Count;
                BigInteger power = 1;
                for (int e = 1; e <= pair.Value; e++)
                {
                    power *= pair.Key;
                    for (int i = 0; i < count; i++) values.Add(values[i] * power);
                }
            }
            return values;
        }
    }
}

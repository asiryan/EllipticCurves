using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    // Dense ascending coefficients over F_p. All returned arrays are canonical and owned by the caller.
    internal static class FiniteFieldPolynomial
    {
        internal static BigInteger[] Normalize(IReadOnlyList<BigInteger> a, BigInteger p, CancellationToken token)
        {
            var result = new BigInteger[a.Count];
            for (int i = 0; i < result.Length; i++) { token.ThrowIfCancellationRequested(); result[i] = Mod(a[i], p); }
            return Trim(result);
        }

        private static BigInteger[] Trim(BigInteger[] a)
        {
            int length = a.Length;
            while (length > 0 && a[length - 1].IsZero) length--;
            if (length != a.Length) Array.Resize(ref a, length);
            return a;
        }

        internal static BigInteger[] Add(BigInteger[] a, BigInteger[] b, BigInteger p, bool subtract = false)
        {
            var result = new BigInteger[Math.Max(a.Length, b.Length)];
            for (int i = 0; i < result.Length; i++)
                result[i] = Mod((i < a.Length ? a[i] : 0) + (subtract ? -1 : 1) * (i < b.Length ? b[i] : 0), p);
            return Trim(result);
        }

        internal static BigInteger[] Remainder(BigInteger[] a, BigInteger[] b, BigInteger p, CancellationToken token)
        {
            if (b.Length == 0) throw new DivideByZeroException();
            a = (BigInteger[])a.Clone();
            var inverse = b[b.Length - 1].IsOne ? BigInteger.One : BigInteger.ModPow(b[b.Length - 1], p - 2, p);
            for (int i = a.Length - 1; i >= b.Length - 1; i--)
            {
                token.ThrowIfCancellationRequested();
                var factor = Mod(a[i] * inverse, p);
                if (factor.IsZero) continue;
                for (int j = 0; j < b.Length; j++)
                { token.ThrowIfCancellationRequested(); a[i - b.Length + 1 + j] = Mod(a[i - b.Length + 1 + j] - factor * b[j], p); }
            }
            return Trim(a);
        }

        internal static BigInteger[] Multiply(BigInteger[] a, BigInteger[] b, BigInteger[] modulus, BigInteger p, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (a.Length == 0 || b.Length == 0) return Array.Empty<BigInteger>();
            var result = new BigInteger[checked(a.Length + b.Length - 1)];
            for (int i = 0; i < a.Length; i++) for (int j = 0; j < b.Length; j++)
            { token.ThrowIfCancellationRequested(); result[i + j] = Mod(result[i + j] + a[i] * b[j], p); }
            return Remainder(result, modulus, p, token);
        }

        internal static BigInteger[] Power(BigInteger[] a, BigInteger exponent, BigInteger[] modulus, BigInteger p, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            var result = new BigInteger[] { 1 };
            while (exponent > 0)
            {
                token.ThrowIfCancellationRequested();
                if (!exponent.IsEven) result = Multiply(result, a, modulus, p, token);
                exponent >>= 1;
                if (exponent > 0) a = Multiply(a, a, modulus, p, token);
            }
            return result;
        }

        internal static bool IsIrreducible(BigInteger[] f, BigInteger p, CancellationToken token)
        {
            int degree = f.Length - 1;
            if (degree == 1) return true;
            var checkpoints = new HashSet<int>(); int remaining = degree;
            for (int divisor = 2; (long)divisor * divisor <= remaining; divisor++)
            {
                token.ThrowIfCancellationRequested();
                if (remaining % divisor != 0) continue;
                checkpoints.Add(degree / divisor);
                do { remaining /= divisor; } while (remaining % divisor == 0);
            }
            if (remaining > 1) checkpoints.Add(degree / remaining);
            var x = new BigInteger[] { 0, 1 }; var h = x;
            for (int i = 1; i <= degree; i++)
            {
                h = Power(h, p, f, p, token);
                if (!checkpoints.Contains(i)) continue;
                var a = f; var b = Add(h, x, p, true);
                while (b.Length > 0)
                { token.ThrowIfCancellationRequested(); var r = Remainder(a, b, p, token); a = b; b = r; }
                if (a.Length != 1) return false;
            }
            return h.SequenceEqual(x);
        }
    }
}

using System;
using System.Linq;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    // Small-degree arithmetic over F_p; gcd(f, X^p-X) counts distinct rational roots.
    internal static class FinitePolynomial
    {
        internal static BigInteger[] Trim(BigInteger[] a) { int n = a.Length; while (n > 0 && a[n - 1].IsZero) n--; return a.Take(n).ToArray(); }
        internal static BigInteger[] Rem(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            a = (BigInteger[])a.Clone(); var inv = BigInteger.ModPow(b[b.Length - 1], p - 2, p);
            for (int i = a.Length - 1; i >= b.Length - 1; i--)
            { var q = Mod(a[i] * inv, p); for (int k = 0; k < b.Length; k++) a[i - b.Length + 1 + k] = Mod(a[i - b.Length + 1 + k] - q * b[k], p); }
            return Trim(a);
        }
        internal static int RootCount(BigInteger[] f, BigInteger p, CancellationToken token)
        {
            f = Trim(f.Select(x => Mod(x, p)).ToArray());
            if (f.Length == 0) throw new ArgumentException("Zero polynomial.");
            BigInteger[] Mul(BigInteger[] a, BigInteger[] b)
            { if (a.Length == 0 || b.Length == 0) return Array.Empty<BigInteger>(); var c = new BigInteger[a.Length + b.Length - 1]; for (int i = 0; i < a.Length; i++) for (int j = 0; j < b.Length; j++) c[i + j] = Mod(c[i + j] + a[i] * b[j], p); return Rem(c, f, p); }
            var power = new BigInteger[] { 0, 1 }; var result = new BigInteger[] { 1 };
            for (var k = p; k > 0; k >>= 1) { token.ThrowIfCancellationRequested(); if (!k.IsEven) result = Mul(result, power); if (k > 1) power = Mul(power, power); }
            if (result.Length < 2) Array.Resize(ref result, 2); result[1] = Mod(result[1] - 1, p); result = Trim(result);
            while (result.Length != 0) { var r = Rem(f, result, p); f = result; result = r; }
            return f.Length - 1;
        }
    }
}

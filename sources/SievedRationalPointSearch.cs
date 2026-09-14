using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Search on an integral Weierstrass model. At x=m/d² the square condition is
    // F(m)=4m³+b2*m²*d²+2b4*m*d⁴+b6*d⁶ = (2d³y+a1*m*d+a3*d³)².
    // Each modular mask rejects only nonsquares; the remaining candidates are checked exactly.
    internal sealed class SievedRationalPointSearch
    {
        private const int WordsPerBlock = 1024;
        private readonly EllipticCurveQ curve;
        private readonly RankLowerBoundSearchOptions options;
        private readonly Stopwatch watch;
        private readonly CancellationToken token;
        internal long SquareTests { get; private set; }
        internal long Blocks { get; private set; }

        internal SievedRationalPointSearch(EllipticCurveQ curve, RankLowerBoundSearchOptions options,
            Stopwatch watch, CancellationToken token)
        { this.curve = curve; this.options = options; this.watch = watch; this.token = token; }

        private bool Expired()
        {
            token.ThrowIfCancellationRequested();
            return watch.Elapsed >= options.TimeLimit;
        }

        internal RankLowerBoundSearchStopReason Run(Func<EllipticCurvePoint, RankLowerBoundSearchStopReason?> accept)
        {
            var words = new ulong[WordsPerBlock];
            var first = options.NumeratorCenter - options.NumeratorRadius;
            var last = options.NumeratorCenter + options.NumeratorRadius;
            for (int d = 1; d <= options.DenominatorRootBound; d++)
            {
                if (Expired()) return RankLowerBoundSearchStopReason.TimeLimit;
                if (SquareTests >= options.MaxSquareTests) return RankLowerBoundSearchStopReason.SquareTestLimit;
                BigInteger d2 = (BigInteger)d * d, d3 = d2 * d;
                BigInteger c2 = curve.B2.Num * d2, c1 = 2 * curve.B4.Num * d2 * d2, c0 = curve.B6.Num * d3 * d3;
                var masks = BuildMasks(c2, c1, c0);
                for (var start = first; start <= last;)
                {
                    if (Expired()) return RankLowerBoundSearchStopReason.TimeLimit;
                    int length = (int)BigInteger.Min(WordsPerBlock * 64, last - start + 1);
                    int count = (length + 63) / 64;
                    for (int w = 0; w < count; w++) words[w] = ulong.MaxValue;
                    if (length % 64 != 0) words[count - 1] = (1UL << (length % 64)) - 1;
                    Blocks++;
                    foreach (var mask in masks)
                    {
                        if (Expired()) return RankLowerBoundSearchStopReason.TimeLimit;
                        int residue = Mod(start, mask.Prime), step = 64 % mask.Prime;
                        for (int w = 0; w < count; w++)
                        {
                            words[w] &= mask.Words[residue];
                            residue += step;
                            if (residue >= mask.Prime) residue -= mask.Prime;
                        }
                    }
                    for (int w = 0; w < count; w++)
                    {
                        ulong bits = words[w];
                        while (bits != 0)
                        {
                            if (Expired()) return RankLowerBoundSearchStopReason.TimeLimit;
                            int bit = TrailingZeros(bits);
                            bits &= bits - 1;
                            BigInteger m = start + 64 * w + bit;
                            if (BigInteger.GreatestCommonDivisor(m, d) != 1) continue;
                            if (SquareTests >= options.MaxSquareTests) return RankLowerBoundSearchStopReason.SquareTestLimit;
                            SquareTests++;
                            BigInteger value = ((4 * m + c2) * m + c1) * m + c0;
                            if (value.Sign < 0) continue;
                            var root = InternalMath.IntegerSqrt(value);
                            if (root * root != value) continue;
                            var point = new EllipticCurvePoint(new BigRational(m, d2),
                                new BigRational(root - curve.A1.Num * m * d - curve.A3.Num * d3, 2 * d3));
                            var reason = accept(point);
                            if (reason.HasValue) return reason.Value;
                        }
                    }
                    start += length;
                }
            }
            return RankLowerBoundSearchStopReason.SearchBoxExhausted;
        }

        private sealed class Mask
        {
            internal int Prime;
            internal ulong[] Words;
        }

        private List<Mask> BuildMasks(BigInteger c2, BigInteger c1, BigInteger c0)
        {
            var masks = new List<Mask>();
            for (int p = 3; p <= 251; p += 2)
            {
                token.ThrowIfCancellationRequested();
                bool prime = true;
                for (int q = 3; q * q <= p; q += 2) if (p % q == 0) { prime = false; break; }
                if (!prime) continue;
                var square = new bool[p];
                for (int x = 0; x < p; x++) square[x * x % p] = true;
                int a = Mod(c2, p), b = Mod(c1, p), c = Mod(c0, p);
                var allowed = new bool[p];
                int survivors = 0;
                for (int x = 0; x < p; x++)
                {
                    allowed[x] = square[(((4 * x + a) * x + b) * x + c) % p];
                    if (allowed[x]) survivors++;
                }
                if (survivors == p) continue;
                var mask = new Mask { Prime = p, Words = new ulong[p] };
                for (int x = 0; x < p; x++)
                {
                    ulong word = 0;
                    for (int bit = 0; bit < 64; bit++) if (allowed[(x + bit) % p]) word |= 1UL << bit;
                    mask.Words[x] = word;
                }
                masks.Add(mask);
            }
            return masks;
        }

        private static int Mod(BigInteger n, int p)
        {
            int r = (int)(n % p);
            return r < 0 ? r + p : r;
        }

        private static int TrailingZeros(ulong bits)
        {
            int result = 0;
            if ((bits & 0xffffffffUL) == 0) { result += 32; bits >>= 32; }
            if ((bits & 0xffffUL) == 0) { result += 16; bits >>= 16; }
            if ((bits & 0xffUL) == 0) { result += 8; bits >>= 8; }
            if ((bits & 0xfUL) == 0) { result += 4; bits >>= 4; }
            if ((bits & 3UL) == 0) { result += 2; bits >>= 2; }
            if ((bits & 1UL) == 0) result++;
            return result;
        }
    }
}

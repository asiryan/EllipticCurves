using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>
        /// Compute unconditional algebraic rank bounds locally. With rational 2-torsion,
        /// use descent by 2-isogeny, exact quartic point searches and local obstructions.
        /// Otherwise return a lower bound of 0 or 1 and an unknown (null) upper bound.
        /// Equality of the bounds certifies the exact rank; no BSD or parity assumption is used.
        /// </summary>
        /// <param name="searchBound">Nonnegative bound on each primitive quartic coordinate.
        /// In the fallback, bounds both numerator and denominator of the searched x coordinates.
        /// Zero disables point search. Increasing this can improve the lower bound.</param>
        /// <param name="maxSquareClasses">Maximum number of signed square classes per isogeny.
        /// Exceeding this limit throws NotSupportedException, rather than truncating the upper bound.</param>
        /// <param name="cancellationToken">Cancels factorization, enumeration and point searches.</param>
        public RankBounds GetRankBounds(int searchBound = 32, int maxSquareClasses = 65536,
            CancellationToken cancellationToken = default)
        {
            if (searchBound < 0 || searchBound == int.MaxValue) throw new ArgumentOutOfRangeException(nameof(searchBound));
            if (maxSquareClasses < 2) throw new ArgumentOutOfRangeException(nameof(maxSquareClasses));
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("A singular curve has no Mordell-Weil rank.");
            var e = GetGlobalMinimalModel(cancellationToken);

            // Complete the square integrally: X=4x, Y=8y+4a1*x+4a3.
            // Y^2 = X^3 + b2 X^2 + 8 b4 X + 16 b6.
            var a = e.B2.Num;
            var b = 8 * e.B4.Num;
            var c = 16 * e.B6.Num;
            if (!TryCubicRoot(a, b, c, cancellationToken, out var root))
                return new RankBounds(FindPositiveRank(e, searchBound, cancellationToken), null, false);

            b += 2 * a * root + 3 * root * root;
            a += 3 * root;
            // Remove pure square scaling to make the quartic searches more effective.
            foreach (var p in Factor(b, cancellationToken).Keys)
            {
                var p2 = p * p;
                var p4 = p2 * p2;
                while (a % p2 == 0 && b % p4 == 0) { a /= p2; b /= p4; }
            }
            var first = IsogenyImageBounds(a, b, searchBound, maxSquareClasses, cancellationToken);
            var second = IsogenyImageBounds(-2 * a, a * a - 4 * b, searchBound, maxSquareClasses, cancellationToken);
            int lower = Math.Max(0, first.lower + second.lower - 2);
            int upper = first.upper + second.upper - 2;
            // A point may be divisible in the isogeny quotients and still have infinite order.
            if (lower == 0 && upper > 0) lower = FindPositiveRank(e, searchBound, cancellationToken);
            if (upper < lower) throw new InvalidOperationException("Inconsistent certified rank bounds.");
            return new RankBounds(lower, upper, true);
        }

        private static bool TryCubicRoot(BigInteger a, BigInteger b, BigInteger c,
            CancellationToken token, out BigInteger root)
        {
            root = 0;
            if (c.IsZero) return true;
            // A rational root of a monic integral polynomial is an integral divisor of c.
            foreach (var d in Divisors(Factor(c, token)))
            {
                token.ThrowIfCancellationRequested();
                if (((d + a) * d + b) * d + c == 0) { root = d; return true; }
                if (((-d + a) * d - b) * d + c == 0) { root = -d; return true; }
            }
            return false;
        }

        private static int FindPositiveRank(EllipticCurveQ e, int bound, CancellationToken token)
        {
            // Mazur's theorem: a rational torsion point has order at most 12.
            // Exact additions suffice to certify one point of infinite order.
            for (int denominator = 1; denominator <= bound; denominator++)
            for (long numerator = -(long)bound; numerator <= bound; numerator++)
            {
                token.ThrowIfCancellationRequested();
                if (BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator), denominator) != 1) continue;
                var x = new BigRational(numerator, denominator);
                var t = e.A1 * x + e.A3;
                var rhs = x * x * x + e.A2 * x * x + e.A4 * x + e.A6 + t * t / 4;
                if (!BigRational.IsSquare(rhs, out var y)) continue;
                var point = new EllipticCurvePoint(x, y - t / 2);
                var multiple = point;
                int order = 1;
                while (order <= 12 && !multiple.IsInfinity) { multiple = e.Add(multiple, point); order++; }
                if (order > 12) return 1;
            }
            return 0;
        }

        private static (int lower, int upper) IsogenyImageBounds(BigInteger a, BigInteger b,
            int bound, int maxClasses, CancellationToken token)
        {
            var factors = Factor(b, token).OrderBy(pair => pair.Key).ToArray();
            long classCount = 2;
            foreach (var unused in factors)
            {
                classCount *= 2;
                if (classCount > maxClasses)
                    throw new NotSupportedException("The 2-isogeny descent exceeds maxSquareClasses. Increase the limit or use a smaller model.");
            }
            int count = (int)classCount;
            var values = new BigInteger[count];
            values[0] = 1;
            values[1] = -1;
            int used = 2, torsionClass = b.Sign < 0 ? 1 : 0;
            for (int i = 0; i < factors.Length; i++)
            {
                if ((factors[i].Value & 1) != 0) torsionClass |= 1 << (i + 1);
                for (int j = 0; j < used; j++) values[used + j] = values[j] * factors[i].Key;
                used *= 2;
            }
            var basis = new int[factors.Length + 1];
            AddSquareClass(basis, torsionClass);
            int survivors = 0;
            // These are necessary tests only. Passing does not assert Q_p-solubility.
            // Therefore omitted primes and finite precision can only weaken the upper bound.
            var sieves = new[] { (2, 256), (3, 81), (5, 25), (7, 49),
                (11, 11), (13, 13), (17, 17), (19, 19), (23, 23), (29, 29), (31, 31) };
            foreach (int mask in Enumerable.Range(0, count))
            {
                token.ThrowIfCancellationRequested();
                var d = values[mask];
                var other = Divide(b, d);
                // F(u,v)=d*u^4+a*u^2*v^2+(b/d)*v^4. If both end coefficients
                // are negative, its maximum as a quadratic in (u/v)^2 must be nonnegative.
                bool possible = !(d < 0 && other < 0 && (a <= 0 || a * a - 4 * b < 0));
                foreach (var sieve in sieves)
                {
                    if (!possible) break;
                    possible = QuarticHasResidue(d, a, other, sieve.Item1, sieve.Item2);
                }
                if (!possible) continue;
                survivors++;
                if (ReducesToZero(basis, mask)) continue;
                if (QuarticHasPoint(d, a, other, bound, token)) AddSquareClass(basis, mask);
            }
            if (survivors == 0) throw new InvalidOperationException("The trivial descent class was eliminated.");
            int upper = 0;
            for (int size = survivors; size > 1; size >>= 1) upper++;
            // The actual image is an F_2-space contained in the survivors, so its dimension
            // is <= floor(log2(survivors)), even if the finite sieve survivors are not a group.
            int lower = basis.Count(x => x != 0);
            if (lower > upper) throw new InvalidOperationException("Inconsistent isogeny image bounds.");
            return (lower, upper);
        }

        private static bool QuarticHasResidue(BigInteger d, BigInteger a, BigInteger other, int p, int modulus)
        {
            var squares = new bool[modulus];
            for (int i = 0; i < modulus; i++) squares[i * i % modulus] = true;
            long D = (long)Mod(d, modulus), A = (long)Mod(a, modulus), B = (long)Mod(other, modulus);
            // Primitive projective charts: v=1, or u=1 and p|v.
            for (int x = 0; x < modulus; x++)
            {
                long x2 = (long)x * x % modulus, x4 = x2 * x2 % modulus;
                if (squares[(int)((D * x4 + A * x2 + B) % modulus)]) return true;
                if (x % p == 0 && squares[(int)((D + A * x2 + B * x4) % modulus)]) return true;
            }
            return false;
        }

        private static bool QuarticHasPoint(BigInteger d, BigInteger a, BigInteger other, int bound, CancellationToken token)
        {
            if (bound == 0) return false;
            for (int u = 0; u <= bound; u++)
            for (int v = 0; v <= bound; v++)
            {
                token.ThrowIfCancellationRequested();
                if (BigInteger.GreatestCommonDivisor(u, v) != 1) continue;
                BigInteger u2 = (BigInteger)u * u, v2 = (BigInteger)v * v;
                var square = d * u2 * u2 + a * u2 * v2 + other * v2 * v2;
                if (square < 0) continue;
                var root = InternalMath.IntegerSqrt(square);
                if (root * root == square) return true;
            }
            return false;
        }

        private static bool ReducesToZero(int[] basis, int value)
        {
            for (int i = basis.Length - 1; i >= 0; i--)
                if ((value & (1 << i)) != 0) value ^= basis[i];
            return value == 0;
        }

        private static void AddSquareClass(int[] basis, int value)
        {
            for (int i = basis.Length - 1; i >= 0; i--)
            {
                if ((value & (1 << i)) == 0) continue;
                if (basis[i] == 0) { basis[i] = value; return; }
                value ^= basis[i];
            }
        }
    }
}

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
        /// <summary>Local invariants at all bad primes, in increasing order.</summary>
        public IReadOnlyList<LocalReductionData> LocalData => GetLocalData();
        /// <summary>Product of local Tamagawa numbers at the finite bad primes.</summary>
        public BigInteger TamagawaProduct => GetLocalData().Aggregate(BigInteger.One, (a, d) => a * d.TamagawaNumber);

        /// <summary>Compute exact invariants at all bad primes.</summary>
        public IReadOnlyList<LocalReductionData> GetLocalData(CancellationToken cancellationToken = default)
        {
            var e = GetGlobalMinimalModel(cancellationToken);
            return Array.AsReadOnly(Factor(e.Discriminant.Num, cancellationToken).OrderBy(x => x.Key)
                .Select(x => ComputeLocalData(e, x.Key, x.Value, cancellationToken)).ToArray());
        }

        /// <summary>Compute exact invariants at a prime, including good reduction.</summary>
        public LocalReductionData GetLocalData(BigInteger prime, CancellationToken cancellationToken = default)
        {
            if (!IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime), "A prime is required.");
            var e = GetGlobalMinimalModel(cancellationToken);
            return ComputeLocalData(e, prime, Valuation(e.Discriminant.Num, prime), cancellationToken);
        }

        // Tate's algorithm, using a globally minimal model. Arithmetic criteria:
        // Cremona, Algorithms for Modular Elliptic Curves, III.2; no foreign source code.
        internal static LocalReductionData ComputeLocalData(EllipticCurveQ model, BigInteger p, int n, CancellationToken token)
        {
            int j = Valuation(model.JInvariant.Den, p), root = LocalRootNumber(model, p);
            LocalReductionData Result(string symbol, int f, int cp, ReductionType type = ReductionType.Additive)
                => new LocalReductionData(p, n, f, j, symbol, type, cp, root);
            if (n == 0) return Result("I0", 0, 1, ReductionType.Good);
            BigInteger Inv(BigInteger a) => BigInteger.ModPow(Mod(a, p), p - 2, p);
            bool HasQuadraticRoot(BigInteger a, BigInteger b, BigInteger c)
            {
                a = Mod(a, p); b = Mod(b, p); c = Mod(c, p);
                if (p == 2) return c.IsZero || Mod(a + b + c, p).IsZero;
                if (a.IsZero) return !b.IsZero || c.IsZero;
                var d = Mod(b * b - 4 * a * c, p);
                return d.IsZero || BigInteger.ModPow(d, (p - 1) / 2, p) == 1;
            }
            var e = model;
            BigInteger r, t, s;
            if (p == 2)
            {
                r = Mod(e.B2.Num, p) == 0 ? Mod(e.A4.Num, p) : Mod(e.A3.Num, p);
                t = Mod(e.B2.Num, p) == 0 ? Mod(r * (1 + e.A2.Num + e.A4.Num) + e.A6.Num, p) : Mod(r + e.A4.Num, p);
            }
            else if (p == 3)
            {
                r = Mod(e.B2.Num, p) == 0 ? Mod(-e.B6.Num, p) : Mod(-e.B2.Num * e.B4.Num, p);
                t = Mod(e.A1.Num * r + e.A3.Num, p);
            }
            else
            {
                r = Mod(e.C4.Num, p) == 0 ? Mod(-e.B2.Num * Inv(12), p) : Mod(-(e.C6.Num + e.B2.Num * e.C4.Num) * Inv(12 * e.C4.Num), p);
                t = Mod(-(e.A1.Num * r + e.A3.Num) * Inv(2), p);
            }
            e = Translate(e, r, 0, t);
            if (Mod(e.C4.Num, p) != 0)
            {
                bool split = HasQuadraticRoot(1, e.A1.Num, -e.A2.Num);
                return Result("I" + n, 1, split ? n : (n % 2 == 0 ? 2 : 1), split ? ReductionType.SplitMultiplicative : ReductionType.NonSplitMultiplicative);
            }
            var p2 = p * p; var p3 = p2 * p;
            if (e.A6.Num % p2 != 0) return Result("II", n, 1);
            if (e.B8.Num % p3 != 0) return Result("III", n - 1, 2);
            if (e.B6.Num % p3 != 0) return Result("IV", n - 2, HasQuadraticRoot(1, Divide(e.A3.Num, p), -Divide(e.A6.Num, p2)) ? 3 : 1);
            s = p == 2 ? Mod(e.A2.Num, 2) : -e.A1.Num * Inv(2);
            t = p == 2 ? 2 * Mod(Divide(e.A6.Num, 4), 2) : -e.A3.Num * Inv(2);
            e = Translate(e, 0, s, t);
            var b = Divide(e.A2.Num, p); var c = Divide(e.A4.Num, p2); var d = Divide(e.A6.Num, p3);
            var disc = b * b * c * c - 4 * c * c * c - 4 * b * b * b * d - 27 * d * d + 18 * b * c * d;
            if (Mod(disc, p) != 0) return Result("I0*", n - 4, 1 + FinitePolynomial.RootCount(new[] { d, c, b, BigInteger.One }, p, token));
            bool triple;
            if (p <= 3)
            {
                r = 0;
                while (r < p && (Mod(((r + b) * r + c) * r + d, p) != 0 || Mod(3 * r * r + 2 * b * r + c, p) != 0)) r++;
                if (r == p) throw new InvalidOperationException("Tate auxiliary cubic has no repeated root.");
                triple = Mod(b + 3 * r, p) == 0 && Mod(c + 2 * b * r + 3 * r * r, p) == 0;
            }
            else
            {
                var h = Mod(3 * c - b * b, p); triple = h.IsZero;
                r = triple ? Mod(-b * Inv(3), p) : Mod((b * c - 9 * d) * Inv(2 * h), p);
            }
            e = Translate(e, p * r, 0, 0);
            if (!triple)
            {
                int m = 1; BigInteger mx = p2, my = p2;
                while (m <= n)
                {
                    token.ThrowIfCancellationRequested();
                    var aa = Divide(e.A3.Num, my); var cc = Divide(e.A6.Num, mx * my);
                    if (Mod(aa * aa + 4 * cc, p) != 0) return Result("I" + m + "*", n - m - 4, HasQuadraticRoot(1, aa, -cc) ? 4 : 2);
                    t = my * Mod(p == 2 ? cc : -aa * Inv(2), p);
                    e = Translate(e, 0, 0, t); my *= p; m++;
                    var a2 = Divide(e.A2.Num, p); var a4 = Divide(e.A4.Num, p * mx); var a6 = Divide(e.A6.Num, mx * my);
                    if (Mod(a4 * a4 - 4 * a2 * a6, p) != 0) return Result("I" + m + "*", n - m - 4, HasQuadraticRoot(a2, a4, a6) ? 4 : 2);
                    r = mx * Mod(p == 2 ? a6 * a2 : -a4 * Inv(2 * a2), p);
                    e = Translate(e, r, 0, 0); mx *= p; m++;
                }
                throw new InvalidOperationException("Tate refinement exceeded the discriminant valuation.");
            }
            var a3 = Divide(e.A3.Num, p2); var a6star = Divide(e.A6.Num, p2 * p2);
            if (Mod(a3 * a3 + 4 * a6star, p) != 0) return Result("IV*", n - 6, HasQuadraticRoot(1, a3, -a6star) ? 3 : 1);
            t = p2 * Mod(p == 2 ? a6star : -a3 * Inv(2), p);
            e = Translate(e, 0, 0, t);
            if (e.A4.Num % (p2 * p2) != 0) return Result("III*", n - 7, 2);
            if (e.A6.Num % (p3 * p3) != 0) return Result("II*", n - 8, 1);
            throw new InvalidOperationException("Unexpected nonminimal model in Tate's algorithm.");
        }
    }
}

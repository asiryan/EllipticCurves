using System;
using System.Linq;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>The reduced global minimal integral model, computed locally over Q.</summary>
        public EllipticCurveQ GlobalMinimalModel => GetGlobalMinimalModel();

        /// <summary>The exact conductor, computed locally, including wild reduction at 2 and 3.</summary>
        public BigInteger Conductor => GetConductor();

        /// <summary>
        /// Compute the reduced global minimal integral model using invariant scaling and
        /// integral reconstruction. Integer factorization can be expensive; cancellation is supported.
        /// </summary>
        public EllipticCurveQ GetGlobalMinimalModel(CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("A singular curve has no elliptic minimal model.");
            var (c4, c6, delta) = InternalMath.IntegralInvariants(this);
            foreach (var p in Factor(delta, cancellationToken).Keys.OrderBy(p => p))
            {
                var p4 = BigInteger.Pow(p, 4);
                var p6 = BigInteger.Pow(p, 6);
                var p12 = p6 * p6;
                while (delta % p12 == 0 && c4 % p4 == 0 && c6 % p6 == 0)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    if (!TryIntegralModel(c4 / p4, c6 / p6, out _)) break;
                    c4 /= p4;
                    c6 /= p6;
                    delta /= p12;
                }
            }
            if (!TryIntegralModel(c4, c6, out var model))
                throw new InvalidOperationException("Unable to reconstruct integral invariants.");
            return model;
        }

        /// <summary>
        /// Compute N = product p^f_p on a minimal model, using Tate's algorithm at 2 and 3.
        /// This method performs no network access and does not require external mathematical software.
        /// </summary>
        public BigInteger GetConductor(CancellationToken cancellationToken = default)
        {
            var model = GetGlobalMinimalModel(cancellationToken);
            BigInteger conductor = 1;
            foreach (var pair in Factor(model.Discriminant.Num, cancellationToken))
            {
                cancellationToken.ThrowIfCancellationRequested();
                var p = pair.Key;
                int exponent = model.C4.Num % p != 0 ? 1
                    : p > 3 ? 2 : WildConductorExponent(model, (int)p, pair.Value, cancellationToken);
                if (exponent < 1) throw new InvalidOperationException("Invalid conductor exponent on a minimal model.");
                conductor *= BigInteger.Pow(p, exponent);
            }
            return conductor;
        }

        // A reduced integral model has a1,a3 in {0,1} and a2 in {-1,0,1}.
        // Trying these 12 possibilities is an exact form of Kraus's integrality test.
        private static bool TryIntegralModel(BigInteger c4, BigInteger c6, out EllipticCurveQ model)
        {
            for (int a1 = 0; a1 <= 1; a1++)
            for (int a2 = -1; a2 <= 1; a2++)
            {
                BigInteger b2 = a1 * a1 + 4 * a2;
                if ((b2 * b2 - c4) % 24 != 0) continue;
                var b4 = (b2 * b2 - c4) / 24;
                var numerator = -b2 * b2 * b2 + 36 * b2 * b4 - c6;
                if (numerator % 216 != 0) continue;
                var b6 = numerator / 216;
                for (int a3 = 0; a3 <= 1; a3++)
                {
                    if ((b4 - a1 * a3) % 2 != 0 || (b6 - a3 * a3) % 4 != 0) continue;
                    model = new EllipticCurveQ(a1, a2, a3,
                        new BigRational((b4 - a1 * a3) / 2), new BigRational((b6 - a3 * a3) / 4));
                    return true;
                }
            }
            model = null;
            return false;
        }

        // Translation x=X+r, y=Y+sX+t (no scaling).
        private static EllipticCurveQ Translate(EllipticCurveQ e, BigInteger r, BigInteger s, BigInteger t)
        {
            BigRational R = new BigRational(r), S = new BigRational(s), T = new BigRational(t);
            return new EllipticCurveQ(e.A1 + 2 * S,
                e.A2 - S * e.A1 + 3 * R - S * S,
                e.A3 + R * e.A1 + 2 * T,
                e.A4 - S * e.A3 + 2 * R * e.A2 - (T + R * S) * e.A1 + 3 * R * R - 2 * S * T,
                e.A6 + R * e.A4 + R * R * e.A2 + R * R * R - T * e.A3 - T * T - R * T * e.A1);
        }

        // Tate's algorithm: component counts give f=v(Delta)+1-components.
        // Only p=2,3 reach here, and the input is already globally minimal.
        // Mathematical reference: Cremona, Algorithms for Modular Elliptic Curves, III.3.2.
        private static int WildConductorExponent(EllipticCurveQ e, int p, int n, CancellationToken token)
        {
            BigInteger r, s, t;
            if (p == 2)
            {
                r = Mod(e.A4.Num, 2);
                t = Mod(r * (1 + e.A2.Num + e.A4.Num) + e.A6.Num, 2);
            }
            else
            {
                r = Mod(e.B2.Num % 3 == 0 ? -e.B6.Num : -e.B2.Num * e.B4.Num, 3);
                t = Mod(e.A1.Num * r + e.A3.Num, 3);
            }
            e = Translate(e, r, 0, t);
            var p2 = new BigInteger(p * p);
            var p3 = p2 * p;
            if (e.A6.Num % p2 != 0) return n;               // II
            if (e.B8.Num % p3 != 0) return n - 1;           // III
            if (e.B6.Num % p3 != 0) return n - 2;           // IV

            s = p == 2 ? Mod(e.A2.Num, 2) : -2 * e.A1.Num;
            // Keep the factor a3 here: reducing t modulo 3 would lose the
            // extra precision needed to make a3 divisible by 9 and a6 by 27.
            t = p == 2 ? 2 * Mod(Divide(e.A6.Num, 4), 2) : -2 * e.A3.Num;
            e = Translate(e, 0, s, t);
            var b = Divide(e.A2.Num, p);
            var c = Divide(e.A4.Num, p2);
            var d = Divide(e.A6.Num, p3);
            var discriminant = b * b * c * c - 4 * c * c * c - 4 * b * b * b * d - 27 * d * d + 18 * b * c * d;
            if (discriminant % p != 0) return n - 4;       // I0*

            // Find the repeated root directly in F_2 or F_3; derivative also detects triple roots.
            r = 0;
            while (r < p && (Mod(r * r * r + b * r * r + c * r + d, p) != 0
                || Mod(3 * r * r + 2 * b * r + c, p) != 0)) r++;
            if (r == p) throw new InvalidOperationException("Missing repeated root in Tate's algorithm.");
            bool triple = Mod(b + 3 * r, p).IsZero && Mod(c + 2 * b * r + 3 * r * r, p).IsZero;
            e = Translate(e, p * r, 0, 0);
            if (!triple)
            {
                int m = 1;
                BigInteger mx = p2, my = p2;
                while (m <= n)
                {
                    token.ThrowIfCancellationRequested();
                    var a3 = Divide(e.A3.Num, my);
                    var a6 = Divide(e.A6.Num, mx * my);
                    if ((a3 * a3 + 4 * a6) % p != 0) return n - m - 4;
                    t = my * Mod(p == 2 ? a6 : -2 * a3, p);
                    e = Translate(e, 0, 0, t);
                    my *= p;
                    m++;
                    var a2 = Divide(e.A2.Num, p);
                    var a4 = Divide(e.A4.Num, p * mx);
                    a6 = Divide(e.A6.Num, mx * my);
                    if ((a4 * a4 - 4 * a2 * a6) % p != 0) return n - m - 4;
                    r = mx * Mod(p == 2 ? a6 * a2 : -2 * a2 * a4, p);
                    e = Translate(e, r, 0, 0);
                    mx *= p;
                    m++;
                }
                throw new InvalidOperationException("Tate refinement exceeded the discriminant valuation.");
            }

            var x3 = Divide(e.A3.Num, p2);
            var x6 = Divide(e.A6.Num, p2 * p2);
            if ((x3 * x3 + 4 * x6) % p != 0) return n - 6; // IV*
            t = p2 * Mod(p == 2 ? x6 : -2 * x3, p);
            e = Translate(e, 0, 0, t);
            if (e.A4.Num % (p2 * p2) != 0) return n - 7;   // III*
            if (e.A6.Num % (p3 * p3) != 0) return n - 8;   // II*
            throw new InvalidOperationException("Tate's algorithm encountered a non-minimal model.");
        }
    }
}

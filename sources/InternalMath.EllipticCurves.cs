using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    internal static partial class InternalMath
    {
        internal readonly struct InvariantsInt
        {
            public readonly BigInteger B2;
            public readonly BigInteger B4;
            public readonly BigInteger B6;
            public readonly BigInteger B8;
            public readonly BigInteger C4;
            public readonly BigInteger C6;
            public readonly BigInteger Delta;

            public InvariantsInt(BigInteger b2, BigInteger b4, BigInteger b6, BigInteger b8, BigInteger c4, BigInteger c6, BigInteger delta)
            {
                B2 = b2;
                B4 = b4;
                B6 = b6;
                B8 = b8;
                C4 = c4;
                C6 = c6;
                Delta = delta;
            }
        }

        internal static InvariantsInt ComputeInvariants(BigInteger a1, BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6)
        {
            // b2 = a1^2 + 4a2
            BigInteger b2 = a1 * a1 + 4 * a2;
            // b4 = a1*a3 + 2a4
            BigInteger b4 = a1 * a3 + 2 * a4;
            // b6 = a3^2 + 4a6
            BigInteger b6 = a3 * a3 + 4 * a6;
            // b8 = a1^2*a6 - a1*a3*a4 + 4a2*a6 + a2*a3^2 - a4^2
            BigInteger b8 = (a1 * a1) * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * (a3 * a3) - a4 * a4;

            // c4 = b2^2 - 24 b4
            BigInteger c4 = b2 * b2 - 24 * b4;
            // c6 = -b2^3 + 36 b2 b4 - 216 b6
            BigInteger c6 = -(b2 * b2 * b2) + 36 * b2 * b4 - 216 * b6;
            // Delta = -b2^2 b8 - 8 b4^3 - 27 b6^2 + 9 b2 b4 b6
            BigInteger delta = -(b2 * b2) * b8 - 8 * (b4 * b4 * b4) - 27 * (b6 * b6) + 9 * b2 * b4 * b6;

            return new InvariantsInt(b2, b4, b6, b8, c4, c6, delta);
        }

        internal static int Valuation(BigInteger n, BigInteger p)
        {
            if (p.Sign <= 0) throw new ArgumentOutOfRangeException(nameof(p), "Prime must be positive.");
            if (p.IsOne) throw new ArgumentOutOfRangeException(nameof(p), "Prime must be >= 2.");
            if (n.IsZero) return int.MaxValue;

            n = BigInteger.Abs(n);
            int v = 0;
            while (!n.IsZero)
            {
                BigInteger rem;
                n = BigInteger.DivRem(n, p, out rem);
                if (!rem.IsZero) break;
                v++;
            }
            return v;
        }

        internal static bool Divides(BigInteger p, BigInteger n) => (n % p).IsZero;

        internal static bool DividesPower(BigInteger p, BigInteger n, int exponent)
        {
            if (exponent <= 0) return true;
            if (p.Sign <= 0 || p.IsOne) throw new ArgumentOutOfRangeException(nameof(p));
            if (n.IsZero) return true;
            n = BigInteger.Abs(n);
            for (int i = 0; i < exponent; i++)
            {
                BigInteger rem;
                n = BigInteger.DivRem(n, p, out rem);
                if (!rem.IsZero) return false;
            }
            return true;
        }

        internal static BigInteger Mod(BigInteger a, BigInteger m)
        {
            if (m.Sign <= 0) throw new ArgumentOutOfRangeException(nameof(m), "Modulus must be positive.");
            BigInteger r = a % m;
            if (r.Sign < 0) r += m;
            return r;
        }

        internal static BigInteger ModInverse(BigInteger a, BigInteger p)
        {
            // Extended Euclid. Assumes gcd(a,p)=1.
            a = Mod(a, p);
            if (a.IsZero) throw new ArithmeticException("0 has no inverse modulo p.");

            BigInteger t = BigInteger.Zero;
            BigInteger newT = BigInteger.One;
            BigInteger r = p;
            BigInteger newR = a;

            while (!newR.IsZero)
            {
                BigInteger q = r / newR;
                (t, newT) = (newT, t - q * newT);
                (r, newR) = (newR, r - q * newR);
            }

            if (r != BigInteger.One) throw new ArithmeticException("Value is not invertible modulo p.");
            return Mod(t, p);
        }

        internal static int LegendreSymbol(BigInteger a, BigInteger p)
        {
            if (p == 2) throw new ArgumentOutOfRangeException(nameof(p), "Legendre symbol is for odd primes.");
            a = Mod(a, p);
            if (a.IsZero) return 0;
            BigInteger ls = BigInteger.ModPow(a, (p - 1) / 2, p);
            if (ls.IsOne) return 1;
            if (ls == p - 1) return -1;
            // Should not happen for prime p.
            return 0;
        }

        internal static bool QuadRoots(BigInteger a, BigInteger b, BigInteger c, BigInteger p)
        {
            a = Mod(a, p);
            b = Mod(b, p);
            c = Mod(c, p);

            if (p == 2)
            {
                // Brute force x in {0,1}.
                for (int x = 0; x < 2; x++)
                {
                    BigInteger xx = x;
                    BigInteger v = (a * xx * xx + b * xx + c) % p;
                    if (v.IsZero) return true;
                }
                return false;
            }

            if (a.IsZero)
            {
                // Linear congruence b x + c == 0.
                if (b.IsZero) return c.IsZero;
                return true;
            }

            BigInteger disc = Mod(b * b - 4 * a * c, p);
            int ls = LegendreSymbol(disc, p);
            return ls >= 0;
        }

        internal static BigInteger ExactDiv(BigInteger numerator, BigInteger denominator, string? context = null)
        {
            BigInteger rem;
            BigInteger q = BigInteger.DivRem(numerator, denominator, out rem);
            if (!rem.IsZero)
                throw new ArithmeticException(context is null
                    ? $"Non-exact division: {numerator} / {denominator}"
                    : $"Non-exact division in {context}: {numerator} / {denominator}");
            return q;
        }

        internal static (BigInteger a1, BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6)
            ReducedMinimalModelFromInvariants(BigInteger c4, BigInteger c6)
        {
            // Implements the Kraus–Laska–Connell algorithm (Cremona, Algorithms for EC, §3.2)
            // to compute the unique reduced global minimal model from integral invariants c4,c6
            // satisfying Kraus’s conditions.

            BigInteger delta = ExactDiv(BigInteger.Pow(c4, 3) - c6 * c6, 1728, context: "Δ=(c4^3-c6^2)/1728");
            if (delta.IsZero) throw new ArgumentException("Discriminant must be non-zero.", nameof(c4));

            BigInteger u = BigInteger.One;
            BigInteger g = BigInteger.GreatestCommonDivisor(c6 * c6, delta);
            g = BigInteger.Abs(g);

            // Prime divisors of g.
            foreach (var kv in FactorAbs(g))
            {
                BigInteger p = kv.Key;
                int vp = Valuation(g, p);
                int d = vp / 12;
                if (d <= 0) continue;

                if (p == 2)
                {
                    BigInteger pow4d = BigInteger.One << (4 * d); // 2^(4d)
                    BigInteger pow6d = BigInteger.One << (6 * d); // 2^(6d)

                    BigInteger a = Mod(ExactDiv(c4, pow4d, "c4/2^(4d)"), 16);
                    BigInteger b = Mod(ExactDiv(c6, pow6d, "c6/2^(6d)"), 32);
                    // If (b mod 4 != -1) and not(a=0 and (b=0 or b=8)) then d=d-1
                    if (Mod(b, 4) != 3 && !(a.IsZero && (b.IsZero || b == 8))) d -= 1;
                }
                else if (p == 3)
                {
                    int v3c6 = Valuation(c6, 3);
                    if (v3c6 == 6 * d + 2) d -= 1;
                }

                if (d <= 0) continue;
                u *= BigInteger.Pow(p, d);
            }

            // Scale invariants down.
            BigInteger u2 = u * u;
            BigInteger u4 = u2 * u2;
            BigInteger u6 = u4 * u2;
            c4 = ExactDiv(c4, u4, "c4/u^4");
            c6 = ExactDiv(c6, u6, "c6/u^6");

            // Recover reduced integral model from invariants.
            // b2 = -c6 mod 12 in {-5,...,6}
            BigInteger b2 = Mod(-c6, 12);
            if (b2 > 6) b2 -= 12;

            BigInteger b4 = ExactDiv(b2 * b2 - c4, 24, "b4=(b2^2-c4)/24");
            BigInteger b6 = ExactDiv(-b2 * b2 * b2 + 36 * b2 * b4 - c6, 216, "b6=(-b2^3+36*b2*b4-c6)/216");

            BigInteger a1 = Mod(b2, 2);
            BigInteger a3 = Mod(b6, 2);

            BigInteger a2 = ExactDiv(b2 - a1, 4, "a2=(b2-a1)/4");
            BigInteger a4 = ExactDiv(b4 - a1 * a3, 2, "a4=(b4-a1*a3)/2");
            BigInteger a6 = ExactDiv(b6 - a3, 4, "a6=(b6-a3)/4");

            return (a1, a2, a3, a4, a6);
        }

        private static void TransCoord(ref BigInteger a1, ref BigInteger a2, ref BigInteger a3, ref BigInteger a4, ref BigInteger a6,
            BigInteger r, BigInteger s, BigInteger t, BigInteger u)
        {
            // Transformation T(r,s,t,u):
            // x = u^2 x' + r
            // y = u^3 y' + s u^2 x' + t
            // See (3.1.3) in Cremona, Algorithms for EC.
            if (u.IsZero) throw new ArgumentOutOfRangeException(nameof(u), "u must be non-zero.");

            BigInteger ua1 = a1 + 2 * s;
            BigInteger u2a2 = a2 - s * a1 + 3 * r - s * s;
            BigInteger u3a3 = a3 + r * a1 + 2 * t;
            BigInteger u4a4 = a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t;
            BigInteger u6a6 = a6 + r * a4 + r * r * a2 + r * r * r - t * a3 - t * t - r * t * a1;

            if (u == BigInteger.One)
            {
                a1 = ua1;
                a2 = u2a2;
                a3 = u3a3;
                a4 = u4a4;
                a6 = u6a6;
                return;
            }

            BigInteger u2 = u * u;
            BigInteger u3 = u2 * u;
            BigInteger u4 = u2 * u2;
            BigInteger u6 = u3 * u3;

            a1 = ExactDiv(ua1, u, "a1' = (a1+2s)/u");
            a2 = ExactDiv(u2a2, u2, "a2' = (...) / u^2");
            a3 = ExactDiv(u3a3, u3, "a3' = (...) / u^3");
            a4 = ExactDiv(u4a4, u4, "a4' = (...) / u^4");
            a6 = ExactDiv(u6a6, u6, "a6' = (...) / u^6");
        }

        internal static (string Kodaira, int ConductorExponent, int LocalIndex) TateAlgorithm(
            BigInteger a1, BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6,
            BigInteger p)
        {
            if (p.Sign <= 0 || p.IsOne) throw new ArgumentOutOfRangeException(nameof(p), "p must be a prime >= 2.");

            // Port of the pseudocode in Cremona (Algorithms for EC, §3.2: Tate’s algorithm).
            while (true)
            {
                var inv = ComputeInvariants(a1, a2, a3, a4, a6);
                int n = Valuation(inv.Delta, p);

                // Type I0.
                if (n == 0) return ("I0", 0, 1);

                // Change coordinates so that p | a3, a4, a6.
                BigInteger r, s, t;
                if (p == 2)
                {
                    if (Divides(p, inv.B2))
                    {
                        r = Mod(a4, p);
                        t = Mod(r * (1 + a2 + a4) + a6, p);
                    }
                    else
                    {
                        r = Mod(a3, p);
                        t = Mod(r + a4, p);
                    }
                }
                else if (p == 3)
                {
                    r = Divides(p, inv.B2) ? Mod(-inv.B6, p) : Mod(-inv.B2 * inv.B4, p);
                    t = Mod(a1 * r + a3, p);
                }
                else
                {
                    if (Divides(p, inv.C4))
                        r = Mod(-ModInverse(12, p) * inv.B2, p);
                    else
                        r = Mod(-ModInverse(12 * inv.C4, p) * (inv.C6 + inv.B2 * inv.C4), p);

                    t = Mod(-ModInverse(2, p) * (a1 * r + a3), p);
                }

                TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, r, 0, t, 1);

                // Recompute invariants (Δ unchanged under u=1 transforms, but b/c4/b6 etc can change).
                inv = ComputeInvariants(a1, a2, a3, a4, a6);

                // Test for types In (multiplicative reduction).
                if (!Divides(p, inv.C4))
                {
                    int cp;
                    if (QuadRoots(1, a1, -a2, p)) cp = n;
                    else if ((n & 1) == 0) cp = 2;
                    else cp = 1;
                    return ($"I{n}", 1, cp);
                }

                // Types II, III, IV.
                if (!DividesPower(p, a6, 2)) return ("II", n, 1);
                if (!DividesPower(p, inv.B8, 3)) return ("III", n - 1, 2);
                if (!DividesPower(p, inv.B6, 3))
                {
                    BigInteger p2 = p * p;
                    int cp = QuadRoots(1, a3 / p, -a6 / p2, p) ? 3 : 1;
                    return ("IV", n - 2, cp);
                }

                // Change coordinates so that p|a1,a2; p^2|a3,a4; p^3|a6.
                if (p == 2)
                {
                    s = Mod(a2, 2);
                    t = 2 * Mod(a6 / 4, 2);
                }
                else
                {
                    s = Mod(-a1 * ModInverse(2, p), p);
                    t = Mod(-a3 * ModInverse(2, p), p);
                }
                TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, 0, s, t, 1);
                inv = ComputeInvariants(a1, a2, a3, a4, a6);

                // Auxiliary cubic T^3 + b T^2 + c T + d.
                BigInteger p2b = p * p;
                BigInteger p3 = p2b * p;
                BigInteger b = a2 / p;
                BigInteger c = a4 / p2b;
                BigInteger d = a6 / p3;

                BigInteger w = 27 * d * d - inv.B2 * c * c + 4 * b * b * b * d - 18 * b * c * d + 4 * c * c * c;
                BigInteger x = 3 * c - inv.B2;

                // Distinct roots: type I0*.
                if (!Divides(p, w))
                {
                    int cp = 1 + CountCubicRootsModP(b, c, d, p);
                    return ("I0*", n - 4, cp);
                }

                // Double root: type I_m*.
                if (!Divides(p, x))
                {
                    // Change coordinates so that the double root is T ≡ 0.
                    if (p == 2) r = c;
                    else if (p == 3) r = b * c;
                    else r = (b * c - 9 * d) * ModInverse(2 * x, p);
                    r = p * Mod(r, p);
                    TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, r, 0, 0, 1);
                    inv = ComputeInvariants(a1, a2, a3, a4, a6);

                    int m = 1;
                    BigInteger mx = p * p;
                    BigInteger my = p * p;
                    int cp = 0;

                    while (cp == 0)
                    {
                        BigInteger xa2 = a2 / p;
                        BigInteger xa3 = a3 / my;
                        BigInteger xa4 = a4 / (p * mx);
                        BigInteger xa6 = a6 / (mx * my);

                        if (!Divides(p, xa3 * xa3 + 4 * xa6))
                        {
                            cp = QuadRoots(1, xa3, -xa6, p) ? 4 : 2;
                        }
                        else
                        {
                            if (p == 2) t = my * xa6;
                            else t = my * Mod(-xa3 * ModInverse(2, p), p);
                            TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, 0, 0, t, 1);

                            my *= p;
                            m++;

                            xa2 = a2 / p;
                            xa3 = a3 / my;
                            xa4 = a4 / (p * mx);
                            xa6 = a6 / (mx * my);

                            if (!Divides(p, xa4 * xa4 - 4 * xa2 * xa6))
                            {
                                cp = QuadRoots(xa2, xa4, xa6, p) ? 4 : 2;
                            }
                            else
                            {
                                if (p == 2) r = mx * Mod(xa6 * xa2, 2);
                                else r = mx * Mod(-xa4 * ModInverse(2 * xa2, p), p);
                                TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, r, 0, 0, 1);
                                mx *= p;
                                m++;
                            }
                        }
                    }

                    return ($"I{m}*", n - m - 4, cp);
                }

                // Triple root case: types IV*, III*, II* or non-minimal.
                BigInteger rp = (p == 3) ? -d : -b * ModInverse(3, p);
                r = p * Mod(rp, p);
                TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, r, 0, 0, 1);

                BigInteger p4 = p2b * p2b;
                BigInteger x3 = a3 / p2b;
                BigInteger x6 = a6 / p4;

                // Type IV*.
                if (!Divides(p, x3 * x3 + 4 * x6))
                {
                    int cp = QuadRoots(1, x3, -x6, p) ? 3 : 1;
                    return ("IV*", n - 6, cp);
                }

                // Change coordinates so that p^3 | a3, p^5 | a6.
                t = (p == 2) ? x6 : x3 * ModInverse(2, p);
                t = -p2b * Mod(t, p);
                TransCoord(ref a1, ref a2, ref a3, ref a4, ref a6, 0, 0, t, 1);

                // Types III*, II*.
                if (!DividesPower(p, a4, 4)) return ("III*", n - 7, 2);
                if (!DividesPower(p, a6, 6)) return ("II*", n - 8, 1);

                // Equation non-minimal: divide each a_i by p^i and start again.
                // (Scaling transform with u=p.)
                a1 = ExactDiv(a1, p, "a1/p");
                a2 = ExactDiv(a2, p2b, "a2/p^2");
                a3 = ExactDiv(a3, p3, "a3/p^3");
                a4 = ExactDiv(a4, p4, "a4/p^4");
                BigInteger p6 = p3 * p3;
                a6 = ExactDiv(a6, p6, "a6/p^6");
            }
        }

        private static int CountCubicRootsModP(BigInteger b, BigInteger c, BigInteger d, BigInteger p)
        {
            // Count roots of T^3 + b T^2 + c T + d mod p.
            // For small p (fits in int and not too large) brute force; otherwise use
            // a deterministic method based on gcd(f, T^p - T).
            if (p <= 2000)
            {
                int pi = (int)p;
                int roots = 0;
                for (int x = 0; x < pi; x++)
                {
                    BigInteger xx = x;
                    BigInteger v = (xx * xx * xx + b * xx * xx + c * xx + d) % p;
                    if (v.IsZero) roots++;
                }
                return roots;
            }

            // Polynomial arithmetic over F_p for degree <=3.
            // f(T) = T^3 + b T^2 + c T + d
            // Compute g = gcd(f, T^p - T) in F_p[T]; deg(g) is number of distinct roots.
            var f = new BigInteger[4];
            f[0] = Mod(d, p);
            f[1] = Mod(c, p);
            f[2] = Mod(b, p);
            f[3] = BigInteger.One;

            // Compute T^p mod f.
            var xpoly = new BigInteger[2];
            xpoly[0] = BigInteger.Zero;
            xpoly[1] = BigInteger.One;

            var tp = PolyPowMod(xpoly, p, f, p);
            // h(T) = T^p - T
            var h = PolySub(tp, xpoly, p);
            var g = PolyGcd(f, h, p);
            return PolyDegree(g);
        }

        private static int PolyDegree(BigInteger[] poly)
        {
            for (int i = poly.Length - 1; i >= 0; i--)
            {
                if (!poly[i].IsZero) return i;
            }
            return -1;
        }

        private static BigInteger[] PolyTrim(BigInteger[] poly)
        {
            int deg = PolyDegree(poly);
            if (deg < 0) return new BigInteger[1];
            if (deg == poly.Length - 1) return poly;
            var r = new BigInteger[deg + 1];
            Array.Copy(poly, r, deg + 1);
            return r;
        }

        private static BigInteger[] PolyAdd(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            int n = Math.Max(a.Length, b.Length);
            var r = new BigInteger[n];
            for (int i = 0; i < n; i++)
            {
                BigInteger ai = i < a.Length ? a[i] : BigInteger.Zero;
                BigInteger bi = i < b.Length ? b[i] : BigInteger.Zero;
                r[i] = Mod(ai + bi, p);
            }
            return PolyTrim(r);
        }

        private static BigInteger[] PolySub(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            int n = Math.Max(a.Length, b.Length);
            var r = new BigInteger[n];
            for (int i = 0; i < n; i++)
            {
                BigInteger ai = i < a.Length ? a[i] : BigInteger.Zero;
                BigInteger bi = i < b.Length ? b[i] : BigInteger.Zero;
                r[i] = Mod(ai - bi, p);
            }
            return PolyTrim(r);
        }

        private static BigInteger[] PolyMul(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            var r = new BigInteger[a.Length + b.Length - 1];
            for (int i = 0; i < a.Length; i++)
            {
                if (a[i].IsZero) continue;
                for (int j = 0; j < b.Length; j++)
                {
                    if (b[j].IsZero) continue;
                    r[i + j] = Mod(r[i + j] + a[i] * b[j], p);
                }
            }
            return PolyTrim(r);
        }

        private static BigInteger[] PolyMod(BigInteger[] a, BigInteger[] modPoly, BigInteger p)
        {
            // Reduce a modulo modPoly (assumed monic) in F_p[T].
            a = PolyTrim(a);
            int degM = PolyDegree(modPoly);
            if (degM <= 0) throw new ArgumentException("modPoly must have positive degree", nameof(modPoly));
            BigInteger invLead = BigInteger.One; // monic

            var r = (BigInteger[])a.Clone();
            while (true)
            {
                int degR = PolyDegree(r);
                if (degR < degM) break;
                BigInteger lead = r[degR];
                if (!lead.IsZero)
                {
                    BigInteger coeff = Mod(lead * invLead, p);
                    int shift = degR - degM;
                    // r -= coeff * T^shift * modPoly
                    for (int i = 0; i <= degM; i++)
                    {
                        int idx = i + shift;
                        r[idx] = Mod(r[idx] - coeff * modPoly[i], p);
                    }
                }
                r = PolyTrim(r);
            }
            return r;
        }

        private static BigInteger[] PolyPowMod(BigInteger[] a, BigInteger exp, BigInteger[] modPoly, BigInteger p)
        {
            // Exponentiation by squaring: a^exp mod modPoly.
            var result = new BigInteger[1];
            result[0] = BigInteger.One;
            var basePoly = PolyMod(a, modPoly, p);

            BigInteger e = exp;
            while (e > 0)
            {
                if (!e.IsEven)
                {
                    result = PolyMod(PolyMul(result, basePoly, p), modPoly, p);
                }
                e >>= 1;
                if (e.IsZero) break;
                basePoly = PolyMod(PolyMul(basePoly, basePoly, p), modPoly, p);
            }
            return result;
        }

        private static BigInteger[] PolyGcd(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            a = PolyTrim(a);
            b = PolyTrim(b);
            while (PolyDegree(b) >= 0)
            {
                var r = PolyRemainder(a, b, p);
                a = b;
                b = r;
            }

            // Make monic.
            int degA = PolyDegree(a);
            if (degA < 0) return new BigInteger[1];
            BigInteger lead = a[degA];
            if (!lead.IsOne)
            {
                BigInteger inv = ModInverse(lead, p);
                for (int i = 0; i <= degA; i++) a[i] = Mod(a[i] * inv, p);
            }
            return PolyTrim(a);
        }

        private static BigInteger[] PolyRemainder(BigInteger[] a, BigInteger[] b, BigInteger p)
        {
            a = PolyTrim(a);
            b = PolyTrim(b);
            int degB = PolyDegree(b);
            if (degB < 0) throw new ArgumentException("Division by zero polynomial", nameof(b));

            var r = (BigInteger[])a.Clone();
            int degR;
            while ((degR = PolyDegree(r)) >= degB)
            {
                BigInteger leadR = r[degR];
                if (!leadR.IsZero)
                {
                    BigInteger invLeadB = ModInverse(b[degB], p);
                    BigInteger coeff = Mod(leadR * invLeadB, p);
                    int shift = degR - degB;
                    for (int i = 0; i <= degB; i++)
                    {
                        r[i + shift] = Mod(r[i + shift] - coeff * b[i], p);
                    }
                }
                r = PolyTrim(r);
            }
            return r;
        }

        /// <summary>
        /// Returns true iff n is a perfect square in ℤ.
        /// </summary>
        public static bool IsSquare(BigInteger n)
        {
            if (n.Sign < 0) return false;
            if (n.IsZero || n.IsOne) return true;

            // Fast modular filters.
            // Squares mod 16 are 0,1,4,9.
            int mod16 = (int)(n & 0xF);
            if (mod16 != 0 && mod16 != 1 && mod16 != 4 && mod16 != 9) return false;

            // Squares mod 3 are 0,1.
            int mod3 = (int)(n % 3);
            if (mod3 < 0) mod3 += 3;
            if (mod3 == 2) return false;

            var r = IntegerSqrt(n);
            return r * r == n;
        }

        /// <summary>
        /// Integer square root: ⌊sqrt(n)⌋ for n ≥ 0.
        /// </summary>
        

        private static int BitLength(BigInteger n)
        {
            if (n.Sign == 0) return 0;
            n = BigInteger.Abs(n);

            // Two's-complement, little-endian.
            // For non-negative numbers, the most-significant byte may be 0x00 to avoid a sign bit.
            byte[] bytes = n.ToByteArray();
            int len = bytes.Length;
            while (len > 1 && bytes[len - 1] == 0) len--;

            byte msb = bytes[len - 1];
            int bits = (len - 1) * 8;

            int b = msb;
            int msbBits = 0;
            while (b != 0) { msbBits++; b >>= 1; }
            return bits + msbBits;
        }
    }
}

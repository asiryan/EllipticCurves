using System;
using System.Numerics;
using System.Threading;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>The exact sign of the functional equation of L(E,s).</summary>
        public int RootNumber => GetRootNumber();

        /// <summary>Compute the global root number, including the local signs at 2 and 3, without network access.</summary>
        public int GetRootNumber(CancellationToken cancellationToken = default)
        {
            var e = GetGlobalMinimalModel(cancellationToken);
            int sign = -1; // The real place.
            foreach (var p in Factor(e.Discriminant.Num, cancellationToken).Keys)
            {
                cancellationToken.ThrowIfCancellationRequested();
                sign *= LocalRootNumber(e, p);
            }
            return sign;
        }

        internal static int LocalRootNumber(EllipticCurveQ e, BigInteger p)
        {
            if (e.Discriminant.Num % p != 0) return 1;
            if (p == 2) return RootAtTwo(e);
            if (p == 3) return RootAtThree(e);
            var c4 = e.C4.Num;
            var c6 = e.C6.Num;
            if (c4 % p != 0) return -QuadraticCharacter(-c6, p);
            int d = Valuation(e.Discriminant.Num, p), a = Valuation(c4, p);
            if (3L * a < d) return QuadraticCharacter(-1, p); // Potentially multiplicative.
            return QuadraticCharacter(d == 3 || d == 9 ? -2 : d == 4 || d == 8 ? -3 : -1, p);
        }

        private static int QuadraticCharacter(BigInteger value, BigInteger p)
        {
            var r = BigInteger.ModPow(Mod(value, p), (p - 1) / 2, p);
            return r.IsZero ? 0 : r.IsOne ? 1 : -1;
        }

        // Infinity must remain larger than all finite valuations after normalization.
        private static int Valuation(BigInteger value, BigInteger p)
        {
            if (value.IsZero) return int.MaxValue / 4;
            int v = 0;
            while (value % p == 0) { value /= p; v++; }
            return v;
        }

        private sealed class RootInvariants
        {
            internal readonly int A, B, D;
            internal readonly BigInteger U4, U6, UD, C4, C6;
            private readonly int p;
            internal RootInvariants(EllipticCurveQ e, int prime)
            {
                p = prime;
                int a = Valuation(e.C4.Num, p), b = Valuation(e.C6.Num, p), d = Valuation(e.Discriminant.Num, p);
                int m = Math.Min(a / 4, Math.Min(b / 6, d / 12));
                A = a - 4 * m; B = b - 6 * m; D = d - 12 * m;
                C4 = e.C4.Num / BigInteger.Pow(p, 4 * m);
                C6 = e.C6.Num / BigInteger.Pow(p, 6 * m);
                U4 = e.C4.Num.IsZero ? BigInteger.Zero : e.C4.Num / BigInteger.Pow(p, a);
                U6 = e.C6.Num.IsZero ? BigInteger.Zero : e.C6.Num / BigInteger.Pow(p, b);
                UD = e.Discriminant.Num / BigInteger.Pow(p, d);
            }
            internal BigInteger C4At(int exponent) => Divide(C4, BigInteger.Pow(p, exponent));
            internal BigInteger C6At(int exponent) => Divide(C6, BigInteger.Pow(p, exponent));
        }

        private static bool Residue(BigInteger n, int modulus, params int[] values)
        {
            int r = (int)Mod(n, modulus);
            foreach (int value in values) if (r == value) return true;
            return false;
        }
        private static int SignIf(bool positive) => positive ? 1 : -1;

        // Cowland Kellock--Dokchitser, Root numbers and parity phenomena (2023), Appendix A.
        // Columns are (v(Delta),v(c6),v(c4)); one original Rizzo row is retained below.
        private static int RootAtTwo(EllipticCurveQ e)
        {
            var v = new RootInvariants(e, 2);
            int a = v.A, b = v.B, d = v.D;
            var u = v.U4; var w = v.U6;
            if (a == 0 && b == 0)
                return SignIf(Residue(w, 4, 3) && (d == 0 || Residue(w, 8, 3)));
            if (d == 0 && b == 3 && a == 3)
                return SignIf(Residue(u, 4, 1) ? Residue(w, 8, 1, 7) : Residue(w, 8, 1, 3));
            if (d == 0 && b == 3 && a >= 4) return SignIf(Residue(w, 4, 1));
            if (d == 0 && b >= 4 && a == 2)
            {
                if (Residue(u, 4, 3)) return SignIf(b == 4);
                if (b == 4) return SignIf(Residue(u + 4 * w, 16, 9, 13));
                // Rizzo, Table III: the same condition applies for every b>=5.
                // Keep this row from the original table (the 2023 table reverses the b=5 case).
                return SignIf(Residue(u + 4 * v.C6At(4), 16, 5, 9));
            }
            if (d == 1 && b == 3 && a == 2) return SignIf(Residue(u + 4 * w, 16, 3) || Residue(u, 16, 11));
            if (d == 2 && b == 3 && a == 2) return SignIf(Residue(v.UD - w, 4, 0));
            if (d == 2 && b == 4 && a == 3) return SignIf(Residue(u + w, 8, 0, 6));
            if (d == 2 && b == 4 && a >= 4) return SignIf(Residue(w, 4, 1));
            if (d == 3 && b == 3 && a == 2) return SignIf(Residue(v.UD, 4, 3));
            if (d == 3 && b == 5 && a == 3) return SignIf(Residue(2 * w + u, 8, 1, 3));
            if (d == 3 && b >= 6 && a == 3) return SignIf(Residue(u, 8, 5, 7));
            if (d >= 4 && b == 3 && a == 2) return SignIf(Residue(w, 4, 3));
            if (d == 4 && b == 5 && a == 4)
                return SignIf(Residue(u - w, 4, 0) ? Residue(u, 4, 1) : Residue(u, 4, 1) && Residue(u * w, 8, 3));
            if (d == 4 && b == 5 && a >= 5)
                return SignIf(Residue(w, 4, 3) ? a == 5 : a == 5 && Residue(w, 8, 5));
            if (d == 6 && b == 6 && a == 5) return SignIf(Residue(u, 4, 3));
            if (d == 6 && b == 6 && a >= 6) return SignIf(Residue(w, 4, 1));
            if (d == 6 && b >= 7 && a == 4)
                return SignIf(Residue(u, 4, 1) ? b == 7 : Residue(u - 4 * v.C6At(7), 16, 7, 11));
            if (d == 7 && b == 6 && a == 4) return SignIf(Residue(w, 8, 5) || Residue(w - 5 * u, 8, 0));
            if (d == 8 && b == 6 && a == 4) return SignIf(Residue(2 * w + u, 16, 3) || Residue(2 * w + u, 32, 23));
            if (d == 8 && b == 7 && a == 5) return SignIf(Residue(2 * u + w, 8, 7) || Residue(w, 8, 3));
            if (d == 8 && b == 7 && a >= 6)
                return SignIf(a == 6 && (Residue(w, 4, 3) || Residue(2 * u + w, 8, 3)));
            if (d == 9 && b == 6 && a == 4) return SignIf(Residue(2 * w + u, 32, 11) || Residue(w, 8, 7));
            if (d == 9 && b == 8 && a == 5) return SignIf(Residue(2 * w + u, 8, 1, 7));
            if (d == 9 && b >= 9 && a == 5) return SignIf(Residue(u, 8, 1, 3));
            if (d == 10 && b == 6 && a == 4)
                return SignIf(Residue(w, 4, 1) || Residue(u - 2 * w, 64, 3, 19));
            if (d == 10 && b == 8 && a == 6) return SignIf(Residue(u * w, 4, 3));
            if (d == 10 && b == 8 && a >= 7) return SignIf(Residue(w, 4, 1));
            if (d == 11 && b == 6 && a == 4) return SignIf(Residue(w, 8, 1, 3, 5));
            throw new InvalidOperationException($"Unrecognized root-number invariants at 2: ({a},{b},{d}).");
        }

        // Rizzo, Average Root Numbers for a Nonconstant Family of Elliptic Curves (2003), Table II.
        private static int RootAtThree(EllipticCurveQ e)
        {
            var v = new RootInvariants(e, 3);
            int a = v.A, b = v.B, d = v.D;
            var u = v.U4; var w = v.U6;
            if (a == 0 && b == 0) return d == 0 ? 1 : SignIf(Residue(w, 3, 1));
            if (a == 1 && b >= 2 && d == 0) return 1;
            if (a >= 2 && b == 2 && d == 1) return SignIf(Residue(w, 3, 1));
            if (a >= 2 && b == 3 && d == 3)
                return SignIf(Residue(w * w + 2 - 3 * v.C4At(2), 9, 0) || Residue(w, 9, 4, 7, 8));
            if (a == 2 && b == 4 && d == 3) return SignIf(!Residue(u - w, 3, 0));
            if (a == 2 && b >= 5 && d == 3) return 1;
            if (a == 2 && b == 3 && d == 4) return 1;
            if (a == 2 && b == 3 && d == 5) return SignIf(Residue(v.UD - w, 3, 0));
            if (a >= 3 && b == 4 && d == 5) return SignIf(Residue(w, 3, 2));
            if (a == 2 && b == 3 && d >= 6) return -1;
            if (a == 3 && b == 5 && d == 6) return SignIf(Residue(u, 3, 2));
            if (a == 3 && b >= 6 && d == 6) return -1;
            if (a >= 4 && b == 5 && d == 7) return SignIf(Residue(w, 3, 2));
            if (a >= 4 && b == 6 && d == 9)
                return SignIf(Residue(w * w + 2 - 3 * v.C4At(4), 9, 0)
                    || (a == 4 ? Residue(w, 9, 4, 8) : Residue(w, 9, 1, 2)));
            if (a == 4 && b == 7 && d == 9) return SignIf(Residue(w, 3, 2));
            if (a == 4 && b >= 8 && d == 9) return 1;
            if (a == 4 && b == 6 && d == 10) return SignIf(Residue(w, 9, 2, 7));
            if (a == 4 && b == 6 && d == 11) return SignIf(Residue(w, 3, 1));
            if (a >= 5 && b == 7 && d == 11) return SignIf(Residue(w, 3, 1));
            throw new InvalidOperationException($"Unrecognized root-number invariants at 3: ({a},{b},{d}).");
        }
    }
}

using System;
using System.Numerics;
using static EllipticCurves.NativeNumberTheory;

namespace EllipticCurves
{
    internal sealed class RootNumberInvariants
    {
        internal readonly int A, B, D;
        internal readonly BigInteger U4, U6, UD, C4, C6;
        private readonly int p;
        internal RootNumberInvariants(EllipticCurveQ e, int prime)
        {
            p = prime;
            int a = EllipticCurveQ.Valuation(e.C4.Num, p), b = EllipticCurveQ.Valuation(e.C6.Num, p), d = EllipticCurveQ.Valuation(e.Discriminant.Num, p);
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
}

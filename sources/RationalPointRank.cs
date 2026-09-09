using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    // Exact lower-bound certificates from good reduction and Kummer characters.
    // For every simple root r of f mod p, P -> (x(P)-r | p) is a homomorphism
    // E(F_p) -> F_2, with f'(r) at (r,0) and 1 at infinity.
    // The combined image dimension is <= rank(E(Q)) + dim E(Q)[2].
    internal sealed class RationalPointRank
    {
        private readonly EllipticCurveQ curve;
        private readonly DescentBudget budget;
        private readonly int torsionDimension;
        private readonly List<(int prime, int root, int derivative)> characters = new List<(int, int, int)>();
        private readonly List<BigInteger> basis = new List<BigInteger>();
        private readonly HashSet<EllipticCurvePoint> points = new HashSet<EllipticCurvePoint>();
        private int imageDimension;
        internal int LowerBound { get; private set; }
        internal int ImageDimension => imageDimension;

        internal RationalPointRank(EllipticCurveQ curve, int torsionDimension, DescentBudget budget)
        {
            this.curve = curve; this.torsionDimension = torsionDimension; this.budget = budget;
            var composite = new bool[budget.Options.ReductionPrimeBound + 1];
            var delta = curve.Discriminant.Num;
            // X=4x, Y=8y+4a1*x+4a3 gives an integral monic cubic.
            var a = curve.B2.Num; var b = 8 * curve.B4.Num; var c = 16 * curve.B6.Num;
            for (int p = 2; p < composite.Length; p++)
            {
                budget.Token.ThrowIfCancellationRequested();
                if (composite[p]) continue;
                for (int multiple = p * 2; multiple < composite.Length; multiple += p) composite[multiple] = true;
                if (p == 2 || delta % p == 0) continue;
                long aa = (long)NativeNumberTheory.Mod(a, p), bb = (long)NativeNumberTheory.Mod(b, p), cc = (long)NativeNumberTheory.Mod(c, p);
                for (int r = 0; r < p; r++)
                    if ((((r + aa) * r + bb) * r + cc) % p == 0)
                        characters.Add((p, r, (int)((3L * r * r + 2 * aa * r + bb) % p)));
            }
            for (int i = 0; i < characters.Count; i++) basis.Add(0);
        }

        internal void Add(EllipticCurvePoint point)
        {
            budget.Token.ThrowIfCancellationRequested();
            if (point.IsInfinity || !points.Add(point)) return;
            if (!curve.IsOnCurve(point)) throw new InvalidOperationException("Rank witness is not on the curve.");
            BigInteger image = 0;
            for (int i = 0; i < characters.Count; i++)
            {
                var ch = characters[i]; int p = ch.prime;
                if (point.X.Den % p == 0) continue; // good reduction sends P to infinity
                var x = NativeNumberTheory.Mod(4 * point.X.Num * BigInteger.ModPow(point.X.Den % p, p - 2, p), p);
                var value = NativeNumberTheory.Mod(x - ch.root, p);
                if (value.IsZero) value = ch.derivative;
                if (BigInteger.ModPow(value, (p - 1) / 2, p) == p - 1) image |= BigInteger.One << i;
            }
            for (int i = basis.Count - 1; i >= 0; i--)
            {
                if ((image & (BigInteger.One << i)).IsZero) continue;
                if (basis[i].IsZero) { basis[i] = image; imageDimension++; break; }
                image ^= basis[i];
            }
            LowerBound = Math.Max(LowerBound, imageDimension - torsionDimension);
            // An even multiple can vanish in every tested Kummer quotient. Preserve
            // the independent infinite-order certificate in that situation.
            if (LowerBound == 0)
            {
                var multiple = point;
                int n = 1;
                while (n <= 12 && !multiple.IsInfinity)
                { budget.Token.ThrowIfCancellationRequested(); multiple = curve.Add(multiple, point); n++; }
                if (n > 12) LowerBound = 1;
            }
        }

        internal void Search()
        {
            int bound = budget.Options.SearchBound;
            // On an integral Weierstrass equation every rational x denominator
            // in lowest terms is a square. Keep the existing numerator/denominator box.
            for (long denominatorRoot = 1; denominatorRoot * denominatorRoot <= bound; denominatorRoot++)
            for (long numerator = -(long)bound; numerator <= bound; numerator++)
            {
                if (!budget.PointStep()) return;
                if (BigInteger.GreatestCommonDivisor(numerator, denominatorRoot) != 1) continue;
                var x = new BigRational(numerator, denominatorRoot * denominatorRoot);
                var linear = curve.A1 * x + curve.A3;
                var square = x * x * x + curve.A2 * x * x + curve.A4 * x + curve.A6 + linear * linear / 4;
                if (BigRational.IsSquare(square, out var y)) Add(new EllipticCurvePoint(x, y - linear / 2));
            }
        }
    }
}

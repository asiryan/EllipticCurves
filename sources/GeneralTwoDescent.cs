using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EllipticCurves
{
    // General 2-descent by reduced binary quartics. This implementation follows
    // the mathematical construction in Cremona, Algorithms for Modular Elliptic
    // Curves, III.3.6 (invariants, reduction regions and rational equivalence).
    internal sealed class GeneralTwoDescent
    {
        private sealed class Covering
        {
            internal readonly BinaryQuartic Form;
            internal readonly int[] Signature;
            internal bool HasPoint;
            internal Covering(BinaryQuartic form, int[] signature) { Form = form; Signature = signature; }
        }
        private readonly EllipticCurveQ curve;
        private readonly DescentBudget budget;
        private readonly RationalPointRank points;
        private readonly int torsionDimension;
        private readonly List<Covering> coverings = new List<Covering>();
        private readonly BigInteger[] badPrimes;
        private BigRational modelScale;
        private int solubleClasses = 1; // the identity class
        internal int LowerBound => Math.Max(points.LowerBound, CeilingLog2(solubleClasses) - torsionDimension);

        internal GeneralTwoDescent(EllipticCurveQ curve, int torsionDimension, RationalPointRank points, DescentBudget budget)
        {
            this.curve = curve; this.torsionDimension = torsionDimension; this.points = points; this.budget = budget;
            badPrimes = NativeNumberTheory.Factor(2 * curve.Discriminant.Num, budget.Token).Keys.OrderBy(p => p).ToArray();
        }

        internal int ComputeSelmerDimension()
        {
            var i = curve.C4.Num; var j = 2 * curve.C6.Num;
            modelScale = BigRational.One;
            // A reduced minimal elliptic model requires at most these two pairs.
            if (i % 16 == 0 && j % 64 == 0) { i /= 16; j /= 64; modelScale /= 2; }
            bool largePair = !(i % 4 == 0 && j % 8 == 0 && (2 * i + j) % 16 == 0);
            if (largePair) modelScale *= 2;
            Enumerate(i, j, largePair ? 4 : 1);
            if (largePair) Enumerate(16 * i, 64 * j, 1);
            int count = checked(coverings.Count + 1);
            if ((count & (count - 1)) != 0)
                throw new InvalidOperationException("The complete set of locally soluble quartic classes is not a power of two.");
            return CeilingLog2(count);
        }

        private static int CeilingLog2(int n)
        {
            int bits = 0;
            for (int m = n - 1; m > 0; m >>= 1) bits++;
            return bits;
        }

        private void Enumerate(BigInteger i, BigInteger j, int coefficientScale)
        {
            var roots = DescentPolynomial.RealRoots(new BigRational[] { new BigRational(j), new BigRational(-3 * i), 0, 1 },
                new BigRational(1, BigInteger.One << 96), budget);
            RationalInterval ii = i;
            // Regions can overlap. Keeping all boundary representatives and then
            // testing rational equivalence avoids tie-breaking and rounding assumptions.
            if (roots.Count == 3)
            {
                var small = roots[0]; var middle = roots[1]; var large = roots[2];
                var k = (4 * ii - large.Square()) / 3;
                var sqrtK = k.Sqrt();
                var typeOneBound = (k + sqrtK * large) / (3 * sqrtK + large + 2 * middle);
                Region(1, DescentPolynomial.Floor(typeOneBound.Upper),
                    (a, b) => (middle / 2 + RationalInterval.Exact(new BigRational(3 * b * b, 8 * a)),
                               large / 2 + RationalInterval.Exact(new BigRational(3 * b * b, 8 * a))));
                var numerator = ii - middle.Square();
                var positiveBound = numerator / (3 * (middle - small));
                Region(1, DescentPolynomial.Floor(positiveBound.Upper), (a, b) =>
                    ((4 * a * middle - 4 * numerator / 3 + (RationalInterval)(3 * b * b)) / (8 * a),
                     (4 * a * small + (RationalInterval)(3 * b * b)) / (8 * a)));
                var negativeBound = numerator / (3 * (large - middle));
                Region(-DescentPolynomial.Floor(negativeBound.Upper), -1, (a, b) =>
                    ((4 * a * large + (RationalInterval)(3 * b * b)) / (8 * a),
                     (4 * a * middle - 4 * numerator / 3 + (RationalInterval)(3 * b * b)) / (8 * a)));
            }
            else if (roots.Count == 1)
            {
                var phi = roots[0];
                var radius = (4 * (phi.Square() - ii) / 27).Sqrt();
                Region(DescentPolynomial.Ceiling((phi / 3 - radius).Lower), DescentPolynomial.Floor((phi / 3 + radius).Upper), (a, b) =>
                {
                    var low = ((RationalInterval)(9 * a * a + 3 * b * b) - 2 * a * phi + (4 * ii - phi.Square()) / 3) / (8 * BigInteger.Abs(a));
                    var high = (4 * a * phi + (RationalInterval)(3 * b * b)) / (8 * BigInteger.Abs(a));
                    return a.Sign > 0 ? (low, high) : (-high, -low);
                });
            }
            else throw new InvalidOperationException("A nonsingular resolvent cubic must have one or three distinct real roots.");

            void Region(BigInteger minA, BigInteger maxA, Func<BigInteger, BigInteger, (RationalInterval low, RationalInterval high)> cBounds)
            {
                for (var a = minA; a <= maxA; a++)
                {
                    budget.Step();
                    if (a.IsZero) continue; // a zero leading coefficient is the identity class
                    for (var b = -2 * BigInteger.Abs(a) + 1; b <= 2 * BigInteger.Abs(a); b++)
                    {
                        budget.Step();
                        var bounds = cBounds(a, b);
                        var minC = DescentPolynomial.Ceiling(bounds.low.Lower);
                        var maxC = DescentPolynomial.Floor(bounds.high.Upper);
                        for (var c = minC; c <= maxC; c++)
                        {
                            budget.Step();
                            var p = 3 * b * b - 8 * a * c;
                            var syzygy = p * p * p - 48 * i * a * a * p - 64 * j * a * a * a;
                            if (syzygy < 0 || syzygy % 27 != 0) continue;
                            var square = syzygy / 27;
                            var r = InternalMath.IntegerSqrt(square);
                            if (r * r != square) continue;
                            // Enumerate both signs; this makes reflection boundary cases explicit.
                            Candidate(r);
                            if (!r.IsZero) Candidate(-r);
                            void Candidate(BigInteger rr)
                            {
                                var dn = rr - b * b * b + 4 * a * b * c;
                                if (dn % (8 * a * a) != 0) return;
                                var d = dn / (8 * a * a);
                                var en = i + 3 * b * d - c * c;
                                if (en % (12 * a) != 0) return;
                                var form = new BinaryQuartic(a, b, c, d, en / (12 * a));
                                if (form.I != i || form.J != j) throw new InvalidOperationException("Quartic reconstruction failed its invariant check.");
                                Consider(form.Scale(coefficientScale));
                            }
                        }
                    }
                }
            }
        }

        private void Consider(BinaryQuartic q)
        {
            budget.Step();
            var signature = Signature(q);
            if (!signature.Contains(0) && q.HasRationalRoot(budget)) return;
            foreach (var previous in coverings)
            {
                if (!signature.SequenceEqual(previous.Signature) || !q.Equivalent(previous.Form, budget)) continue;
                if (!previous.HasPoint) TryPoint(q, previous);
                return;
            }
            if (!QuarticLocalSolubility.Everywhere(q, badPrimes, budget)) return;
            if (coverings.Count >= budget.Options.MaxSquareClasses - 1)
                throw new DescentLimitException("MaxSquareClasses was reached while enumerating the 2-Selmer group.");
            var covering = new Covering(q, signature);
            coverings.Add(covering);
            TryPoint(q, covering);
        }

        private void TryPoint(BinaryQuartic q, Covering covering)
        {
            if (!q.TryPoint(budget, out var u, out var v, out var y)) return;
            if (y.IsZero) throw new InvalidOperationException("A nontrivial covering unexpectedly has a rational branch point.");
            var image = q.MapPoint(u, v, y);
            var x = image.x / (36 * modelScale * modelScale) - curve.B2 / 12;
            var yy = image.y / (216 * BigRational.Pow(modelScale, 3)) - (curve.A1 * x + curve.A3) / 2;
            points.Add(new EllipticCurvePoint(x, yy));
            covering.HasPoint = true;
            solubleClasses++;
        }

        private static int[] Signature(BinaryQuartic q)
        {
            var result = new int[5]; int index = 0;
            var discriminant = 4 * q.I * q.I * q.I - q.J * q.J;
            foreach (int p in new[] { 5, 7, 11, 13, 17 })
            {
                if (discriminant % p == 0) { result[index++] = -1; continue; }
                int count = q.A % p == 0 ? 1 : 0;
                for (int x = 0; x < p; x++) if (q.Evaluate(x) % p == 0) count++;
                result[index++] = count;
            }
            return result;
        }
    }
}

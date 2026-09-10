using System;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Compute a period basis, the BSD real period and lattice area, using certified roots and AGM.</summary>
        public PeriodResult GetPeriods(RealComputationOptions options = null, CancellationToken cancellationToken = default)
        {
            var c = RealContext(options, cancellationToken); var e = GetGlobalMinimalModel(cancellationToken);
            // After completing the square: Y² = X³+b2 X²/4+b4 X/2+b6/4.
            var f = new BigRational[] { e.B6 / 4, e.B4 / 2, e.B2 / 4, 1 };
            var budget = new DescentBudget(new RankComputationOptions { MaxDescentWork = c.Options.MaxRootWork }, cancellationToken);
            System.Collections.Generic.List<RationalInterval> roots;
            try { roots = DescentPolynomial.RealRoots(f, new BigRational(1, c.Scale), budget); }
            catch (DescentLimitException ex) { throw new ArithmeticException("Period root-isolation work limit reached.", ex); }
            RealEnclosure Root(int i) => new RealEnclosure(roots[i].Lower, roots[i].Upper);
            RealEnclosure w, im;
            if (roots.Count == 3)
            {
                var a = c.Sqrt(c.Sub(Root(2), Root(0)));
                w = c.Div(c.Pi(), c.Agm(a, c.Sqrt(c.Sub(Root(2), Root(1)))));
                im = c.Div(c.Pi(), c.Agm(a, c.Sqrt(c.Sub(Root(1), Root(0)))));
            }
            else if (roots.Count == 1)
            {
                var r = Root(0);
                var separation = c.Div(c.Add(c.Mul(c.I(3), r), c.I(e.B2 / 4)), c.I(2));
                var distance = c.Sqrt(c.Add(c.Add(c.Mul(c.I(3), c.Mul(r, r)), c.Mul(c.I(e.B2 / 2), r)), c.I(e.B4 / 2)));
                var a = c.Sqrt(distance);
                w = c.Div(c.Pi(), c.Agm(a, c.Sqrt(c.Div(c.Add(distance, separation), c.I(2)))));
                im = c.Div(c.Pi(), c.Mul(c.I(2), c.Agm(a, c.Sqrt(c.Div(c.Sub(distance, separation), c.I(2))))));
            }
            else throw new ArithmeticException("Unexpected real-root count for a nonsingular cubic.");
            var re = roots.Count == 3 ? c.I(0) : c.Div(w, c.I(2));
            var real = c.Mul(w, c.I(e.NumberOfRealComponents)); var area = c.Mul(w, im);
            return new PeriodResult(e, c.Finish(w), c.Finish(re), c.Finish(im), c.Finish(real), c.Finish(area));
        }
    }
}

using System;
using System.Collections.Generic;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Numerical elliptic logarithm of a rational (hence real) point, using the minimal-model period basis.
        /// Handles both real components. Options control period/root preparation and iteration limits;
        /// the returned double logarithm is not certified and does not inherit DecimalDigits as an error bound.</summary>
        public RealEllipticLogarithmResult RealEllipticLogarithm(EllipticCurvePoint point,
            RealComputationOptions options = null, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("An elliptic logarithm requires a nonsingular curve.");
            if (!IsOnCurve(point)) throw new ArgumentException("Point is not on this curve.", nameof(point));
            var c = RealContext(options, cancellationToken); var map = GetMinimalModelIsomorphism(cancellationToken);
            var e = map.Target; var p = map.Map(point); var periods = e.GetPeriods(c.Options, cancellationToken);
            double period = periods.PrimitiveRealPeriod.Approximation;
            if (!(period > 0) || double.IsInfinity(period)) throw new ArithmeticException("Period is outside the double range.");
            if (periods.PrimitiveRealPeriod.LowerBound <= 0 || periods.PrimitiveRealPeriod.Width > periods.PrimitiveRealPeriod.LowerBound / 10000000000L)
                throw new ArithmeticException("Relative period precision is insufficient for the logarithm; increase PrecisionBits and DecimalDigits.");
            if (p.IsInfinity) return new RealEllipticLogarithmResult(e, 0, 0, 0, period);
            var polynomial = new BigRational[] { e.B6 / 4, e.B4 / 2, e.B2 / 4, 1 };
            var budget = new DescentBudget(new RankComputationOptions { MaxDescentWork = c.Options.MaxRootWork }, cancellationToken);
            List<RationalInterval> isolated;
            try { isolated = DescentPolynomial.RealRoots(polynomial, new BigRational(1, c.Scale), budget); }
            catch (DescentLimitException ex) { throw new ArithmeticException("Logarithm root-isolation work limit reached.", ex); }
            var y = p.Y + (e.A1 * p.X + e.A3) / 2;
            var roots = new RealEnclosure[isolated.Count];
            for (int i = 0; i < roots.Length; i++)
                roots[i] = y.IsZero && isolated[i].Lower <= p.X && p.X <= isolated[i].Upper ? c.I(p.X) :
                    new RealEnclosure(isolated[i].Lower, isolated[i].Upper);
            double Number(RealEnclosure value)
            {
                var middle = (value.LowerBound + value.UpperBound) / 2;
                if (value.LowerBound < 0 || (middle.IsZero ? !value.Width.IsZero : value.Width > middle / 10000000000L))
                    throw new ArithmeticException("Root precision is insufficient for the logarithm; increase PrecisionBits.");
                double result = RealEnclosure.ToDouble(middle);
                if (double.IsInfinity(result) || double.IsNaN(result) || (result == 0 && !middle.IsZero))
                    throw new ArithmeticException("Logarithm intermediate is outside the supported double range.");
                return result;
            }
            double Rf(double a, double b, double d) => CarlsonIntegral.Rf(a, b, d, c.Options.MaxIterations, cancellationToken);
            var x = c.I(p.X); double real; int component = 0;
            if (roots.Length == 3 && p.X < roots[2].LowerBound)
            {
                component = 1;
                // x=e3+(e2-e3) sin^2(theta). The imaginary part is omega2/2.
                var span = c.Sub(roots[1], roots[0]); var outer = c.Sub(roots[2], roots[0]);
                double sine = Math.Sqrt(Number(c.Div(c.Sub(x, roots[0]), span)));
                double cosineSquared = Number(c.Div(c.Sub(roots[1], x), span));
                double other = Number(c.Div(c.Sub(roots[2], x), outer));
                real = sine * Rf(cosineSquared, other, 1) / Math.Sqrt(Number(outer));
                if (y.Sign < 0) real = period - real;
            }
            else
            {
                // With r the largest real root and x-r=C^(1/2)s^2, the tail
                // integral becomes half an incomplete Legendre integral after s=tan(theta).
                var r = roots[roots.Length - 1]; var distance = c.Sub(x, r);
                var b = c.Add(c.Mul(c.I(3), r), c.I(e.B2 / 4));
                var derivative = c.Add(c.Add(c.Mul(c.I(3), c.Mul(r, r)), c.Mul(c.I(e.B2 / 2), r)), c.I(e.B4 / 2));
                var squareRoot = c.Sqrt(derivative); var fourthRoot = c.Sqrt(squareRoot);
                var complement = c.Div(c.Add(c.Mul(c.I(2), squareRoot), b), c.Mul(c.I(4), squareRoot));
                double q = Number(complement), h = Number(fourthRoot);
                var ratio = c.Div(distance, squareRoot);
                bool large = ratio.LowerBound + ratio.UpperBound >= 2;
                double t = Math.Sqrt(Number(large ? c.Div(c.I(1), ratio) : ratio));
                double sine = 2 * t / (1 + t * t), cosine = (1 - t * t) / (1 + t * t);
                double incomplete = sine * Rf(cosine * cosine, cosine * cosine + q * sine * sine, 1);
                real = large ? incomplete / (2 * h) : period / 2 - incomplete / (2 * h);
                // For positive z near zero, wp'(z)<0, hence the lower branch has positive tail integral.
                if (y.Sign > 0) real = period - real;
            }
            if (double.IsNaN(real) || double.IsInfinity(real) || real < -period * 1e-10 || real > period * (1 + 1e-10))
                throw new ArithmeticException("Logarithm fell outside its fundamental interval.");
            real = real < 0 ? 0 : real >= period ? 0 : real;
            return new RealEllipticLogarithmResult(e, real, component == 0 ? 0 : Number(periods.SecondPeriodImaginaryPart) / 2, component, period);
        }
    }
}

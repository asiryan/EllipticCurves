using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Certified Faltings height -log(A)/2, where A is the Neron period-lattice area
        /// of the global minimal model. Uses the LMFDB normalization and is invariant under Q-model changes.</summary>
        public RealEnclosure FaltingsHeight(RealComputationOptions options = null, CancellationToken cancellationToken = default)
            => ComputeFaltingsHeight(false, options, cancellationToken);

        /// <summary>Certified stable Faltings height: -log(A)/2 + log(denominator(j)/|Delta|)/12,
        /// with A and Delta on the global minimal model. Uses the LMFDB normalization.</summary>
        public RealEnclosure StableFaltingsHeight(RealComputationOptions options = null, CancellationToken cancellationToken = default)
            => ComputeFaltingsHeight(true, options, cancellationToken);

        private RealEnclosure ComputeFaltingsHeight(bool stable, RealComputationOptions options, CancellationToken token)
        {
            var c = RealContext(options, token);
            var periods = GetPeriods(c.Options, token);
            var result = c.Div(c.Neg(c.Log(periods.Area)), c.I(2));
            if (stable)
            {
                var e = periods.MinimalModel;
                var ratio = new BigRational(e.JInvariant.Den, BigInteger.Abs(e.Discriminant.Num));
                result = c.Add(result, c.Div(c.Log(c.I(ratio)), c.I(12)));
            }
            return c.Finish(result);
        }
    }
}

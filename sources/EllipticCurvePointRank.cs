using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Certify a rank lower bound from supplied rational points using exact good-reduction characters.
        /// No factorization, minimalization, point search, BSD or GRH is required.
        /// A failure to certify all points does not prove dependence.</summary>
        public PointRankCertificate GetRankLowerBound(IReadOnlyList<EllipticCurvePoint> points,
            int reductionPrimeBound = 1009, CancellationToken cancellationToken = default)
        {
            if (points == null) throw new ArgumentNullException(nameof(points));
            var options = new RankComputationOptions { ReductionPrimeBound = reductionPrimeBound, SearchBound = 0 }.Snapshot();
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("A singular curve has no Mordell-Weil rank.");
            var input = new EllipticCurvePoint[points.Count];
            for (int i = 0; i < input.Length; i++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                input[i] = points[i];
                if (!IsOnCurve(input[i])) throw new ArgumentException("A supplied point is not on the curve.", nameof(points));
            }
            // Clear coefficient denominators by x'=d^2*x, y'=d^3*y. Good reduction
            // of this integral model suffices; computing a minimal model is unnecessary.
            BigInteger d = 1;
            foreach (var coefficient in new[] { A1, A2, A3, A4, A6 })
                d = d / BigInteger.GreatestCommonDivisor(d, coefficient.Den) * coefficient.Den;
            var map = ChangeModel(new BigRational(1, d), 0, 0, 0);
            var certificate = new RationalPointRank(map.Target, null, new DescentBudget(options, cancellationToken));
            foreach (var point in input) certificate.Add(map.Map(point));
            return new PointRankCertificate(input.Length, certificate);
        }
    }
}

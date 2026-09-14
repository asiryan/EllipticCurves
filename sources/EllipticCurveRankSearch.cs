using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>
        /// Search for rational points from the equation alone and certify a lower rank bound.
        /// Uses exact modular sieving, without databases, factorization, minimalization or full descent.
        /// A zero result only means this search did not certify positive rank. A finite box and finite
        /// set of reduction characters can miss independent points. Cancellation throws; time/work
        /// limits return the points and bound already established. Progress reports bound increases.
        /// </summary>
        public RankLowerBoundSearchResult SearchRankLowerBound(RankLowerBoundSearchOptions options = null,
            CancellationToken cancellationToken = default, IProgress<RankLowerBoundSearchResult> progress = null)
        {
            options = (options ?? new RankLowerBoundSearchOptions()).Snapshot();
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("A singular curve has no Mordell-Weil rank.");
            var watch = Stopwatch.StartNew();
            BigInteger scale = 1;
            foreach (var a in new[] { A1, A2, A3, A4, A6 })
                scale = scale / BigInteger.GreatestCommonDivisor(scale, a.Den) * a.Den;
            var map = ChangeModel(new BigRational(1, scale), 0, 0, 0);
            var budget = new DescentBudget(new RankComputationOptions
            { ReductionPrimeBound = options.ReductionPrimeBound, SearchBound = 0 }, cancellationToken);
            var certificate = new RationalPointRank(map.Target, null, budget);
            var points = new List<EllipticCurvePoint>();
            var search = new SievedRationalPointSearch(map.Target, options, watch, cancellationToken);
            RankLowerBoundSearchResult Result(RankLowerBoundSearchStopReason reason) =>
                new RankLowerBoundSearchResult(new PointRankCertificate(points.Count, certificate),
                    points, reason, search.SquareTests, search.Blocks, watch.Elapsed, scale);
            var stop = search.Run(point =>
            {
                int before = certificate.LowerBound;
                certificate.Add(point);
                points.Add(map.MapBack(point));
                bool target = options.TargetLowerBound.HasValue && certificate.LowerBound >= options.TargetLowerBound.Value;
                if (certificate.LowerBound > before && progress != null)
                    progress.Report(Result(target ? RankLowerBoundSearchStopReason.TargetReached : RankLowerBoundSearchStopReason.Searching));
                if (target) return RankLowerBoundSearchStopReason.TargetReached;
                if (points.Count >= options.MaxPoints) return RankLowerBoundSearchStopReason.PointLimit;
                return (RankLowerBoundSearchStopReason?)null;
            });
            return Result(stop);
        }
    }
}

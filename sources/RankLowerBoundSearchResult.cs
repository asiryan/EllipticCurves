using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>Why a bounded lower-rank search returned. None of these states asserts an exact rank.</summary>
    public enum RankLowerBoundSearchStopReason
    {
        /// <summary>A progress snapshot; the search has not stopped.</summary>
        Searching,
        /// <summary>The requested coordinate box was exhausted.</summary>
        SearchBoxExhausted,
        /// <summary>The elapsed-time budget was reached.</summary>
        TimeLimit,
        /// <summary>The exact square-test allowance was reached.</summary>
        SquareTestLimit,
        /// <summary>The retained-point limit was reached.</summary>
        PointLimit,
        /// <summary>The requested lower rank bound has been proved.</summary>
        TargetReached
    }

    /// <summary>Proved lower bound and point witnesses returned by an equation-only search.</summary>
    public sealed class RankLowerBoundSearchResult
    {
        /// <summary>A proved lower bound. Zero does not prove rank zero.</summary>
        public int LowerBound => Certificate.LowerBound;
        /// <summary>Exact certificate for the returned points.</summary>
        public PointRankCertificate Certificate { get; }
        /// <summary>Found points in the original input model, one per pair P,-P. They need not all be independent.</summary>
        public IReadOnlyList<EllipticCurvePoint> Points { get; }
        /// <summary>Reason for stopping; exhausting a finite box does not prove an exact rank.</summary>
        public RankLowerBoundSearchStopReason StopReason { get; }
        /// <summary>Number of exact integer-square tests after the sieve.</summary>
        public long SquareTests { get; }
        /// <summary>Number of sieve blocks prepared, including a final partially processed block.</summary>
        public long SieveBlocks { get; }
        /// <summary>Elapsed time including coefficient normalization and certification.</summary>
        public TimeSpan Elapsed { get; }
        /// <summary>The search uses x'=IntegralScale²*x, y'=IntegralScale³*y to clear coefficient denominators.</summary>
        public BigInteger IntegralScale { get; }

        internal RankLowerBoundSearchResult(PointRankCertificate certificate,
            List<EllipticCurvePoint> points, RankLowerBoundSearchStopReason reason,
            long squareTests, long sieveBlocks, TimeSpan elapsed, BigInteger scale)
        {
            Certificate = certificate;
            Points = Array.AsReadOnly(points.ToArray());
            StopReason = reason;
            SquareTests = squareTests;
            SieveBlocks = sieveBlocks;
            Elapsed = elapsed;
            IntegralScale = scale;
        }
    }
}

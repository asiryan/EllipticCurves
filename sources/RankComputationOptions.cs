using System;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>Search and work limits for unconditional algebraic rank bounds.</summary>
    public sealed class RankComputationOptions
    {
        /// <summary>Coordinate bound for rational points and quartic points. Zero disables point searches.</summary>
        public int SearchBound { get; set; } = 32;
        /// <summary>Maximum square classes per 2-isogeny descent, or covering classes including the identity in general descent.</summary>
        public int MaxSquareClasses { get; set; } = 65536;
        /// <summary>Maximum counted enumeration, polynomial and local-lifting steps. An incomplete descent gives no upper bound.</summary>
        public long MaxDescentWork { get; set; } = 5000000;
        /// <summary>Maximum total point-search steps, on the curve and its coverings. Exhaustion preserves proved lower bounds.</summary>
        public long MaxPointSearchWork { get; set; } = 1000000;
        /// <summary>Good odd primes up to this bound are used to certify lower bounds by reduction.</summary>
        public int ReductionPrimeBound { get; set; } = 101;
        /// <summary>Use the general binary-quartic descent even when a rational 2-isogeny is available.</summary>
        public bool PreferGeneralTwoDescent { get; set; }

        internal RankComputationOptions Snapshot()
        {
            var copy = (RankComputationOptions)MemberwiseClone();
            if (copy.SearchBound < 0 || copy.SearchBound == int.MaxValue) throw new ArgumentOutOfRangeException(nameof(SearchBound));
            if (copy.MaxSquareClasses < 2) throw new ArgumentOutOfRangeException(nameof(MaxSquareClasses));
            if (copy.MaxDescentWork < 0) throw new ArgumentOutOfRangeException(nameof(MaxDescentWork));
            if (copy.MaxPointSearchWork < 0) throw new ArgumentOutOfRangeException(nameof(MaxPointSearchWork));
            if (copy.ReductionPrimeBound < 3 || copy.ReductionPrimeBound > 10000) throw new ArgumentOutOfRangeException(nameof(ReductionPrimeBound));
            return copy;
        }
    }

    internal sealed class DescentLimitException : Exception
    {
        internal DescentLimitException(string reason) : base(reason) { }
    }

    internal sealed class DescentBudget
    {
        internal readonly RankComputationOptions Options;
        internal readonly CancellationToken Token;
        internal long Work { get; private set; }
        internal long PointWork { get; private set; }
        internal bool PointSearchExhausted { get; private set; }
        internal DescentBudget(RankComputationOptions options, CancellationToken token)
        { Options = options; Token = token; }
        internal void Step()
        {
            Token.ThrowIfCancellationRequested();
            if (Work >= Options.MaxDescentWork) throw new DescentLimitException("MaxDescentWork was reached; the upper-bound computation is incomplete.");
            Work++;
        }
        internal bool PointStep()
        {
            Token.ThrowIfCancellationRequested();
            if (PointWork >= Options.MaxPointSearchWork) { PointSearchExhausted = true; return false; }
            PointWork++;
            return true;
        }
    }
}

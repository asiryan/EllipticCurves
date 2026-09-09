namespace EllipticCurves
{
    /// <summary>Unconditional bounds for the Mordell-Weil rank over Q.</summary>
    public sealed class RankBounds
    {
        internal RankBounds(int lower, int? upper, bool twoIsogeny, bool general = false,
            int? selmerDimension = null, string reason = null, long descentWork = 0, long pointWork = 0)
        {
            if (lower < 0 || (upper.HasValue && upper.Value < lower))
                throw new System.ArgumentException("Inconsistent rank bounds.");
            LowerBound = lower;
            UpperBound = upper;
            UsedTwoIsogenyDescent = twoIsogeny;
            UsedGeneralTwoDescent = general;
            TwoSelmerDimension = selmerDimension;
            Reason = reason ?? "";
            DescentWork = descentWork;
            PointSearchWork = pointWork;
        }

        /// <summary>A proved lower bound; zero does not assert that the rank is zero.</summary>
        public int LowerBound { get; }

        /// <summary>A proved upper bound, or null when no upper bound was computed.</summary>
        public int? UpperBound { get; }

        /// <summary>True precisely when both proved bounds coincide.</summary>
        public bool IsExact => UpperBound.HasValue && LowerBound == UpperBound.Value;

        /// <summary>The exact algebraic rank when proved, otherwise null.</summary>
        public int? ExactRank => IsExact ? LowerBound : (int?)null;

        /// <summary>Whether rational 2-torsion allowed descent on a pair of 2-isogenous curves.</summary>
        public bool UsedTwoIsogenyDescent { get; }

        /// <summary>Whether general binary-quartic 2-descent was attempted.</summary>
        public bool UsedGeneralTwoDescent { get; }

        /// <summary>Dimension of Sel_2(E/Q), only when general 2-descent completed.</summary>
        public int? TwoSelmerDimension { get; }

        /// <summary>Explanation of the bounds and any exhausted work limit.</summary>
        public string Reason { get; }

        /// <summary>Counted descent steps; integer factorization is not included.</summary>
        public long DescentWork { get; }

        /// <summary>Counted rational-point search steps.</summary>
        public long PointSearchWork { get; }

        /// <inheritdoc />
        public override string ToString() => IsExact ? LowerBound.ToString()
            : $"[{LowerBound}, {(UpperBound.HasValue ? UpperBound.Value.ToString() : "unknown")}]";
    }
}

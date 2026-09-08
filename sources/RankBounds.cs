namespace EllipticCurves
{
    /// <summary>Unconditional bounds for the Mordell-Weil rank over Q.</summary>
    public sealed class RankBounds
    {
        internal RankBounds(int lower, int? upper, bool twoIsogeny)
        {
            LowerBound = lower;
            UpperBound = upper;
            UsedTwoIsogenyDescent = twoIsogeny;
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

        /// <inheritdoc />
        public override string ToString() => IsExact ? LowerBound.ToString()
            : $"[{LowerBound}, {(UpperBound.HasValue ? UpperBound.Value.ToString() : "unknown")}]";
    }
}

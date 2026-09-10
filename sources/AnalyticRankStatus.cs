namespace EllipticCurves
{
    /// <summary>Whether an analytic rank was proved, estimated numerically, or left unresolved.</summary>
    public enum AnalyticRankStatus
    {
        /// <summary>No rank is returned: a computation limit or numerical ambiguity was encountered.</summary>
        Inconclusive,
        /// <summary>A numerical estimate; small derivatives have not been proved to vanish.</summary>
        NumericalEstimate,
        /// <summary>Rank zero or one proved by an outward-rounded interval excluding zero.</summary>
        Certified
    }
}

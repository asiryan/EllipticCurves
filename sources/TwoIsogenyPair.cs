namespace EllipticCurves
{
    /// <summary>A degree-two isogeny and its dual, whose compositions are multiplication by two on both models.</summary>
    public sealed class TwoIsogenyPair
    {
        /// <summary>Isogeny from the input curve to its quotient.</summary>
        public RationalIsogeny Forward { get; }
        /// <summary>Dual isogeny back to the exact input model.</summary>
        public RationalIsogeny Dual { get; }
        internal TwoIsogenyPair(RationalIsogeny forward, RationalIsogeny dual) { Forward = forward; Dual = dual; }
    }
}

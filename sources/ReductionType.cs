namespace EllipticCurves
{
    /// <summary>Reduction of an elliptic curve at a finite prime.</summary>
    public enum ReductionType
    {
        /// <summary>Smooth reduction.</summary>
        Good,
        /// <summary>Split nodal reduction.</summary>
        SplitMultiplicative,
        /// <summary>Nonsplit nodal reduction.</summary>
        NonSplitMultiplicative,
        /// <summary>Additive reduction.</summary>
        Additive
    }
}

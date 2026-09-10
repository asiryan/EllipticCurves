using System.Numerics;

namespace EllipticCurves
{
    /// <summary>Exact local invariants, always for a minimal model.</summary>
    public sealed class LocalReductionData
    {
        /// <summary>Residue characteristic.</summary>
        public BigInteger Prime { get; }
        /// <summary>Valuation of the minimal discriminant.</summary>
        public int DiscriminantValuation { get; }
        /// <summary>Valuation of the conductor.</summary>
        public int ConductorValuation { get; }
        /// <summary>Valuation of the denominator of j.</summary>
        public int JDenominatorValuation { get; }
        /// <summary>Kodaira symbol in ASCII, e.g. I0, I5, IV, I3*.</summary>
        public string KodairaSymbol { get; }
        /// <summary>Reduction type including splitting for multiplicative reduction.</summary>
        public ReductionType ReductionType { get; }
        /// <summary>Local Tamagawa number.</summary>
        public int TamagawaNumber { get; }
        /// <summary>Local root number, +1 or -1.</summary>
        public int RootNumber { get; }
        internal LocalReductionData(BigInteger p, int n, int f, int j, string symbol, ReductionType type, int cp, int root)
        { Prime = p; DiscriminantValuation = n; ConductorValuation = f; JDenominatorValuation = j; KodairaSymbol = symbol; ReductionType = type; TamagawaNumber = cp; RootNumber = root; }
    }
}

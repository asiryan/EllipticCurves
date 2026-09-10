using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>Exact saturation of a subgroup modulo the entire rational torsion subgroup.</summary>
    public sealed class SaturationResult
    {
        /// <summary>Generators of the resulting subgroup, modulo rational torsion, on the input model.</summary>
        public IReadOnlyList<EllipticCurvePoint> Generators { get; }
        /// <summary>Primes at which the returned subgroup was proved saturated.</summary>
        public IReadOnlyList<int> CertifiedPrimes { get; }
        /// <summary>Requested primes whose saturation has not been proved.</summary>
        public IReadOnlyList<int> UnresolvedPrimes { get; }
        /// <summary>True only if saturation was proved at every requested prime. This does not assert a full Mordell-Weil basis.</summary>
        public bool IsComplete => IndependenceCertified && UnresolvedPrimes.Count == 0;
        /// <summary>Whether independence of the input generators modulo torsion was proved.</summary>
        public bool IndependenceCertified { get; }
        /// <summary>Exact enlargement index relative to the input subgroup modulo torsion.</summary>
        public BigInteger IndexGain { get; }
        /// <summary>Counted work, excluding minimalization, torsion and height preparation.</summary>
        public long Work { get; }
        /// <summary>Completion or limit explanation.</summary>
        public string Reason { get; }
        internal SaturationResult(IEnumerable<EllipticCurvePoint> generators, IEnumerable<int> certified, IEnumerable<int> unresolved, bool independent, BigInteger index, long work, string reason)
        { Generators = Array.AsReadOnly(generators.ToArray()); CertifiedPrimes = Array.AsReadOnly(certified.ToArray()); UnresolvedPrimes = Array.AsReadOnly(unresolved.ToArray()); IndependenceCertified = independent; IndexGain = index; Work = work; Reason = reason; }
    }
}

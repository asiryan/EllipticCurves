using System.Numerics;

namespace EllipticCurves.Explorer.Computations;

internal sealed record ConductorCalculationResult(BigInteger Conductor,
    IReadOnlyDictionary<BigInteger, int> Factorization);

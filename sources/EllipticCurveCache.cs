using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        private readonly CachedComputation<EllipticCurveQ> minimalModelCache = new CachedComputation<EllipticCurveQ>();
        private readonly CachedComputation<IReadOnlyDictionary<BigInteger, int>> minimalDiscriminantCache
            = new CachedComputation<IReadOnlyDictionary<BigInteger, int>>();

        private EllipticCurveQ GetGlobalMinimalModelCore(int maxWorkers, CancellationToken token)
            => minimalModelCache.Get(() =>
            {
                var model = ComputeGlobalMinimalModel(maxWorkers, token);
                if (Equals(model)) return this;
                // A computed minimal model is already normalized. Calls through
                // the returned model must share its factorization with the input.
                model.minimalModelCache.Get(() => model, token);
                return model;
            }, token);

        // The absolute minimal discriminant, not the conductor: their prime
        // divisors agree, but the exponents need not. No process-wide cache.
        internal IReadOnlyDictionary<BigInteger, int> GetMinimalDiscriminantFactorization(CancellationToken token)
            => GetMinimalDiscriminantFactorization(0, token);

        internal IReadOnlyDictionary<BigInteger, int> GetMinimalDiscriminantFactorization(int maxWorkers, CancellationToken token)
        {
            var model = GetGlobalMinimalModelCore(maxWorkers, token);
            return model.minimalDiscriminantCache.Get(() =>
                new ReadOnlyDictionary<BigInteger, int>(new SortedDictionary<BigInteger, int>(
                    NativeNumberTheory.Factor(model.Discriminant.Num, maxWorkers, token))), token);
        }
    }
}

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

        private EllipticCurveQ GetGlobalMinimalModelCore(CancellationToken token, int maxWorkers)
            => minimalModelCache.Get(() =>
            {
                var model = ComputeGlobalMinimalModel(token, maxWorkers);
                if (Equals(model)) return this;
                // A computed minimal model is already normalized. Calls through
                // the returned model must share its factorization with the input.
                model.minimalModelCache.Get(() => model, token);
                return model;
            }, token);

        // The absolute minimal discriminant, not the conductor: their prime
        // divisors agree, but the exponents need not. No process-wide cache.
        internal IReadOnlyDictionary<BigInteger, int> GetMinimalDiscriminantFactorization(CancellationToken token, int maxWorkers = 0)
        {
            var model = GetGlobalMinimalModelCore(token, maxWorkers);
            return model.minimalDiscriminantCache.Get(() =>
                new ReadOnlyDictionary<BigInteger, int>(new SortedDictionary<BigInteger, int>(
                    NativeNumberTheory.Factor(model.Discriminant.Num, token, maxWorkers))), token);
        }
    }
}

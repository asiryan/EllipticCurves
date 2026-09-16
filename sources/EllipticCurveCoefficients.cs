using System;
using System.Collections.Generic;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Return a_p on a minimal model. At bad primes this is the local L-series coefficient (0 or +/-1).
        /// Uses direct counting, with work proportional to p; not SEA.</summary>
        public long GetFrobeniusTrace(int prime, long maxPointCountingWork = 20000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            if (maxPointCountingWork < 0) throw new ArgumentOutOfRangeException(nameof(maxPointCountingWork));
            if (prime > maxPointCountingWork) throw new ArithmeticException("Point-counting work limit reached.");
            return AnalyticCoefficients.Trace(GetGlobalMinimalModel(cancellationToken), prime, cancellationToken);
        }

        /// <summary>Return a_0=0, a_1,...,a_count, indexed by n. Coefficients are exact integers.
        /// A work limit throws; an incomplete coefficient list is never returned.</summary>
        public IReadOnlyList<long> GetFourierCoefficients(int count, long maxPointCountingWork = 20000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (count < 0 || count == int.MaxValue) throw new ArgumentOutOfRangeException(nameof(count));
            if (maxPointCountingWork < 0) throw new ArgumentOutOfRangeException(nameof(maxPointCountingWork));
            if (IsSingular) throw new InvalidOperationException("L-series coefficients require a nonsingular curve.");
            if (count <= 1) return Array.AsReadOnly(count == 0 ? new long[] { 0 } : new long[] { 0, 1 });
            if (count > maxPointCountingWork) throw new ArithmeticException("Coefficient allocation exceeds the work limit.");
            var result = AnalyticCoefficients.Compute(GetGlobalMinimalModel(cancellationToken), count, maxPointCountingWork, cancellationToken);
            if (result == null) throw new ArithmeticException("Point-counting work limit reached.");
            return Array.AsReadOnly(result);
        }

        /// <summary>Return a_n for n >= 1 using the prime factorization of n, prime-power recurrences and multiplicativity.
        /// The work limit bounds the sum of distinct prime divisors of n; no coefficient list is allocated.</summary>
        public long GetFourierCoefficient(int index, long maxPointCountingWork = 20000000, CancellationToken cancellationToken = default)
        {
            if (index < 1) throw new ArgumentOutOfRangeException(nameof(index));
            cancellationToken.ThrowIfCancellationRequested();
            if (maxPointCountingWork < 0) throw new ArgumentOutOfRangeException(nameof(maxPointCountingWork));
            if (IsSingular) throw new InvalidOperationException("L-series coefficients require a nonsingular curve.");
            if (index == 1) return 1;

            var factors = NativeNumberTheory.Factor(index, cancellationToken);
            long work = 0;
            foreach (var prime in factors.Keys) work += (int)prime;
            if (work > maxPointCountingWork) throw new ArithmeticException("Point-counting work limit reached.");

            var model = GetGlobalMinimalModel(cancellationToken);
            long result = 1;
            foreach (var factor in factors)
            {
                int prime = (int)factor.Key;
                long trace = AnalyticCoefficients.Trace(model, prime, cancellationToken);
                long previous = 1, current = trace;
                long correction = model.Discriminant.Num % prime == 0 ? 0 : prime;
                for (int power = 2; power <= factor.Value; power++)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    long next = checked(trace * current - correction * previous);
                    previous = current; current = next;
                }
                result = checked(result * current);
            }
            return result;
        }

        /// <summary>Count E(F_p), including infinity, on a global minimal model. Rejects bad reduction.</summary>
        public long CountPoints(int prime, long maxPointCountingWork = 20000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            if (maxPointCountingWork < 0) throw new ArgumentOutOfRangeException(nameof(maxPointCountingWork));
            if (prime > maxPointCountingWork) throw new ArithmeticException("Point-counting work limit reached.");
            var minimal = GetGlobalMinimalModel(cancellationToken);
            if (minimal.Discriminant.Num % prime == 0) throw new ArgumentException("The curve has bad reduction at this prime.", nameof(prime));
            return (long)prime + 1 - AnalyticCoefficients.Trace(minimal, prime, cancellationToken);
        }
    }
}

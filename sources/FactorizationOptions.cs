using System;

namespace EllipticCurves
{
    /// <summary>Execution settings for certified integer factorization in native arithmetic.</summary>
    public sealed class FactorizationOptions
    {
        /// <summary>Maximum sieve workers. Zero selects automatically by input size and available CPUs;
        /// one runs sequentially. Small inputs may use fewer workers. Primality proofs remain exact.</summary>
        public int MaxDegreeOfParallelism { get; set; }

        internal int WorkerLimit()
        {
            int value = MaxDegreeOfParallelism;
            if (value < 0) throw new ArgumentOutOfRangeException(nameof(MaxDegreeOfParallelism));
            return value;
        }
    }
}

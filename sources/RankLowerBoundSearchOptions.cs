using System;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>Bounds for an equation-only rational point search. No upper rank bound is computed.</summary>
    public sealed class RankLowerBoundSearchOptions
    {
        /// <summary>Search x=m/d² with |m-NumeratorCenter| at most this bound on the integral search model.</summary>
        public BigInteger NumeratorRadius { get; set; } = 2000000000;
        /// <summary>Center of the numerator interval; arbitrary-size integers are supported.</summary>
        public BigInteger NumeratorCenter { get; set; } = BigInteger.Zero;
        /// <summary>Search denominators d² for d=1,...,this bound, in that order.</summary>
        public int DenominatorRootBound { get; set; } = 8;
        /// <summary>Elapsed-time budget. Checked between sieve operations and exact point checks; not a hard deadline.</summary>
        public TimeSpan TimeLimit { get; set; } = TimeSpan.FromSeconds(10);
        /// <summary>Maximum exact square tests after modular sieving.</summary>
        public long MaxSquareTests { get; set; } = 1000000;
        /// <summary>Maximum retained points (one from each pair P,-P).</summary>
        public int MaxPoints { get; set; } = 1024;
        /// <summary>Stop once this lower bound has been proved; null means no rank target.</summary>
        public int? TargetLowerBound { get; set; }
        /// <summary>Good odd primes up to this bound are used for exact independence certificates.</summary>
        public int ReductionPrimeBound { get; set; } = 1009;

        internal RankLowerBoundSearchOptions Snapshot()
        {
            var copy = (RankLowerBoundSearchOptions)MemberwiseClone();
            if (copy.NumeratorRadius.Sign < 0) throw new ArgumentOutOfRangeException(nameof(NumeratorRadius));
            if (copy.DenominatorRootBound < 1 || copy.DenominatorRootBound == int.MaxValue)
                throw new ArgumentOutOfRangeException(nameof(DenominatorRootBound));
            if (copy.TimeLimit <= TimeSpan.Zero) throw new ArgumentOutOfRangeException(nameof(TimeLimit));
            if (copy.MaxSquareTests < 0) throw new ArgumentOutOfRangeException(nameof(MaxSquareTests));
            if (copy.MaxPoints < 1) throw new ArgumentOutOfRangeException(nameof(MaxPoints));
            if (copy.TargetLowerBound.HasValue && copy.TargetLowerBound.Value < 1)
                throw new ArgumentOutOfRangeException(nameof(TargetLowerBound));
            if (copy.ReductionPrimeBound < 3 || copy.ReductionPrimeBound > 10000)
                throw new ArgumentOutOfRangeException(nameof(ReductionPrimeBound));
            return copy;
        }
    }
}

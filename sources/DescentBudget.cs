using System.Threading;

namespace EllipticCurves
{
    internal sealed class DescentBudget
    {
        internal readonly RankComputationOptions Options;
        internal readonly CancellationToken Token;
        private sealed class Counters
        {
            internal long Work, PointWork;
            internal int PointSearchExhausted;
        }
        private readonly Counters counters;
        internal long Work => Interlocked.Read(ref counters.Work);
        internal long PointWork => Interlocked.Read(ref counters.PointWork);
        internal bool PointSearchExhausted => Volatile.Read(ref counters.PointSearchExhausted) != 0;
        internal DescentBudget(RankComputationOptions options, CancellationToken token)
            : this(options, token, new Counters()) { }
        private DescentBudget(RankComputationOptions options, CancellationToken token, Counters counters)
        { Options = options; Token = token; this.counters = counters; }

        // A worker may use a linked stop token, but never receives a fresh work allowance.
        internal DescentBudget WithToken(CancellationToken token) => new DescentBudget(Options, token, counters);

        internal void Step()
        {
            Token.ThrowIfCancellationRequested();
            if (!Take(ref counters.Work, Options.MaxDescentWork))
                throw new DescentLimitException("MaxDescentWork was reached; the upper-bound computation is incomplete.");
        }
        internal bool PointStep()
        {
            Token.ThrowIfCancellationRequested();
            if (Take(ref counters.PointWork, Options.MaxPointSearchWork)) return true;
            Volatile.Write(ref counters.PointSearchExhausted, 1);
            return false;
        }

        private bool Take(ref long counter, long limit)
        {
            if (Options.MaxDegreeOfParallelism == 1)
            {
                if (counter >= limit) return false;
                counter++;
                return true;
            }
            // Reserve before doing the work. Neither contention nor long.MaxValue can
            // overshoot/wrap the allowance, and failed reservations do not count as work.
            while (true)
            {
                Token.ThrowIfCancellationRequested();
                long current = Interlocked.Read(ref counter);
                if (current >= limit) return false;
                if (Interlocked.CompareExchange(ref counter, current + 1, current) == current) return true;
            }
        }
    }
}

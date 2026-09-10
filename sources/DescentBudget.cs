using System.Threading;

namespace EllipticCurves
{
    internal sealed class DescentBudget
    {
        internal readonly RankComputationOptions Options;
        internal readonly CancellationToken Token;
        internal long Work { get; private set; }
        internal long PointWork { get; private set; }
        internal bool PointSearchExhausted { get; private set; }
        internal DescentBudget(RankComputationOptions options, CancellationToken token)
        { Options = options; Token = token; }
        internal void Step()
        {
            Token.ThrowIfCancellationRequested();
            if (Work >= Options.MaxDescentWork) throw new DescentLimitException("MaxDescentWork was reached; the upper-bound computation is incomplete.");
            Work++;
        }
        internal bool PointStep()
        {
            Token.ThrowIfCancellationRequested();
            if (PointWork >= Options.MaxPointSearchWork) { PointSearchExhausted = true; return false; }
            PointWork++;
            return true;
        }
    }
}

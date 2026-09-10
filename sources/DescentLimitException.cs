using System;

namespace EllipticCurves
{
    internal sealed class DescentLimitException : Exception
    {
        internal DescentLimitException(string reason) : base(reason) { }
    }
}

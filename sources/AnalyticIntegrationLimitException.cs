using System;

namespace EllipticCurves
{
    internal sealed class AnalyticIntegrationLimitException : Exception
    {
        internal AnalyticIntegrationLimitException(string message) : base(message) { }
    }
}

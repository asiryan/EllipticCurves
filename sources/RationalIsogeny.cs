using System;
using System.Collections.Generic;
using System.Linq;

namespace EllipticCurves
{
    /// <summary>An exact Q-isogeny with a pointwise rational kernel. Maps are evaluated using exact rational arithmetic.</summary>
    public sealed class RationalIsogeny
    {
        private readonly Func<EllipticCurvePoint, EllipticCurvePoint> map;
        /// <summary>Domain model.</summary>
        public EllipticCurveQ Source { get; }
        /// <summary>Codomain model.</summary>
        public EllipticCurveQ Target { get; }
        /// <summary>Degree of the isogeny, including degree one for the trivial kernel.</summary>
        public int Degree => Kernel.Count;
        /// <summary>All kernel points on Source, including infinity.</summary>
        public IReadOnlyList<EllipticCurvePoint> Kernel { get; }
        internal RationalIsogeny(EllipticCurveQ source, EllipticCurveQ target,
            IEnumerable<EllipticCurvePoint> kernel, Func<EllipticCurvePoint, EllipticCurvePoint> map)
        { Source = source; Target = target; Kernel = Array.AsReadOnly(kernel.ToArray()); this.map = map; }
        /// <summary>Evaluate the isogeny. Kernel points map to infinity; points outside Source are rejected.</summary>
        public EllipticCurvePoint Map(EllipticCurvePoint point)
        {
            if (!Source.IsOnCurve(point)) throw new ArgumentException("Point is not on the source curve.", nameof(point));
            var result = map(point);
            if (!Target.IsOnCurve(result)) throw new ArithmeticException("Isogeny image is not on the target curve.");
            return result;
        }
    }
}

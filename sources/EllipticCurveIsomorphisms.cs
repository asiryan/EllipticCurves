using System;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Construct a coordinate change with its exact point maps.</summary>
        public WeierstrassIsomorphism ChangeModel(BigRational u, BigRational r, BigRational s, BigRational t)
            => new WeierstrassIsomorphism(this, u, r, s, t);

        /// <summary>Find a rational isomorphism to the given model, if one exists.</summary>
        public bool TryGetIsomorphism(EllipticCurveQ target, out WeierstrassIsomorphism isomorphism)
        {
            if (target == null) throw new ArgumentNullException(nameof(target));
            if (IsSingular || target.IsSingular) throw new InvalidOperationException("Singular models are not supported.");
            isomorphism = null;
            if (!InternalMath.IsQIsomorphic(C4, C6, Discriminant, target.C4, target.C6, target.Discriminant, out var u)) return false;
            var s = (u * target.A1 - A1) / 2;
            var r = (u * u * target.A2 - A2 + s * A1 + s * s) / 3;
            var t = (BigRational.Pow(u, 3) * target.A3 - A3 - r * A1) / 2;
            var map = ChangeModel(u, r, s, t);
            if (!map.Target.Equals(target)) return false;
            isomorphism = map;
            return true;
        }

        /// <summary>Return the reduced global minimal model together with both point maps.</summary>
        public WeierstrassIsomorphism GetMinimalModelIsomorphism(CancellationToken cancellationToken = default)
        {
            var minimal = GetGlobalMinimalModel(cancellationToken);
            if (!TryGetIsomorphism(minimal, out var map)) throw new InvalidOperationException("Minimal model is not isomorphic to the input.");
            return map;
        }
    }
}

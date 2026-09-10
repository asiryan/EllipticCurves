using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Reduce the global minimal model over F_p. Bad reduction is rejected.</summary>
        public EllipticCurveFp ReduceModuloPrime(BigInteger prime, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            var e = GetGlobalMinimalModel(cancellationToken);
            return new EllipticCurveFp(prime, e.A1.Num, e.A2.Num, e.A3.Num, e.A4.Num, e.A6.Num, cancellationToken);
        }

        /// <summary>Map a rational point to the global minimal model and reduce at a good prime.
        /// A pole of its minimal x-coordinate reduces to infinity.</summary>
        public EllipticCurvePointFp ReducePointModuloPrime(EllipticCurvePoint point, BigInteger prime, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!IsOnCurve(point)) throw new ArgumentException("Point is not on the curve.", nameof(point));
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            var map = GetMinimalModelIsomorphism(cancellationToken); var e = map.Target;
            var reduction = new EllipticCurveFp(prime, e.A1.Num, e.A2.Num, e.A3.Num, e.A4.Num, e.A6.Num, cancellationToken);
            var p = map.Map(point);
            if (p.IsInfinity || p.X.Den % prime == 0) return EllipticCurvePointFp.Infinity;
            BigInteger Reduce(BigRational value) => NativeNumberTheory.Mod(value.Num * BigInteger.ModPow(value.Den % prime, prime - 2, prime), prime);
            return reduction.CreatePoint(Reduce(p.X), Reduce(p.Y));
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Construct a Velu isogeny from generators of a finite subgroup of E(Q).
        /// All kernel points must be rational; this does not discover nonrational Galois-stable kernels.
        /// The target is a short model and is not necessarily minimal.</summary>
        public RationalIsogeny CreateIsogeny(IReadOnlyList<EllipticCurvePoint> kernelGenerators, CancellationToken cancellationToken = default)
        {
            if (kernelGenerators == null) throw new ArgumentNullException(nameof(kernelGenerators));
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("Isogenies require a nonsingular curve.");
            var generators = kernelGenerators.ToArray();
            foreach (var generator in generators)
                if (!IsOnCurve(generator)) throw new ArgumentException("Kernel generator is not on the curve.", nameof(kernelGenerators));
            var kernel = new HashSet<EllipticCurvePoint> { EllipticCurvePoint.Infinity };
            foreach (var generator in generators)
            {
                cancellationToken.ThrowIfCancellationRequested();
                // Over Q a torsion point has order <= 12 (Mazur). This exact test
                // avoids factoring the discriminant just to validate a supplied kernel.
                var multiples = new List<EllipticCurvePoint>(); var p = EllipticCurvePoint.Infinity;
                do
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    multiples.Add(p); p = Add(p, generator);
                    if (multiples.Count == 12 && !p.IsInfinity) throw new ArgumentException("Kernel generator has infinite order.", nameof(kernelGenerators));
                } while (!p.IsInfinity);
                var previous = kernel.ToArray();
                foreach (var q in previous) foreach (var t in multiples) kernel.Add(Add(q, t));
            }
            return VeluFromKernel(kernel, cancellationToken);
        }

        private RationalIsogeny VeluFromKernel(IEnumerable<EllipticCurvePoint> kernel, CancellationToken token)
        {
            var sourceKernel = kernel.Distinct().OrderBy(p => p.IsInfinity ? 0 : 1).ThenBy(p => p.X).ThenBy(p => p.Y).ToArray();
            var change = GetShortModelIsomorphism(); var shortModel = change.Target;
            var shortKernel = sourceKernel.Where(p => !p.IsInfinity).Select(change.Map).ToArray();
            BigRational t = 0, w = 0;
            foreach (var p in shortKernel)
            {
                token.ThrowIfCancellationRequested();
                var tq = 3 * p.X * p.X + shortModel.A4;
                t += tq; w += 2 * p.Y * p.Y + p.X * tq;
            }
            var target = new EllipticCurveQ(0, 0, 0, shortModel.A4 - 5 * t, shortModel.A6 - 7 * w);
            var lookup = new HashSet<EllipticCurvePoint>(sourceKernel);
            EllipticCurvePoint Evaluate(EllipticCurvePoint point)
            {
                if (lookup.Contains(point)) return EllipticCurvePoint.Infinity;
                var p = change.Map(point); var x = p.X; var y = p.Y;
                foreach (var q in shortKernel)
                { var sum = shortModel.Add(p, q); x += sum.X - q.X; y += sum.Y - q.Y; }
                return new EllipticCurvePoint(x, y);
            }
            return new RationalIsogeny(this, target, sourceKernel, Evaluate);
        }

        /// <summary>Construct a 2-isogeny from a nonzero rational point of order two, together with its dual.</summary>
        public TwoIsogenyPair CreateTwoIsogeny(EllipticCurvePoint kernelPoint, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("Isogenies require a nonsingular curve.");
            if (!IsOnCurve(kernelPoint) || kernelPoint.IsInfinity || !Double(kernelPoint).IsInfinity)
                throw new ArgumentException("A nonzero rational point of order two is required.", nameof(kernelPoint));
            var shortMap = GetShortModelIsomorphism(); var root = shortMap.Map(kernelPoint).X;
            var forward = VeluFromKernel(new[] { EllipticCurvePoint.Infinity, kernelPoint }, cancellationToken);
            var dualKernel = new EllipticCurvePoint(-2 * root, 0);
            var second = forward.Target.VeluFromKernel(new[] { EllipticCurvePoint.Infinity, dualKernel }, cancellationToken);
            var scaling = second.Target.ChangeModel(2, 0, 0, 0);
            if (!scaling.Target.Equals(shortMap.Target)) throw new ArithmeticException("Dual isogeny model mismatch.");
            var dual = new RationalIsogeny(forward.Target, this, second.Kernel,
                p => shortMap.MapBack(scaling.Map(second.Map(p))));
            return new TwoIsogenyPair(forward, dual);
        }

        /// <summary>All 2-isogenies over Q, one for each nonzero rational 2-torsion point, with duals.</summary>
        public IReadOnlyList<TwoIsogenyPair> GetTwoIsogenies(CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (IsSingular) throw new InvalidOperationException("Isogenies require a nonsingular curve.");
            var result = new List<TwoIsogenyPair>();
            foreach (var p in TorsionPoints)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (!p.IsInfinity && Double(p).IsInfinity) result.Add(CreateTwoIsogeny(p, cancellationToken));
            }
            return result.AsReadOnly();
        }
    }
}

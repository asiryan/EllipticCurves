using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class IsogenyTests
{
    [Fact]
    public void KnownVeluQuotientsAndKernels()
    {
        var e = new EllipticCurveQ(0, 0, 0, 0, 1);
        var two = e.CreateIsogeny(new[] { new EllipticCurvePoint(-1, 0) });
        Assert.Equal(new EllipticCurveQ(0, 0, 0, -15, 22), two.Target); Assert.Equal(2, two.Degree);
        var three = e.CreateIsogeny(new[] { new EllipticCurvePoint(0, 1) });
        Assert.Equal(new EllipticCurveQ(0, 0, 0, 0, -27), three.Target); Assert.Equal(3, three.Degree);
        var six = e.CreateIsogeny(new[] { new EllipticCurvePoint(0, 1), new EllipticCurvePoint(-1, 0) });
        Assert.Equal(6, six.Degree);
        foreach (var isogeny in new[] { two, three, six })
        {
            Assert.All(isogeny.Kernel, p => Assert.True(isogeny.Map(p).IsInfinity));
            foreach (var p in e.TorsionPoints) foreach (var q in e.TorsionPoints)
                Assert.Equal(isogeny.Map(e.Add(p, q)), isogeny.Target.Add(isogeny.Map(p), isogeny.Map(q)));
            foreach (int prime in new[] { 5, 7, 11, 13, 17, 19 }) Assert.Equal(e.CountPoints(prime), isogeny.Target.CountPoints(prime));
        }
    }

    [Fact]
    public void DualCompositionsAreDoublingOnBothModels()
    {
        var original = new EllipticCurveQ(0, 0, 0, -25, 0); var p = new EllipticCurvePoint(new BigRational(25, 4), new BigRational(75, 8));
        var change = original.ChangeModel(new BigRational(-3, 2), 7, -2, 5); var e = change.Target; p = change.Map(p);
        var pairs = e.GetTwoIsogenies(); Assert.Equal(3, pairs.Count);
        foreach (var pair in pairs)
        {
            Assert.Equal(e, pair.Dual.Target); Assert.Equal(pair.Forward.Target, pair.Dual.Source);
            for (int n = -4; n <= 4; n++) foreach (var t in e.TorsionPoints)
            {
                var point = e.Add(e.Multiply(p, n), t);
                Assert.Equal(e.Double(point), pair.Dual.Map(pair.Forward.Map(point)));
            }
            foreach (var q in pair.Forward.Target.RationalPoints(30, 4))
                Assert.Equal(pair.Forward.Target.Double(q), pair.Forward.Map(pair.Dual.Map(q)));
        }
    }

    [Fact]
    public void NoncyclicKernelAndTrivialKernelAreSupported()
    {
        var e = new EllipticCurveQ(0, -17, 0, 72, 0);
        var full = e.CreateIsogeny(e.TorsionPoints.ToArray()); Assert.Equal(8, full.Degree);
        Assert.All(e.TorsionPoints, p => Assert.True(full.Map(p).IsInfinity));
        var identity = e.CreateIsogeny(Array.Empty<EllipticCurvePoint>());
        Assert.Equal(1, identity.Degree); Assert.Equal(e.ShortWeierstrass, identity.Target);
        Assert.All(e.TorsionPoints, p => Assert.True(identity.Target.IsOnCurve(identity.Map(p))));
        foreach (int prime in new[] { 5, 7, 11, 13, 17 }) Assert.Equal(e.CountPoints(prime), full.Target.CountPoints(prime));
    }

    [Fact]
    public void InvalidKernelsAndPointsAreRejected()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        Assert.Throws<ArgumentException>(() => e.CreateIsogeny(new[] { p }));
        Assert.Throws<ArgumentException>(() => e.CreateIsogeny(new[] { new EllipticCurvePoint(0, 1) }));
        Assert.Throws<ArgumentException>(() => e.CreateTwoIsogeny(p));
        Assert.Throws<ArgumentException>(() => e.CreateTwoIsogeny(EllipticCurvePoint.Infinity));
        Assert.Empty(e.GetTwoIsogenies());
        Assert.Throws<OperationCanceledException>(() => e.CreateIsogeny(Array.Empty<EllipticCurvePoint>(), new CancellationToken(true)));
    }
}

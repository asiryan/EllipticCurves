using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class PointDivisionTests
{
    [Theory]
    [InlineData(2), InlineData(3), InlineData(4), InlineData(5), InlineData(6), InlineData(-2), InlineData(-3)]
    public void DividesKnownMultiplesExactlyOnChangedModels(int n)
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        var map = e.ChangeModel(new BigRational(-2, 3), 5, -2, 7); e = map.Target; p = map.Map(p);
        var multiple = e.Multiply(p, n);
        Assert.Equal(new[] { p }, e.GetDivisionPoints(multiple, n));
        Assert.True(e.TryDividePoint(multiple, n, out var divided)); Assert.Equal(p, divided);
    }

    [Fact]
    public void AllPreimagesIncludeTorsionTranslatesAndDivisionOfInfinity()
    {
        var e = new EllipticCurveQ(0, 0, 0, -25, 0); var p = new EllipticCurvePoint(new BigRational(25, 4), new BigRational(75, 8));
        var expected = e.TorsionPoints.Select(t => e.Add(p, t)).ToHashSet();
        Assert.Equal(4, expected.Count); Assert.True(expected.SetEquals(e.GetDivisionPoints(e.Double(p), 2)));
        Assert.True(e.TorsionPoints.ToHashSet().SetEquals(e.GetDivisionPoints(EllipticCurvePoint.Infinity, 2)));
        e = new EllipticCurveQ(0, 0, 0, 0, 1); p = new EllipticCurvePoint(0, 1);
        Assert.True(e.TorsionPoints.Where(q => e.Double(q).Equals(p)).ToHashSet().SetEquals(e.GetDivisionPoints(p, 2)));
        Assert.Equal(3, e.GetDivisionPoints(new EllipticCurvePoint(-1, 0), 3).Count);
    }

    [Fact]
    public void NondivisibilityIsDistinctFromResourceExhaustion()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        Assert.False(e.TryDividePoint(p, 2, out var missing)); Assert.True(missing.IsInfinity);
        Assert.Equal(new[] { e.Negate(p) }, e.GetDivisionPoints(p, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetDivisionPoints(p, 0));
        Assert.Throws<ArgumentException>(() => e.GetDivisionPoints(new EllipticCurvePoint(0, 1), 2));
        Assert.Throws<ArithmeticException>(() => e.GetDivisionPoints(e.Double(p), 2, new PointDivisionOptions { MaxWork = 0 }));
        Assert.Throws<ArithmeticException>(() => e.TryDividePoint(p, 100, out _));
        Assert.Throws<OperationCanceledException>(() => e.GetDivisionPoints(p, 1, cancellationToken: new CancellationToken(true)));
    }
}

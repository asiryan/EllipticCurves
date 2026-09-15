using Xunit;

namespace EllipticCurves.Tests;

public class PointRankTests
{
    private static EllipticCurveQ Curve() => new(0, 0, 1, -7, 6);
    private static EllipticCurvePoint[] Points() => new EllipticCurvePoint[]
    {
        new(0, 2), new(1, 0), new(2, 0)
    };

    [Fact]
    public void SuppliedPointsAreCertifiedOnOriginalAndRationalModels()
    {
        var curve = Curve();
        var points = Points();
        var certificate = curve.GetRankLowerBound(points);
        Assert.Equal(3, certificate.LowerBound);
        Assert.Equal(3, certificate.PointCount);
        Assert.True(certificate.IndependenceCertified);
        Assert.Equal(0, certificate.TwoTorsionDimensionUpperBound);
        Assert.True(certificate.GoodPrimeCount > 0);
        Assert.NotNull(certificate.NoTwoTorsionPrime);

        var map = curve.ChangeModel(6, new BigRational(1, 7), new BigRational(1, 5), new BigRational(1, 3));
        var mapped = map.Target.GetRankLowerBound(points.Select(map.Map).ToArray());
        Assert.Equal(3, mapped.LowerBound);
        Assert.True(mapped.IndependenceCertified);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void DependentReplacementDoesNotCertifyAllPoints(bool useSum)
    {
        var curve = Curve();
        var points = Points();
        points[2] = useSum ? curve.Add(points[0], points[1]) : points[0];
        var certificate = curve.GetRankLowerBound(points);
        Assert.Equal(2, certificate.LowerBound);
        Assert.Equal(3, certificate.PointCount);
        Assert.False(certificate.IndependenceCertified);
    }

    [Fact]
    public void TorsionAndInfinityDoNotGivePositiveRank()
    {
        var torsion = new EllipticCurveQ(0, 0, 0, -1, 0);
        var certificate = torsion.GetRankLowerBound(new[]
        {
            new EllipticCurvePoint(0, 0), new(1, 0), new(-1, 0), EllipticCurvePoint.Infinity
        });
        Assert.Equal(0, certificate.LowerBound);
        Assert.Equal(4, certificate.PointCount);
        Assert.False(certificate.IndependenceCertified);
        Assert.Equal(0, Curve().GetRankLowerBound(Array.Empty<EllipticCurvePoint>()).LowerBound);
    }

    [Fact]
    public void InvalidPointsAndCancellationAreRejected()
    {
        var curve = Curve();
        Assert.Throws<ArgumentException>(() => curve.GetRankLowerBound(new[] { new EllipticCurvePoint(0, 3) }));
        Assert.Throws<ArgumentNullException>(() => curve.GetRankLowerBound(null));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).GetRankLowerBound(Points()));
        Assert.Throws<OperationCanceledException>(() => curve.GetRankLowerBound(Points(), cancellationToken: new(true)));
    }
}

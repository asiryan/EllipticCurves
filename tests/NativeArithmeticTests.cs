using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class NativeArithmeticTests
{
    [Theory]
    [InlineData(0, -17, 0, 72, 0, 48)]
    [InlineData(0, -1, 1, -10, -20, 11)]
    [InlineData(0, 0, 1, -1, 0, 37)]
    [InlineData(0, 1, 1, -2, 0, 389)]
    [InlineData(0, 0, 1, 0, 0, 27)]
    [InlineData(0, 0, 0, -1, 0, 32)]
    [InlineData(0, 0, 0, 1, 0, 64)]
    [InlineData(0, 0, 0, 0, 1, 36)]
    [InlineData(0, 0, 0, 0, -1, 144)]
    [InlineData(0, 0, 0, 0, -2, 1728)]
    [InlineData(0, 0, 0, 0, -4, 432)]
    public void ConductorAndMinimalModelAreInvariantUnderCoordinateChanges(
        int a1, int a2, int a3, int a4, int a6, int conductor)
    {
        var e = new EllipticCurveQ(a1, a2, a3, a4, a6);
        Assert.Equal(new BigInteger(conductor), e.Conductor);
        var minimal = e.GlobalMinimalModel;
        Assert.True(e.IsIsomorphic(minimal));
        Assert.Equal(minimal, minimal.GlobalMinimalModel);
        foreach (var u in new BigRational[] { 6, new BigRational(1, 6), -2 })
        {
            var changed = ChangeCoordinates(e, u, 5, -3, 7);
            Assert.Equal(new BigInteger(conductor), changed.Conductor);
            Assert.Equal(minimal, changed.GlobalMinimalModel);
        }
    }

    [Fact]
    public void ReadmeCurveHasExactNativeRankZero()
    {
        var e = new EllipticCurveQ(0, -17, 0, 72, 0);
        Assert.Equal(new EllipticCurveQ(0, 1, 0, -24, 36), e.GlobalMinimalModel);
        var bounds = e.GetRankBounds();
        Assert.True(bounds.UsedTwoIsogenyDescent);
        Assert.True(bounds.IsExact);
        Assert.Equal(0, bounds.ExactRank);
    }

    [Theory]
    [InlineData(1, 0)]
    [InlineData(2, 0)]
    [InlineData(3, 0)]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(34, 2)]
    public void TwoIsogenyDescentProvesKnownRanks(int n, int rank)
    {
        var e = new EllipticCurveQ(0, 0, 0, -n * n, 0);
        var result = e.GetRankBounds();
        Assert.Equal(rank, result.ExactRank);
        Assert.Equal(rank, ChangeCoordinates(e, new BigRational(2, 3), 4, -2, 3).GetRankBounds().ExactRank);
    }

    [Fact]
    public void UnsearchedPointsDoNotProveRankZero()
    {
        var e = new EllipticCurveQ(0, 0, 0, -25, 0);
        var result = e.GetRankBounds(searchBound: 0);
        Assert.Equal(0, result.LowerBound);
        Assert.True(result.UpperBound >= 1);
        Assert.False(result.IsExact);
        Assert.Null(result.ExactRank);
    }

    [Fact]
    public void WithoutRationalTwoTorsionUpperBoundIsExplicitlyUnknown()
    {
        var rankOne = new EllipticCurveQ(0, 0, 1, -1, 0).GetRankBounds();
        Assert.Equal(1, rankOne.LowerBound);
        Assert.Null(rankOne.UpperBound);
        Assert.False(rankOne.UsedTwoIsogenyDescent);
        var rankZero = new EllipticCurveQ(0, -1, 1, -10, -20).GetRankBounds();
        Assert.Equal(0, rankZero.LowerBound);
        Assert.Null(rankZero.ExactRank);
    }

    [Fact]
    public void RejectsSingularCurvesInvalidBoundsAndTruncatedDescent()
    {
        var singular = new EllipticCurveQ(0, 0, 0, 0, 0);
        Assert.Throws<InvalidOperationException>(() => singular.GetConductor());
        Assert.Throws<InvalidOperationException>(() => singular.GetGlobalMinimalModel());
        Assert.Throws<InvalidOperationException>(() => singular.GetRankBounds());
        var e = new EllipticCurveQ(0, 0, 0, -25, 0);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetRankBounds(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetRankBounds(maxSquareClasses: 1));
        Assert.Throws<NotSupportedException>(() => e.GetRankBounds(maxSquareClasses: 2));
    }

    [Fact]
    public void ComputationsHonorCancellation()
    {
        var e = new EllipticCurveQ(0, -17, 0, 72, 0);
        var token = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => e.GetConductor(token));
        Assert.Throws<OperationCanceledException>(() => e.GetGlobalMinimalModel(token));
        Assert.Throws<OperationCanceledException>(() => e.GetRankBounds(cancellationToken: token));
    }

    internal static EllipticCurveQ ChangeCoordinates(EllipticCurveQ e, BigRational u,
        BigRational r, BigRational s, BigRational t) => new(
        (e.A1 + 2 * s) / u,
        (e.A2 - s * e.A1 + 3 * r - s * s) / BigRational.Pow(u, 2),
        (e.A3 + r * e.A1 + 2 * t) / BigRational.Pow(u, 3),
        (e.A4 - s * e.A3 + 2 * r * e.A2 - (t + r * s) * e.A1 + 3 * r * r - 2 * s * t) / BigRational.Pow(u, 4),
        (e.A6 + r * e.A4 + r * r * e.A2 + r * r * r - t * e.A3 - t * t - r * t * e.A1) / BigRational.Pow(u, 6));
}

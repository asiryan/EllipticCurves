using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class RankLowerBoundSearchTests
{
    [Theory]
    [InlineData(0, 0, 0, -1, 0)]
    [InlineData(0, 0, 1, -7, 6)]
    [InlineData(1, 1, 1, -5, 2)]
    [InlineData(-2, 3, -1, -8, 4)]
    public void SieveMatchesIndependentRationalEnumeration(int a1, int a2, int a3, int a4, int a6)
    {
        var curve = new EllipticCurveQ(a1, a2, a3, a4, a6);
        var result = curve.SearchRankLowerBound(new()
        {
            NumeratorRadius = 80, DenominatorRootBound = 4,
            TimeLimit = TimeSpan.FromSeconds(20)
        });
        var expected = curve.RationalPoints(80, 16)
            .Where(p => !p.IsInfinity).Select(p => p.X).ToHashSet();
        Assert.Equal(RankLowerBoundSearchStopReason.SearchBoxExhausted, result.StopReason);
        Assert.True(expected.SetEquals(result.Points.Select(p => p.X)));
        Assert.Equal(expected.Count, result.Points.Count);
        Assert.All(result.Points, p => Assert.True(curve.IsOnCurve(p)));
        Assert.Equal(curve.GetRankLowerBound(result.Points).LowerBound, result.LowerBound);
    }

    [Fact]
    public void SmallCurveCertifiesRankThreeWithoutKnownPoints()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -7, 6);
        var result = curve.SearchRankLowerBound(new()
        { NumeratorRadius = 20, DenominatorRootBound = 3, TargetLowerBound = 3 });
        Assert.Equal(3, result.LowerBound);
        Assert.Equal(RankLowerBoundSearchStopReason.TargetReached, result.StopReason);
    }

    [Fact]
    public void RationalCoefficientNormalizationReturnsOriginalCoordinates()
    {
        var curve = new EllipticCurveQ(0, 0, 0, new BigRational(-1, 16), new BigRational(1, 64));
        var result = curve.SearchRankLowerBound(new() { NumeratorRadius = 0, DenominatorRootBound = 1 });
        Assert.Equal(new BigInteger(64), result.IntegralScale);
        Assert.Contains(new EllipticCurvePoint(0, new BigRational(1, 8)), result.Points);
        Assert.All(result.Points, p => Assert.True(curve.IsOnCurve(p)));
    }

    [Fact]
    public void HugeNumeratorCenterAndBothIntervalEndpointsAreSupported()
    {
        BigInteger x = BigInteger.Pow(10, 50) + 17;
        var curve = new EllipticCurveQ(1, 1, 1, 0, new BigRational(2 + x - x*x*x - x*x));
        foreach (int offset in new[] { -31, 31 })
        {
            var result = curve.SearchRankLowerBound(new()
            { NumeratorCenter = x + offset, NumeratorRadius = 31, DenominatorRootBound = 1 });
            Assert.Contains(result.Points, p => p.X == new BigRational(x));
            Assert.All(result.Points, p => Assert.True(curve.IsOnCurve(p)));
        }
    }

    [Fact]
    public void SieveWordAndBlockBoundariesMatchBruteForce()
    {
        var curve = new EllipticCurveQ(1, 1, 1, -5, 2);
        var result = curve.SearchRankLowerBound(new()
        { NumeratorRadius = 32770, DenominatorRootBound = 1, TimeLimit = TimeSpan.FromSeconds(20) });
        var expected = curve.RationalPoints(32770, 1).Where(p => !p.IsInfinity).Select(p => p.X).ToHashSet();
        Assert.Equal(2, result.SieveBlocks);
        Assert.True(expected.SetEquals(result.Points.Select(p => p.X)));
    }

    [Fact]
    public void TorsionPointsAndAnEmptyBoxDoNotProvePositiveOrExactRank()
    {
        var torsion = new EllipticCurveQ(0, 0, 0, -1, 0);
        var result = torsion.SearchRankLowerBound(new() { NumeratorRadius = 1, DenominatorRootBound = 1 });
        Assert.Equal(0, result.LowerBound);
        Assert.Equal(3, result.Points.Count);
        var positiveRank = new EllipticCurveQ(0, 0, 0, -2, 0);
        var empty = positiveRank.SearchRankLowerBound(new()
        { NumeratorCenter = -10, NumeratorRadius = 0, DenominatorRootBound = 1 });
        Assert.Equal(0, empty.LowerBound);
        Assert.Empty(empty.Points);
        Assert.Equal(RankLowerBoundSearchStopReason.SearchBoxExhausted, empty.StopReason);
    }

    [Fact]
    public void ResourceLimitsAndProgressPreserveCertificates()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -7, 6);
        var squareLimit = curve.SearchRankLowerBound(new() { MaxSquareTests = 0 });
        Assert.Equal(RankLowerBoundSearchStopReason.SquareTestLimit, squareLimit.StopReason);
        Assert.Equal(0, squareLimit.SquareTests);
        var partial = curve.SearchRankLowerBound(new() { NumeratorRadius = 10, MaxSquareTests = 1 });
        Assert.Equal(RankLowerBoundSearchStopReason.SquareTestLimit, partial.StopReason);
        Assert.Equal(1, partial.SquareTests);
        Assert.Equal(1, partial.LowerBound);
        Assert.Single(partial.Points);
        var timeLimit = curve.SearchRankLowerBound(new() { TimeLimit = TimeSpan.FromTicks(1) });
        Assert.Equal(RankLowerBoundSearchStopReason.TimeLimit, timeLimit.StopReason);
        var pointLimit = curve.SearchRankLowerBound(new() { NumeratorRadius = 10, MaxPoints = 1 });
        Assert.Equal(RankLowerBoundSearchStopReason.PointLimit, pointLimit.StopReason);
        Assert.Single(pointLimit.Points);
        var updates = new CaptureProgress();
        var target = curve.SearchRankLowerBound(new() { NumeratorRadius = 10, TargetLowerBound = 2 }, progress: updates);
        Assert.Equal(2, target.LowerBound);
        Assert.NotEmpty(updates.Items);
        Assert.Equal(2, updates.Items.Last().LowerBound);
        Assert.All(updates.Items, r => Assert.Equal(curve.GetRankLowerBound(r.Points).LowerBound, r.LowerBound));
    }

    [Fact]
    public void InvalidInputsAndCancellationAreRejected()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -7, 6);
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.SearchRankLowerBound(new() { NumeratorRadius = -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.SearchRankLowerBound(new() { DenominatorRootBound = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.SearchRankLowerBound(new() { TargetLowerBound = 0 }));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).SearchRankLowerBound());
        Assert.Throws<OperationCanceledException>(() => curve.SearchRankLowerBound(cancellationToken: new(true)));
    }

    private sealed class CaptureProgress : IProgress<RankLowerBoundSearchResult>
    {
        internal List<RankLowerBoundSearchResult> Items { get; } = new();
        public void Report(RankLowerBoundSearchResult value) => Items.Add(value);
    }
}

using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class RankLimitTests
{
    private static EllipticCurveQ RankOneCurve() => new(0, 0, 0, -25, 0);

    private static void AssertIncomplete(RankBounds result, int lower)
    {
        Assert.Equal(lower, result.LowerBound);
        Assert.Null(result.UpperBound);
        Assert.Null(result.TwoSelmerDimension);
        Assert.False(result.IsExact);
        Assert.Null(result.ExactRank);
        Assert.True(result.UsedTwoIsogenyDescent);
        Assert.False(result.UsedGeneralTwoDescent);
        Assert.Contains("MaxSquareClasses", result.Reason);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void FirstImageClassLimitPreservesEarlyPointProof(int maxClasses)
    {
        // The first image has four signed classes, so no descent steps are taken.
        var result = RankOneCurve().GetRankBounds(maxSquareClasses: maxClasses);
        AssertIncomplete(result, 1);
        Assert.Equal(0, result.DescentWork);
        Assert.True(result.PointSearchWork > 0);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(7)]
    public void SecondImageClassLimitPreservesAlreadyProvedLowerBound(int maxClasses)
    {
        // The first image fits four classes; the second requires eight.
        var result = RankOneCurve().GetRankBounds(new RankComputationOptions { MaxSquareClasses = maxClasses });
        AssertIncomplete(result, 1);
        Assert.True(result.DescentWork > 0);
    }

    [Fact]
    public void CompletedFirstImageProofSurvivesWithoutEarlyPointWitnesses()
    {
        var curve = new EllipticCurveQ(0, -38, 0, 74, 0);
        // The minimal-model search box at bound one contains only infinity.
        Assert.Equal(EllipticCurvePoint.Infinity, Assert.Single(curve.GlobalMinimalModel.RationalPoints(1, 1)));
        // The d=2 covering has (u,v,Y)=(1,1,1), giving (2,2) on the input
        // equation. Its first image has eight classes, but the second has sixteen.
        var result = curve.GetRankBounds(new RankComputationOptions { SearchBound = 1, MaxSquareClasses = 8 });
        AssertIncomplete(result, 1);
        Assert.True(result.DescentWork > 0);
    }

    [Fact]
    public void ExactClassThresholdAllowsBothImagesToComplete()
    {
        var curve = RankOneCurve();
        var result = curve.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 8 });
        var baseline = curve.GetRankBounds();
        Assert.Equal(1, result.ExactRank);
        Assert.Equal(baseline.LowerBound, result.LowerBound);
        Assert.Equal(baseline.UpperBound, result.UpperBound);
        Assert.Equal(baseline.DescentWork, result.DescentWork);
        Assert.Equal(baseline.PointSearchWork, result.PointSearchWork);
        Assert.DoesNotContain("MaxSquareClasses", result.Reason);
    }

    [Fact]
    public void GeneralDescentUsesTheSamePartialResultAtItsClassThreshold()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -7, 6);
        var partial = curve.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 7 });
        Assert.Equal(3, partial.LowerBound);
        Assert.Null(partial.UpperBound);
        Assert.Null(partial.TwoSelmerDimension);
        Assert.False(partial.IsExact);
        Assert.Null(partial.ExactRank);
        Assert.True(partial.UsedGeneralTwoDescent);
        Assert.False(partial.UsedTwoIsogenyDescent);
        Assert.Contains("MaxSquareClasses", partial.Reason);

        var complete = curve.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 8 });
        Assert.Equal(3, complete.ExactRank);
        Assert.Equal(3, complete.TwoSelmerDimension);
        Assert.DoesNotContain("MaxSquareClasses", complete.Reason);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ClassLimitWithoutPointProofDoesNotClaimRankZero(bool exhaustPointWork)
    {
        var options = new RankComputationOptions { MaxSquareClasses = 2 };
        if (exhaustPointWork) options.MaxPointSearchWork = 0;
        else options.SearchBound = 0;
        var result = RankOneCurve().GetRankBounds(options);
        AssertIncomplete(result, 0);
        Assert.Equal(0, result.PointSearchWork);
        if (exhaustPointWork) Assert.Contains("MaxPointSearchWork", result.Reason);
    }

    [Fact]
    public void PointWorkLimitKeepsCompletedIsogenyUpperBound()
    {
        var result = RankOneCurve().GetRankBounds(new RankComputationOptions
            { MaxSquareClasses = 8, MaxPointSearchWork = 0 });
        Assert.Equal(0, result.LowerBound);
        Assert.Equal(1, result.UpperBound);
        Assert.Null(result.ExactRank);
        Assert.Equal(0, result.PointSearchWork);
        Assert.Contains("MaxPointSearchWork", result.Reason);
        Assert.DoesNotContain("MaxSquareClasses", result.Reason);
    }

    [Fact]
    public void DescentWorkLimitStillPreservesEarlyPointProof()
    {
        var result = RankOneCurve().GetRankBounds(new RankComputationOptions
            { MaxSquareClasses = 8, MaxDescentWork = 0 });
        Assert.Equal(1, result.LowerBound);
        Assert.Null(result.UpperBound);
        Assert.Null(result.ExactRank);
        Assert.Equal(0, result.DescentWork);
        Assert.Contains("MaxDescentWork", result.Reason);
        Assert.DoesNotContain("MaxSquareClasses", result.Reason);
    }

    [Fact]
    public void InvalidClassLimitsAndCancellationStillThrow()
    {
        var curve = RankOneCurve();
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetRankBounds(maxSquareClasses: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 1 }));
        var token = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => curve.GetRankBounds(maxSquareClasses: 2, cancellationToken: token));
        Assert.Throws<OperationCanceledException>(() => curve.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 2 }, token));
    }
}

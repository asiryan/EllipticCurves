using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class NativeArithmeticTests
{
    [Fact]
    public void ConductorFactorizationReturnsConductorExponentsAndCannotBeModified()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var conductor = curve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = 1 }, out var factors);
        Assert.Equal(new BigInteger(32), conductor);
        var factor = Assert.Single(factors);
        Assert.Equal(new BigInteger(2), factor.Key);
        Assert.Equal(5, factor.Value); // The discriminant is 64 = 2^6, not the conductor.
        Assert.Throws<NotSupportedException>(() => { ((IDictionary<BigInteger, int>)factors)[2] = 99; });
        Assert.Throws<ArgumentNullException>(() => curve.GetConductor(null, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = -1 }, out _));
        Assert.Throws<OperationCanceledException>(() => curve.GetConductor(new FactorizationOptions(), out _, new CancellationToken(true)));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(12)]
    public void ConductorFactorizationPreservesPowersUnderCoordinateChanges(int workers)
    {
        var curve = ChangeCoordinates(new EllipticCurveQ(0, -17, 0, 72, 0), new BigRational(2, 3), 5, -3, 7);
        var conductor = curve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = workers }, out var factors);
        Assert.Equal(new BigInteger(48), conductor);
        Assert.Equal(new BigInteger[] { 2, 3 }, factors.Keys);
        Assert.Equal(new[] { 4, 1 }, factors.Values);
        Assert.Equal(conductor, factors.Aggregate(BigInteger.One, (n, factor) => n * BigInteger.Pow(factor.Key, factor.Value)));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(8)]
    public void ConductorWorkerOptionsPreserveExactResultsAcrossModels(int workers)
    {
        var options = new FactorizationOptions { MaxDegreeOfParallelism = workers };
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.Equal(new BigInteger(32), curve.GetConductor(options, default));
        var changed = ChangeCoordinates(curve, new BigRational(2, 3), 5, -3, 7);
        Assert.Equal(new BigInteger(32), changed.GetConductor(options, default));
        Assert.Equal(new BigInteger(32), curve.GetConductor(default));
        Assert.Equal(workers, options.MaxDegreeOfParallelism);
        Assert.Throws<OperationCanceledException>(() => curve.GetConductor(options, new CancellationToken(true)));
    }

    [Fact]
    public void ConductorRejectsInvalidWorkerOptions()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.Throws<ArgumentNullException>(() => curve.GetConductor((FactorizationOptions)null, default));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = -1 }, default));
    }

    [Fact]
    public void LargeDiscriminantConductorMatchesPari()
    {
        var curve = new EllipticCurveQ(0, 1, 0,
            new BigRational(BigInteger.Parse("-221556180740323405132844117936")),
            new BigRational(BigInteger.Parse("35386140191724122461245294467670188433973860")));
        // PARI/GP ellglobalred with default(factor_proven, 1).
        var expected = BigInteger.Parse("1103561624055499058867562340698878392772504928025988266715523317532246643920");
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        Assert.Equal(expected, curve.GetConductor(timeout.Token));
        var local = curve.GetLocalData(timeout.Token);
        Assert.Equal(new[] { "2", "3", "5", "11", "31", "157", "670606297099",
            "1575838430456954508271967", "81274068710384465721193186106423" },
            local.Select(data => data.Prime.ToString()).ToArray());
        Assert.Equal(new[] { 10, 8, 2, 4, 3, 2, 1, 1, 1 }, local.Select(data => data.DiscriminantValuation).ToArray());
        Assert.Equal(new[] { 4, 1, 1, 1, 1, 1, 1, 1, 1 }, local.Select(data => data.ConductorValuation).ToArray());
        Assert.Equal(expected, local.Aggregate(BigInteger.One,
            (product, data) => product * BigInteger.Pow(data.Prime, data.ConductorValuation)));
    }

    [Fact]
    public void LargeDiscriminantDoesNotPreventMinimalModelReduction()
    {
        var minimal = new EllipticCurveQ(0, 1, 0,
            new BigRational(BigInteger.Parse("-221556180740323405132844117936")),
            new BigRational(BigInteger.Parse("35386140191724122461245294467670188433973860")));
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(5));
        Assert.Equal(minimal, minimal.GetGlobalMinimalModel(timeout.Token));
        // Exercise denominator clearing and scaling at 2, 3, 5 and a larger prime.
        foreach (var scale in new BigRational[] { 6060, new(1, 6060), -6060 })
        {
            var changed = ChangeCoordinates(minimal, scale, 5, -3, 7);
            Assert.Equal(minimal, changed.GetGlobalMinimalModel(timeout.Token));
        }
    }

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
    public void WithoutRationalTwoTorsionGeneralDescentProvesRanks()
    {
        var rankOne = new EllipticCurveQ(0, 0, 1, -1, 0).GetRankBounds();
        Assert.Equal(1, rankOne.LowerBound);
        Assert.Equal(1, rankOne.ExactRank);
        Assert.True(rankOne.UsedGeneralTwoDescent);
        Assert.False(rankOne.UsedTwoIsogenyDescent);
        var rankZero = new EllipticCurveQ(0, -1, 1, -10, -20).GetRankBounds();
        Assert.Equal(0, rankZero.LowerBound);
        Assert.Equal(0, rankZero.ExactRank);
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

using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class PariReferenceTests
{
    public static IEnumerable<object[]> Conductors() => ReadRows("conductors.csv");
    public static IEnumerable<object[]> Ranks() => ReadRows("ranks.csv");

    private static IEnumerable<object[]> ReadRows(string name)
    {
        foreach (var line in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", name)).Distinct())
            if (!line.StartsWith("#") && !string.IsNullOrWhiteSpace(line))
                yield return new object[] { line };
    }

    [Theory]
    [MemberData(nameof(Conductors))]
    public void MatchesIndependentlyComputedPariConductorAndMinimalModel(string row)
    {
        var a = row.Split(',').Select(BigInteger.Parse).ToArray();
        var e = new EllipticCurveQ(new(a[0]), new(a[1]), new(a[2]), new(a[3]), new(a[4]));
        Assert.Equal(a[5], e.Conductor);
        Assert.Equal(new EllipticCurveQ(new(a[6]), new(a[7]), new(a[8]), new(a[9]), new(a[10])), e.GlobalMinimalModel);
    }

    [Theory]
    [MemberData(nameof(Ranks))]
    public void BoundsContainIndependentlyProvedPariRank(string row)
    {
        var a = row.Split(',').Select(int.Parse).ToArray();
        var result = new EllipticCurveQ(0, a[0], 0, a[1], 0).GetRankBounds();
        Assert.True(result.LowerBound <= a[2], $"Computed {result}, actual rank {a[2]}");
        Assert.True(result.UpperBound >= a[2], $"Computed {result}, actual rank {a[2]}");
    }
}

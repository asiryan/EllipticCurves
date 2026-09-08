using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class RootNumberTests
{
    [Fact]
    public void LocalSignsCoverWeightedShortModelsAndZeroInvariants()
    {
        foreach (var row in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "local-roots.csv")))
        {
            var a = row.Split(',').Select(BigInteger.Parse).ToArray();
            var e = new EllipticCurveQ(new(a[0]), new(a[1]), new(a[2]), new(a[3]), new(a[4]));
            Assert.True((int)a[5] == EllipticCurveQ.LocalRootNumber(e, 2), $"Root at 2: {row}");
            Assert.True((int)a[6] == EllipticCurveQ.LocalRootNumber(e, 3), $"Root at 3: {row}");
        }
    }

    public static IEnumerable<object[]> Rows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "root-numbers.csv"))
        .Select(row => new object[] { row });

    [Theory]
    [MemberData(nameof(Rows))]
    public void MatchesIndependentPariLocalAndGlobalSigns(string row)
    {
        var a = row.Split(',').Select(BigInteger.Parse).ToArray();
        var e = new EllipticCurveQ(new(a[0]), new(a[1]), new(a[2]), new(a[3]), new(a[4]));
        Assert.Equal((int)a[5], e.RootNumber);
        var m = e.GlobalMinimalModel;
        Assert.Equal((int)a[6], EllipticCurveQ.LocalRootNumber(m, 2));
        Assert.Equal((int)a[7], EllipticCurveQ.LocalRootNumber(m, 3));
    }
}

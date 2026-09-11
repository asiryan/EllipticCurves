using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class TorsionRegressionTests
{
    public static IEnumerable<object[]> ReferenceCurves() =>
        File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "torsion.csv"))
            .Select(line => new object[] { line });

    private static BigRational Rational(string text)
    {
        var parts = text.Split('/');
        return new BigRational(BigInteger.Parse(parts[0]), parts.Length == 1 ? BigInteger.One : BigInteger.Parse(parts[1]));
    }

    [Theory]
    [MemberData(nameof(ReferenceCurves))]
    public void FullTorsionSetsMatchPariOnOriginalAndChangedModels(string row)
    {
        var fields = row.Split(',');
        var a = fields.Take(5).Select(Rational).ToArray();
        var curve = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        var expected = fields[6].Split(';').Select(text =>
        {
            if (text == "O") return EllipticCurvePoint.Infinity;
            var xy = text.Split(':');
            return new EllipticCurvePoint(Rational(xy[0]), Rational(xy[1]));
        }).ToHashSet();
        var actual = curve.TorsionPoints.ToArray();
        Assert.Equal(expected.Count, actual.Length);
        Assert.True(expected.SetEquals(actual), row);
        Assert.Equal(fields[5], curve.TorsionStructure);
        Assert.All(actual, point =>
        {
            Assert.True(curve.IsOnCurve(point));
            var order = curve.TorsionOrder(point);
            Assert.NotNull(order);
            Assert.True(curve.Multiply(point, order.Value).IsInfinity);
        });
    }

    [Fact]
    public void IntegerCubicRootsMatchExhaustiveSearchIncludingRepeatedRoots()
    {
        for (int a = -30; a <= 30; a++)
            for (int b = -30; b <= 30; b++)
            {
                // All roots lie in [-31,31] by Cauchy's bound.
                var expected = Enumerable.Range(-31, 63).Where(x => x * x * x + a * x + b == 0)
                    .Select(x => new BigInteger(x)).ToArray();
                Assert.Equal(expected, InternalMath.IntegralShortCubicRoots(a, b).ToArray());
            }
    }

    [Fact]
    public void IntegerCubicRootSearchHandlesHugeAndCloselySpacedRoots()
    {
        var n = BigInteger.Pow(10, 100);
        // Roots n, n+1 and -2n-1 straddle a stationary point near n+1/2.
        var a = -3 * n * n - 3 * n - 1;
        var b = n * (n + 1) * (2 * n + 1);
        Assert.Equal(new[] { -2 * n - 1, n, n + 1 }, InternalMath.IntegralShortCubicRoots(a, b).ToArray());
        Assert.Equal(new[] { n }, InternalMath.IntegralShortCubicRoots(1, -n * n * n - n).ToArray());
    }
}

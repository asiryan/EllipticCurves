using System.Numerics;
using System.Reflection;
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
        // Exercise single-point queries before preparing the full subgroup.
        Assert.All(expected, point =>
        {
            Assert.True(curve.IsOnCurve(point));
            var order = curve.TorsionOrder(point);
            Assert.NotNull(order);
            Assert.True(curve.Multiply(point, order.Value).IsInfinity);
            for (int n = 1; n < order.Value; n++) Assert.False(curve.Multiply(point, n).IsInfinity);
            Assert.True(curve.IsTorsionPoint(point));
        });
        AssertTorsionCacheUnprepared(curve);
        var actual = curve.TorsionPoints.ToArray();
        Assert.Equal(expected.Count, actual.Length);
        Assert.True(expected.SetEquals(actual), row);
        Assert.Equal(fields[5], curve.TorsionStructure);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void SinglePointQueriesAvoidCostlyColdSubgroupPreparation(int order)
    {
        // Factoring the discriminant would include this large composite. The
        // regression checks cache independence directly, without a time limit.
        var large = BigInteger.Parse("21785392458764315483988614758901932764423833726404496768609056369259382575196330");
        var curve = order == 2 ? new EllipticCurveQ(0, 0, 0, new BigRational(large), 0)
            : new EllipticCurveQ(0, 0, 0, 0, new BigRational(large * large));
        var point = new EllipticCurvePoint(0, order == 2 ? BigRational.Zero : new BigRational(large));
        Assert.Equal(order, curve.TorsionOrder(point));
        Assert.True(curve.IsTorsionPoint(point));
        AssertTorsionCacheUnprepared(curve);
    }

    [Fact]
    public void SinglePointQueriesRejectLargeInfiniteOrderCoordinates()
    {
        var curve = new EllipticCurveQ(0, 0, 0, 0, -2);
        // (3,5) has infinite order by Lutz-Nagell: 25 does not divide 1728.
        var point = curve.Multiply(new(3, 5), 32);
        Assert.True(point.X.Num.ToByteArray().Length > 100);
        Assert.Null(curve.TorsionOrder(point));
        Assert.False(curve.IsTorsionPoint(point));
        AssertTorsionCacheUnprepared(curve);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ExactPointOrderStillWorksWhenSmallReductionsAreUnavailable(bool rationalModel)
    {
        var product = new BigInteger(5 * 7 * 11 * 13);
        var y = rationalModel ? new BigRational(1, product) : new BigRational(product);
        var curve = new EllipticCurveQ(0, 0, 0, 0, y * y);
        Assert.Equal(3, curve.TorsionOrder(new(0, y)));
        Assert.True(curve.IsTorsionPoint(new(0, y)));
        AssertTorsionCacheUnprepared(curve);
    }

    [Fact]
    public void SinglePointQueriesPreserveValidationAndSingularInfinityBehavior()
    {
        var singular = new EllipticCurveQ(0, 0, 0, 0, 0);
        Assert.Equal(1, singular.TorsionOrder(EllipticCurvePoint.Infinity));
        Assert.Throws<InvalidOperationException>(() => singular.TorsionOrder(new(0, 0)));
        Assert.Throws<InvalidOperationException>(() => singular.IsTorsionPoint(new(0, 0)));
        Assert.Throws<InvalidOperationException>(() => singular.IsTorsionPoint(EllipticCurvePoint.Infinity));
        Assert.Throws<ArgumentException>(() => singular.TorsionOrder(new(1, 0)));
        Assert.Throws<ArgumentException>(() => singular.IsTorsionPoint(new(1, 0)));
    }

    private static void AssertTorsionCacheUnprepared(EllipticCurveQ curve)
    {
        var cache = typeof(EllipticCurveQ).GetField("torsionPointsCache", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(curve)!;
        Assert.Null(cache.GetType().GetField("value", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(cache));
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
    public void NearbyIntegerRootsAllContributeTwoTorsionPoints()
    {
        // x^3 - 7x - 6 = (x+2)(x+1)(x-3). Splitting at ceil(sqrt(2))
        // incorrectly puts both negative roots on the same search interval.
        var curve = new EllipticCurveQ(0, 0, 0, -7, -6);
        var expected = new HashSet<EllipticCurvePoint>
        {
            EllipticCurvePoint.Infinity, new(-2, 0), new(-1, 0), new(3, 0)
        };
        Assert.True(expected.SetEquals(curve.TorsionPoints));
        Assert.Equal("Z/2Z x Z/2Z", curve.TorsionStructure);
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

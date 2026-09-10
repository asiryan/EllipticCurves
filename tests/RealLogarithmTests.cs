using EllipticCurves;
using Xunit;
using static EllipticCurves.Tests.ExtendedReferenceTests;

namespace EllipticCurves.Tests;

public class RealLogarithmTests
{
    public static IEnumerable<object[]> LogRows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "real-logarithms.csv")).Select(row => new object[] { row });

    [Theory, MemberData(nameof(LogRows))]
    public void LogarithmsMatchIndependentPariValues(string row)
    {
        var v = row.Split(','); var e = new EllipticCurveQ(Parse(v[0]), Parse(v[1]), Parse(v[2]), Parse(v[3]), Parse(v[4]));
        var p = new EllipticCurvePoint(Parse(v[5]), Parse(v[6])); var result = e.RealEllipticLogarithm(p);
        double reference = RealEnclosure.ToDouble(Parse(v[7])); double imaginary = RealEnclosure.ToDouble(Parse(v[8]));
        double period = RealEnclosure.ToDouble(Parse(v[9]));
        Assert.True(CircleDistance(reference, result.RealPart, period) < 1e-10 * period, $"{row}: got {result.RealPart:G17}");
        Assert.True(Math.Abs(imaginary - result.ImaginaryPart) < 1e-10 * period);
        Assert.True(Math.Abs(period - result.PrimitiveRealPeriod) < 1e-10 * period);
        Assert.Equal(imaginary == 0 ? 0 : 1, result.ComponentIndex);
    }

    [Fact]
    public void LogarithmsRespectAdditionNegationAndMinimalModelChanges()
    {
        var e = new EllipticCurveQ(0, 0, 0, -25, 0); var g = new EllipticCurvePoint(new BigRational(25, 4), new BigRational(75, 8));
        var change = e.ChangeModel(new BigRational(-2, 3), 7, 2, -5);
        var points = e.TorsionPoints.Concat(new[] { g, e.Negate(g), e.Double(g), e.Add(g, new EllipticCurvePoint(0, 0)) }).ToArray();
        var logs = points.Select(p => e.RealEllipticLogarithm(p)).ToArray();
        for (int i = 0; i < points.Length; i++)
        {
            var changed = change.Target.RealEllipticLogarithm(change.Map(points[i]));
            // Minimal-model maps can differ by [-1]. Compare with the actual map used by each call.
            var nativeMap = e.GetMinimalModelIsomorphism(); var changedMap = change.Target.GetMinimalModelIsomorphism();
            var a = nativeMap.Map(points[i]); var b = changedMap.Map(change.Map(points[i]));
            double expected = a.Equals(b) ? logs[i].RealPart : logs[i].PrimitiveRealPeriod - logs[i].RealPart;
            Assert.True(CircleDistance(expected, changed.RealPart, logs[i].PrimitiveRealPeriod) < 1e-10);
            Assert.Equal(logs[i].ComponentIndex, changed.ComponentIndex);
            for (int j = 0; j < points.Length; j++)
            {
                var sum = e.RealEllipticLogarithm(e.Add(points[i], points[j]));
                Assert.True(CircleDistance(sum.RealPart, logs[i].RealPart + logs[j].RealPart, sum.PrimitiveRealPeriod) < 1e-10);
                Assert.Equal(logs[i].ComponentIndex ^ logs[j].ComponentIndex, sum.ComponentIndex);
            }
        }
    }

    [Fact]
    public void CarlsonIntegralPreservesHomogeneityAndNearSingularValues()
    {
        foreach (double scale in new[] { 1e-300, 1.0, 1e300 })
            Assert.True(Math.Abs(CarlsonIntegral.Rf(0, scale, scale, 64, default) * Math.Sqrt(scale) - Math.PI / 2) < 1e-14);
        Assert.True(Math.Abs(CarlsonIntegral.Rf(0, 1e-280, 1, 64, default) - (Math.Log(4) + 140 * Math.Log(10))) < 1e-12);
    }

    [Fact]
    public void InfinityAndFailureConditionsAreExplicit()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        Assert.Equal(0, e.RealEllipticLogarithm(EllipticCurvePoint.Infinity).RealPart);
        Assert.Throws<ArgumentException>(() => e.RealEllipticLogarithm(new EllipticCurvePoint(0, 1)));
        Assert.Throws<OperationCanceledException>(() => e.RealEllipticLogarithm(EllipticCurvePoint.Infinity, cancellationToken: new CancellationToken(true)));
        Assert.Throws<ArithmeticException>(() => e.RealEllipticLogarithm(new EllipticCurvePoint(0, 0), new RealComputationOptions { MaxRootWork = 0 }));
        Assert.Throws<ArithmeticException>(() => CarlsonIntegral.Rf(0, 1, 2, 0, default));
        Assert.Throws<ArithmeticException>(() => CarlsonIntegral.Rf(0, 0, 1, 64, default));
        Assert.Equal(1.0, CarlsonIntegral.Rf(1, 1, 1, 64, default));
    }
    private static double CircleDistance(double a, double b, double period)
    { double difference = (a - b) / period; return Math.Abs(difference - Math.Round(difference)) * period; }
}

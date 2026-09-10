using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class ExtensionCurveTests
{
    public static IEnumerable<object[]> Rows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "extension-curves.csv")).Select(s => new object[] { s });

    [Theory, MemberData(nameof(Rows))]
    public void InvariantsCountingAndPointMapsAgreeWithPari(string row)
    {
        var v = row.Split(','); var field = FiniteFieldTests.Field(int.Parse(v[0]), v[1].Split(':').Select(int.Parse).ToArray());
        FiniteFieldElement Element(int i) => FiniteFieldTests.Decode(field, int.Parse(v[i]));
        var e = new EllipticCurveFq(field, Element(2), Element(3), Element(4), Element(5), Element(6));
        Assert.Equal(Element(8), e.Discriminant); Assert.Equal(Element(9), e.C4);
        Assert.Equal(Element(10), e.C6); Assert.Equal(Element(11), e.JInvariant);
        EllipticCurvePointFq Point(int i) => v[i] == "-1" ? EllipticCurvePointFq.Infinity : e.CreatePoint(Element(i), Element(i + 1));
        var p = Point(12); var q = Point(14);
        Assert.Equal(Point(16), e.Add(p, q)); Assert.Equal(Point(19), e.Multiply(p, int.Parse(v[18])));
        var points = e.Points().ToArray();
        Assert.Equal(int.Parse(v[7]), points.Length); Assert.Equal(points.Length, points.Distinct().Count());
        Assert.Equal(new BigInteger(points.Length), e.CountPoints());
        Assert.Contains(p, points); Assert.Contains(q, points);
        foreach (var point in points.Take(15))
        {
            Assert.True(e.IsOnCurve(point)); Assert.True(e.Multiply(point, points.Length).IsInfinity);
            Assert.True(e.Add(point, e.Negate(point)).IsInfinity);
            Assert.Equal(e.Add(e.Add(point, p), q), e.Add(point, e.Add(p, q)));
            Assert.Equal(point, e.Subtract(e.Add(point, p), p));
            Assert.Equal(e.Negate(e.Double(point)), e.Multiply(point, -2));
        }
    }

    [Theory, MemberData(nameof(FiniteFieldTests.Fields), MemberType = typeof(FiniteFieldTests))]
    public void BaseChangeCountsMatchFrobeniusRecurrence(int prime, int[] polynomial)
    {
        var field = FiniteFieldTests.Field(prime, polynomial);
        var baseCurve = new EllipticCurveFp(prime, 0, 0, 1, -1, 0);
        var extended = new EllipticCurveFq(field, 0, 0, 1, -1, 0);
        var ap = prime + 1 - baseCurve.CountPoints();
        BigInteger previous = 2, trace = ap;
        for (int i = 2; i <= field.Degree; i++) { var next = ap * trace - prime * previous; previous = trace; trace = next; }
        Assert.Equal(field.Order + 1 - trace, extended.CountPoints());
        Assert.Equal(field.CreateElement(baseCurve.JInvariant), extended.JInvariant);
    }

    [Theory]
    [InlineData(2), InlineData(3), InlineData(5), InlineData(7)]
    public void DegreeOneAgreesWithExistingPrimeFieldApi(int prime)
    {
        var field = FiniteFieldTests.Field(prime, new[] { 0, 1 });
        var e = new EllipticCurveFq(field, 0, 0, 1, -1, 0);
        var baseline = new EllipticCurveFp(prime, 0, 0, 1, -1, 0);
        EllipticCurvePointFq Lift(EllipticCurvePointFp p) => p.IsInfinity ? EllipticCurvePointFq.Infinity : e.CreatePoint(p.X, p.Y);
        Assert.Equal(baseline.CountPoints(), e.CountPoints());
        foreach (var p in baseline.Points()) foreach (var q in baseline.Points())
            Assert.Equal(Lift(baseline.Add(p, q)), e.Add(Lift(p), Lift(q)));
    }

    [Fact]
    public void InvalidPointsFieldsSingularCurvesAndLimitsAreRejected()
    {
        var field = FiniteFieldTests.Field(3, new[] { 1, 0, 1 });
        var other = FiniteFieldTests.Field(3, new[] { 2, 1, 1 });
        var e = new EllipticCurveFq(field, 0, 0, 1, -1, 0); var f = new EllipticCurveFq(other, 0, 0, 1, -1, 0);
        Assert.Throws<ArgumentNullException>(() => new EllipticCurveFq(null, 0, 0, 1, -1, 0));
        Assert.Throws<ArgumentException>(() => new EllipticCurveFq(field, 0, 0, 0, 0, 0));
        Assert.Throws<ArgumentException>(() => e.CreatePoint(0, 1)); Assert.False(e.IsOnCurve(default));
        Assert.Throws<ArgumentException>(() => e.Add(e.CreatePoint(0, 0), f.CreatePoint(0, 0)));
        Assert.Throws<ArgumentException>(() => e.CreatePoint(field.Zero, other.Zero));
        Assert.Throws<ArithmeticException>(() => e.CountPoints(80));
        Assert.Throws<ArithmeticException>(() => e.Points(80).First());
        Assert.True(e.CountPoints(81) > 0);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.CountPoints(-1));
        var cancelled = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => e.CountPoints(cancellationToken: cancelled));
        Assert.Throws<OperationCanceledException>(() => e.Multiply(EllipticCurvePointFq.Infinity, 0, cancelled));
        Assert.Equal(e, new EllipticCurveFq(field, 3, 0, 4, 2, 0));
        Assert.Equal(e.GetHashCode(), new EllipticCurveFq(field, 3, 0, 4, 2, 0).GetHashCode());
        Assert.NotEqual(e, f);
        var large = new EllipticCurveFq(new FiniteField(65537, new BigInteger[] { 0, 1 }), 0, 0, 1, -1, 0);
        Assert.Throws<ArithmeticException>(() => large.CountPoints());
    }
}

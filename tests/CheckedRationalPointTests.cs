using System.Numerics;
using EllipticCurves.Explorer.Computations;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class CheckedRationalPointTests
{
    private static EllipticCurveQ Curve() => new(0, 0, 1, -1, 0);

    [Fact]
    public void FactoriesCheckTheEquationAndCanonicalRationalCoordinates()
    {
        var curve = Curve();
        var expected = new EllipticCurvePoint(new BigRational(1, 4), new BigRational(-5, 8));
        Assert.Equal(expected, curve.CreatePoint(new BigRational(2, 8), new BigRational(-10, 16)));
        Assert.True(curve.TryCreatePoint(expected.X, expected.Y, out var point));
        Assert.Equal(expected, point);
        Assert.Throws<ArgumentException>(() => curve.CreatePoint(0, 1));
        Assert.False(curve.TryCreatePoint(0, 1, out var rejected));
        Assert.True(rejected.IsInfinity);
    }

    [Fact]
    public void SingularModelsRejectCheckedPointsAndEveryCheckedGroupOperation()
    {
        var curve = new EllipticCurveQ(0, 0, 0, 0, 0);
        var infinity = EllipticCurvePoint.Infinity;
        Assert.Throws<InvalidOperationException>(() => curve.CreatePoint(0, 0));
        Assert.False(curve.TryCreatePoint(0, 0, out var rejected));
        Assert.True(rejected.IsInfinity);
        Assert.Throws<InvalidOperationException>(() => curve.AddChecked(infinity, infinity));
        Assert.Throws<InvalidOperationException>(() => curve.DoubleChecked(infinity));
        Assert.Throws<InvalidOperationException>(() => curve.NegateChecked(infinity));
        Assert.Throws<InvalidOperationException>(() => curve.SubtractChecked(infinity, infinity));
        Assert.Throws<InvalidOperationException>(() => curve.MultiplyChecked(infinity, 0));
        // Representing a singular equation and testing its equation remain useful.
        Assert.True(curve.IsOnCurve(new EllipticCurvePoint(0, 0)));
    }

    [Fact]
    public void IdentityAndZeroScalarShortcutsCannotBypassMembershipChecks()
    {
        var curve = Curve();
        var invalid = new EllipticCurvePoint(0, 1);
        var infinity = EllipticCurvePoint.Infinity;
        Assert.Equal("p", Assert.Throws<ArgumentException>(() => curve.AddChecked(invalid, infinity)).ParamName);
        Assert.Equal("q", Assert.Throws<ArgumentException>(() => curve.AddChecked(infinity, invalid)).ParamName);
        Assert.Throws<ArgumentException>(() => curve.DoubleChecked(invalid));
        Assert.Throws<ArgumentException>(() => curve.NegateChecked(invalid));
        Assert.Throws<ArgumentException>(() => curve.SubtractChecked(invalid, infinity));
        Assert.Throws<ArgumentException>(() => curve.SubtractChecked(infinity, invalid));
        foreach (int n in new[] { 0, 1, -1, 2 })
            Assert.Throws<ArgumentException>(() => curve.MultiplyChecked(invalid, n));
    }

    [Fact]
    public void DefaultPointIsCheckedByItsCoordinatesAndNotRejectedUnconditionally()
    {
        var curve = Curve();
        var point = default(EllipticCurvePoint);
        Assert.False(point.IsInfinity);
        Assert.Equal(new EllipticCurvePoint(0, 0), point);
        Assert.Equal(curve.Double(point), curve.DoubleChecked(point));
        var other = new EllipticCurveQ(0, 0, 0, 0, 1);
        Assert.Throws<ArgumentException>(() => other.DoubleChecked(point));
    }

    [Theory]
    [InlineData(-7)]
    [InlineData(-1)]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(7)]
    public void CheckedOperationsKeepExactGroupLawResultsOnGeneralRationalModels(int n)
    {
        var original = Curve();
        var change = original.ChangeModel(new BigRational(2, 3), 5, -2, 7);
        var curve = change.Target;
        var p = change.Map(original.CreatePoint(0, 0));
        var q = change.Map(original.Multiply(new EllipticCurvePoint(0, 0), 3));
        foreach (var point in new[] { p, q, EllipticCurvePoint.Infinity })
        {
            Assert.Equal(curve.Negate(point), curve.NegateChecked(point));
            Assert.Equal(curve.Double(point), curve.DoubleChecked(point));
            Assert.Equal(curve.Multiply(point, n), curve.MultiplyChecked(point, n));
            Assert.Equal(change.Map(original.Multiply(change.MapBack(point), n)), curve.MultiplyChecked(point, n));
            foreach (var other in new[] { p, q, EllipticCurvePoint.Infinity })
            {
                var sum = curve.AddChecked(point, other);
                Assert.Equal(curve.Add(point, other), sum);
                Assert.True(curve.IsOnCurve(sum));
                Assert.Equal(curve.Subtract(point, other), curve.SubtractChecked(point, other));
            }
        }
    }

    [Fact]
    public void CheckedScalarMultiplicationHonorsCancellationBeforeValidationAndTrivialResults()
    {
        var token = new CancellationToken(true);
        foreach (var curve in new[] { Curve(), new EllipticCurveQ(0, 0, 0, 0, 0) })
        foreach (var point in new[] { new EllipticCurvePoint(0, 1), EllipticCurvePoint.Infinity })
        {
            var error = Assert.Throws<OperationCanceledException>(() => curve.MultiplyChecked(point, BigInteger.Zero, token));
            Assert.Equal(token, error.CancellationToken);
        }
    }

    [Fact]
    public void ExplorerPreservesItsExistingCheckedInputActionsAndExposesPointFactories()
    {
        Assert.DoesNotContain(CalculationCatalog.All, operation =>
            operation.Member?.DeclaringType == typeof(EllipticCurveQ) && operation.Member.Name.EndsWith("Checked", StringComparison.Ordinal));
        Assert.Contains(CalculationCatalog.All, operation => operation.Id == "EllipticCurveQ.Add(EllipticCurvePoint,EllipticCurvePoint)");
        Assert.Contains(CalculationCatalog.All, operation => operation.Id == "EllipticCurveQ.CreatePoint(BigRational,BigRational)");
        Assert.Contains(CalculationCatalog.All, operation => operation.Id == "EllipticCurveQ.TryCreatePoint(BigRational,BigRational,EllipticCurvePoint&)");
    }
}

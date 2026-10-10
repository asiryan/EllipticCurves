using System.Numerics;
using EllipticCurves.Explorer.Computations;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class CancellationRegressionTests
{
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void TorsionOperationsHonorPreCancellationAndAllowRetry(bool warmCache)
    {
        var curve = new EllipticCurveQ(0, 0, 0, 0, 1);
        var point = new EllipticCurvePoint(0, 1);
        if (warmCache) Assert.Equal(6, curve.GetTorsionPoints().Count);
        var token = new CancellationToken(true);

        Cancelled(token, () => curve.GetTorsionPoints(token));
        Cancelled(token, () => curve.GetTorsionStructure(token));
        Cancelled(token, () => curve.IsTorsionPoint(point, token));
        Cancelled(token, () => curve.TorsionOrder(point, token));
        Cancelled(token, () => curve.IsTorsionPoint(EllipticCurvePoint.Infinity, token));
        Cancelled(token, () => curve.TorsionOrder(EllipticCurvePoint.Infinity, token));

        // Cancellation must not publish an empty/partial set or poison the cache.
        var torsion = curve.GetTorsionPoints();
        Assert.Equal(6, torsion.Count);
        Assert.Contains(point, torsion);
        Assert.All(torsion, p => Assert.True(curve.IsOnCurve(p)));
        Assert.Equal("Z/6Z", curve.GetTorsionStructure());
        Assert.True(curve.IsTorsionPoint(point, default));
        Assert.Equal(3, curve.TorsionOrder(point, default));
        Assert.Same(torsion, curve.GetTorsionPoints());
        Assert.Equal(torsion, curve.TorsionPoints);
        Assert.Equal(curve.GetTorsionStructure(), curve.TorsionStructure);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(-1)]
    public void ScalarMultiplicationHonorsPreCancellationEvenForTrivialResults(int scalar)
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        var token = new CancellationToken(true);
        Cancelled(token, () => curve.Multiply(new EllipticCurvePoint(0, 0), scalar, token));
        Cancelled(token, () => curve.Multiply(EllipticCurvePoint.Infinity, scalar, token));
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(4)]
    public void ExactSquareTestingHonorsPreCancellationBeforeItsFastPaths(int value)
    {
        var token = new CancellationToken(true);
        Cancelled(token, () => BigRational.IsSquare(value, out _, token));
    }

    [Fact]
    public void TorsionArithmeticHelpersHonorPreCancellationIncludingTrivialInputs()
    {
        var token = new CancellationToken(true);
        var curve = new EllipticCurveQ(0, 0, 0, 0, 1);
        Cancelled(token, () => InternalMath.FactorAbs(0, token));
        Cancelled(token, () => InternalMath.FactorAbs(1, token));
        Cancelled(token, () => InternalMath.FactorAbs(36, token));
        Cancelled(token, () => InternalMath.IntegerSqrt(0, token));
        Cancelled(token, () => InternalMath.IntegerSqrt(1, token));
        Cancelled(token, () => InternalMath.IntegerSqrt(36, token));
        Cancelled(token, () => InternalMath.IntegerKthRoot(0, 3, token));
        Cancelled(token, () => InternalMath.TryIntegerKthRoot(1, 3, out _, token));
        Cancelled(token, () => InternalMath.TryRationalKthRoot(-1, 2, out _, token));
        Cancelled(token, () => InternalMath.ModPow(1, 0, 5, token));
        Cancelled(token, () => InternalMath.Legendre(0, 5, token));
        Cancelled(token, () => InternalMath.CountPointsFpShort(0, 1, 5, token));
        Cancelled(token, () => InternalMath.OrderDivides(curve, EllipticCurvePoint.Infinity, 0, token));
        Cancelled(token, () => InternalMath.IsTorsionWithCandidates(curve, EllipticCurvePoint.Infinity, new(), token));
        Cancelled(token, () => InternalMath.IntegralShortCubicRoots(0, 0, token).ToArray());
        Cancelled(token, () => InternalMath.EnumerateSquareDivisors(new(), token).ToArray());
        Cancelled(token, () => InternalMath.EnumerateDivisorsAbs(0, token).ToArray());
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void DivisorEnumerationCanBeCancelledBetweenResultsAndStartedAgain(bool squareOnly)
    {
        IEnumerable<BigInteger> Divisors(CancellationToken token) => squareOnly
            ? InternalMath.EnumerateSquareDivisors(new() { [2] = 2, [3] = 2 }, token)
            : InternalMath.EnumerateDivisorsAbs(36, token);
        using var cancellation = new CancellationTokenSource();
        using var divisors = Divisors(cancellation.Token).GetEnumerator();
        Assert.True(divisors.MoveNext());
        Assert.Equal(BigInteger.One, divisors.Current);
        cancellation.Cancel();
        Cancelled(cancellation.Token, () => divisors.MoveNext());

        var expected = squareOnly ? new[] { 1, 4, 9, 36 } : new[] { 1, 2, 3, 4, 6, 9, 12, 18, 36 };
        Assert.Equal(expected.Select(n => new BigInteger(n)), Divisors(default).OrderBy(n => n));
    }

    [Fact]
    public void ZeroDivisorConventionStillHonorsCancellationAtEnumerationCompletion()
    {
        using var cancellation = new CancellationTokenSource();
        using var divisors = InternalMath.EnumerateDivisorsAbs(0, cancellation.Token).GetEnumerator();
        Assert.True(divisors.MoveNext());
        Assert.Equal(BigInteger.Zero, divisors.Current);
        cancellation.Cancel();
        Cancelled(cancellation.Token, () => divisors.MoveNext());
        Assert.Equal(new[] { BigInteger.Zero }, InternalMath.EnumerateDivisorsAbs(0).ToArray());
    }

    [Fact]
    public void IntervalHelpersHonorPreCancellationOnExactAndTrivialInputs()
    {
        var token = new CancellationToken(true);
        Cancelled(token, () => DyadicInterval.SqrtInteger(0, token));
        Cancelled(token, () => DyadicInterval.ExpNegative(DyadicInterval.Integer(0), token));
        Cancelled(token, () => DyadicInterval.Pi(token));
        Cancelled(token, () => RationalInterval.Exact(0).Sqrt(token));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void PointSearchHonorsPreCancellation(bool integral)
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        var token = new CancellationToken(true);
        Cancelled(token, () => Search(curve, integral, token).ToArray());
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void PointSearchCanBeCancelledBetweenAffineResultsAndStartedAgain(bool integral)
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        using var cancellation = new CancellationTokenSource();
        using var points = Search(curve, integral, cancellation.Token).GetEnumerator();
        Assert.True(points.MoveNext());
        Assert.True(points.Current.IsInfinity);
        Assert.True(points.MoveNext());
        Assert.False(points.Current.IsInfinity);
        Assert.True(curve.IsOnCurve(points.Current));

        cancellation.Cancel();
        Cancelled(cancellation.Token, () => points.MoveNext());

        var retry = Search(curve, integral, default).ToArray();
        Assert.Contains(new EllipticCurvePoint(0, 0), retry);
        Assert.All(retry, p => Assert.True(curve.IsOnCurve(p)));
    }

    [Fact]
    public void CancellableArithmeticAndSearchKeepTheirExactResults()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        var point = new EllipticCurvePoint(0, 0);
        Assert.Equal(new EllipticCurvePoint(1, 0), curve.Multiply(point, 2, default));
        Assert.Equal(new EllipticCurvePoint(0, -1), curve.Multiply(point, -1, default));
        Assert.Equal(EllipticCurvePoint.Infinity, curve.Multiply(point, BigInteger.Zero, default));

        // x=-1, 0, 1 each has y=0 and y=-1, plus the point at infinity.
        var rational = curve.RationalPoints(1, 1, default).ToArray();
        Assert.Equal(7, rational.Length);
        Assert.Equal(rational, curve.IntegralPoints(1, default).ToArray());
        Assert.Equal(rational, curve.RationalPoints(1, 1).ToArray());
        Assert.Equal(rational, curve.IntegralPoints(1).ToArray());
    }

    [Theory]
    [InlineData("Multiply", "EllipticCurveQ.Multiply(EllipticCurvePoint,BigInteger)")]
    [InlineData("RationalPoints", "EllipticCurveQ.RationalPoints(Int32,Int32)")]
    [InlineData("IntegralPoints", "EllipticCurveQ.IntegralPoints(Int32)")]
    [InlineData("IsTorsionPoint", "EllipticCurveQ.IsTorsionPoint(EllipticCurvePoint)")]
    [InlineData("TorsionOrder", "EllipticCurveQ.TorsionOrder(EllipticCurvePoint)")]
    public void ExplorerKeepsOneLegacyActionAndItsSavedId(string name, string id)
    {
        var operation = Assert.Single(CalculationCatalog.All, o =>
            o.Context == CalculationContext.RationalCurve && o.Member?.Name == name);
        Assert.Equal(id, operation.Id);
        Assert.DoesNotContain(CalculationCatalog.All, o => o.Member?.Name is "GetTorsionPoints" or "GetTorsionStructure");
        Assert.Single(CalculationCatalog.All, o => o.Id == "Q.TorsionPoints");
        Assert.Single(CalculationCatalog.All, o => o.Id == "Q.TorsionStructure");
    }

    [Theory]
    [InlineData("Q.TorsionPoints")]
    [InlineData("Q.TorsionStructure")]
    public async Task ExplorerTorsionActionsPassCancellationBeforeStartingPreparation(string id)
    {
        // A singular curve would otherwise fail preparation. Cancelling at the
        // computation boundary must take precedence and retain the caller's token.
        using var cancellation = new CancellationTokenSource();
        var request = new CalculationRequest(id, "y^2 = x^3", new());
        var error = await Assert.ThrowsAnyAsync<OperationCanceledException>(() =>
            CalculationEngine.ExecuteAsync(request, update =>
            {
                if (update.Message.StartsWith("Computing", StringComparison.Ordinal)) cancellation.Cancel();
            }, cancellation.Token));
        Assert.True(cancellation.IsCancellationRequested);
        Assert.Equal(cancellation.Token, error.CancellationToken);
    }

    [Theory]
    [InlineData("Multiply")]
    [InlineData("RationalPoints")]
    [InlineData("IntegralPoints")]
    public async Task ExplorerLegacyActionsRemainCancellable(string name)
    {
        var operation = Assert.Single(CalculationCatalog.All, o =>
            o.Context == CalculationContext.RationalCurve && o.Member?.Name == name);
        using var cancellation = new CancellationTokenSource();
        var request = new CalculationRequest(operation.Id, "y^2 + y = x^3 - x", new());
        var error = await Assert.ThrowsAnyAsync<OperationCanceledException>(() =>
            CalculationEngine.ExecuteAsync(request, update =>
            {
                if (update.Message.StartsWith("Computing", StringComparison.Ordinal)) cancellation.Cancel();
            }, cancellation.Token));
        Assert.True(cancellation.IsCancellationRequested);
        Assert.Equal(cancellation.Token, error.CancellationToken);
    }

    [Fact]
    public async Task ExplorerKeepsTheLegacySquareActionAndItsExactOutValue()
    {
        var operation = Assert.Single(CalculationCatalog.All, o =>
            o.Member?.DeclaringType == typeof(BigRational) && o.Member.Name == nameof(BigRational.IsSquare));
        Assert.Equal("BigRational.IsSquare(BigRational,BigRational&)", operation.Id);
        var request = new CalculationRequest(operation.Id, "y^2 = x^3 - x", new() { ["r"] = "9/4" });
        var output = await CalculationEngine.ExecuteAsync(request);
        Assert.Contains("Succeeded: True", output);
        Assert.Contains("s: 3/2", output);

        using var cancellation = new CancellationTokenSource();
        var error = await Assert.ThrowsAnyAsync<OperationCanceledException>(() =>
            CalculationEngine.ExecuteAsync(request, update =>
            {
                if (update.Message.StartsWith("Computing", StringComparison.Ordinal)) cancellation.Cancel();
            }, cancellation.Token));
        Assert.Equal(cancellation.Token, error.CancellationToken);
    }

    [Theory]
    [InlineData("1/2; 2/3; -1; 3/2")]
    [InlineData("1/5")]
    public void ExplorerFieldElementPreparationHonorsPreCancellation(string text)
    {
        var field = new FiniteField(5, new BigInteger[] { 2, 0, 1 });
        var token = new CancellationToken(true);
        Cancelled(token, () => CalculationInput.Element(field, text, token));
        Cancelled(token, () => CalculationInput.Read(typeof(FiniteFieldElement), "a",
            new Dictionary<string, string> { ["a"] = text }, field, token));
        var curve = new EllipticCurveFq(field, 0, 0, 0, 1, 0);
        Cancelled(token, () => CalculationInput.Read(typeof(EllipticCurvePointFq), "P",
            new Dictionary<string, string> { ["P.infinity"] = "False", ["P.x"] = text, ["P.y"] = "0" }, curve, token));
    }

    [Fact]
    public void ExplorerFieldElementPreparationPreservesExactHornerEvaluationAndLegacyCall()
    {
        var field = new FiniteField(5, new BigInteger[] { 2, 0, 1 });
        const string text = "1/2; 2/3; -1; 3/2";
        // In F_25, t^2=3, so this cubic with rational coefficients equals t.
        Assert.Equal(field.Generator, CalculationInput.Element(field, text, default));
        Assert.Equal(field.Generator, CalculationInput.Element(field, text));
        Assert.Equal(field.Generator, CalculationInput.Read(typeof(FiniteFieldElement), "a",
            new Dictionary<string, string> { ["a"] = text }, field, default));
    }

    private static IEnumerable<EllipticCurvePoint> Search(EllipticCurveQ curve, bool integral, CancellationToken token)
        => integral ? curve.IntegralPoints(10, token) : curve.RationalPoints(10, 10, token);

    private static void Cancelled(CancellationToken token, Action operation)
    {
        var error = Assert.ThrowsAny<OperationCanceledException>(operation);
        Assert.Equal(token, error.CancellationToken);
    }
}

using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class SaturationTests
{
    [Theory]
    [InlineData("389a1", 2), InlineData("5077a1", 2), InlineData("389a1", 3)]
    public void EnlargesMixedRankTwoAndThreeSubgroups(string label, int prime)
    {
        var data = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture(label)); var e = data.GlobalMinimalModel; var basis = data.MordellWeilGenerators.ToArray();
        // Determinant p^r and mixed basis vectors force tests of nontrivial F_p combinations.
        var input = basis.Select(p => e.Multiply(p, prime)).ToArray(); input[0] = e.Add(input[0], input[1]);
        var result = e.Saturate(input, new[] { prime });
        Assert.True(result.IsComplete, result.Reason); Assert.Equal(BigInteger.Pow(prime, basis.Length), result.IndexGain);
        var actual = e.Regulator(result.Generators); var expected = e.Regulator(basis);
        Assert.True(actual.LowerBound <= expected.UpperBound && expected.LowerBound <= actual.UpperBound);
    }
    [Theory]
    [InlineData(2), InlineData(3)]
    public void TorsionTranslatesAreIncluded(int prime)
    {
        var e = new EllipticCurveQ(0, 0, 0, -25, 0); var p = new EllipticCurvePoint(new BigRational(-4), new BigRational(6)); var t = new EllipticCurvePoint(0, 0);
        var result = e.Saturate(new[] { e.Add(e.Multiply(p, prime), t) }, new[] { prime });
        Assert.True(result.IsComplete, result.Reason); Assert.Equal(new BigInteger(prime), result.IndexGain);
        var difference = e.Add(result.Generators[0], e.Negate(p)); var sum = e.Add(result.Generators[0], p);
        Assert.True(e.TorsionPoints.Contains(difference) || e.TorsionPoints.Contains(sum));
    }
    [Fact]
    public void SeveralPrimePowersAndNonminimalModelsPreserveIndices()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0); var map = e.ChangeModel(new BigRational(2, 3), -4, 2, 7);
        var result = map.Target.Saturate(new[] { map.Map(e.Multiply(p, 12)) }, new[] { 3, 2, 2 });
        Assert.True(result.IsComplete, result.Reason); Assert.Equal(new BigInteger(12), result.IndexGain); Assert.Equal(new[] { 2, 3 }, result.CertifiedPrimes);
        var q = map.MapBack(result.Generators[0]); Assert.True(q.Equals(p) || q.Equals(e.Negate(p)));
        var atTwoOnly = e.Saturate(new[] { e.Multiply(p, 6) }, new[] { 2 });
        Assert.True(atTwoOnly.IsComplete); Assert.Equal(new BigInteger(2), atTwoOnly.IndexGain); Assert.Equal(new[] { 2 }, atTwoOnly.CertifiedPrimes);
    }
    [Fact]
    public void LimitsPreserveVerifiedEnlargementsAndEarlierPrimes()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        var limited = e.Saturate(new[] { e.Multiply(p, 4) }, new[] { 2 }, new SaturationOptions { MaxEnlargements = 1 });
        Assert.False(limited.IsComplete); Assert.Equal(new BigInteger(2), limited.IndexGain); Assert.Empty(limited.CertifiedPrimes); Assert.Equal(new[] { 2 }, limited.UnresolvedPrimes);
        Assert.True(e.Multiply(limited.Generators[0], 2).Equals(e.Multiply(p, 4)));
        var degree = e.Saturate(new[] { e.Multiply(p, 3) }, new[] { 2, 3 }, new SaturationOptions { MaxDivisionDegree = 4 });
        Assert.False(degree.IsComplete); Assert.Equal(new[] { 2 }, degree.CertifiedPrimes); Assert.Equal(new[] { 3 }, degree.UnresolvedPrimes); Assert.Equal(BigInteger.One, degree.IndexGain);
        Assert.True(e.Saturate(Array.Empty<EllipticCurvePoint>(), new[] { 2, 3 }).IsComplete);
        Assert.True(e.Saturate(new[] { p }, Array.Empty<int>()).IsComplete);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.Saturate(new[] { p }, new[] { 4 }));
        Assert.ThrowsAny<OperationCanceledException>(() => e.Saturate(new[] { p }, new[] { 2 }, cancellationToken: new CancellationToken(true)));
    }
    [Fact]
    public void SevenDivisionRecoversKnownGenerator()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        var result = e.Saturate(new[] { e.Multiply(p, 7) }, new[] { 7 });
        Assert.True(result.IsComplete, result.Reason); Assert.Equal(new BigInteger(7), result.IndexGain);
    }

    [Fact]
    public void DivisionPolynomialCoordinatesAgreeWithExactGroupLaw()
    {
        foreach (var e in new[] { new EllipticCurveQ(0, 0, 0, -1, 1), new EllipticCurveQ(0, 0, 0, -25, 0) })
        {
            var p = e.A6.IsZero ? new EllipticCurvePoint(-4, 6) : new EllipticCurvePoint(0, 1);
            var budget = new DescentBudget(new RankComputationOptions { MaxDescentWork = 2000000 }, default);
            var polynomials = new DivisionPolynomials(e.A4, e.A6, budget);
            foreach (int n in new[] { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13 })
            {
                var equation = polynomials.Multiplication(n);
                var target = e.Multiply(p, n);
                var denominator = DescentPolynomial.Evaluate(equation.denominator, p.X);
                if (target.IsInfinity) Assert.True(denominator.IsZero);
                else Assert.Equal(target.X, DescentPolynomial.Evaluate(equation.numerator, p.X) / denominator);
            }
        }
    }
}

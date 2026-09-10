using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class ExtendedArithmeticTests
{
    [Fact]
    public void ModelMapsPreserveGroupLawAndInvert()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        foreach (var u in new BigRational[] { -2, new BigRational(2, 3), 1 })
        {
            var map = e.ChangeModel(u, 5, -3, 7);
            for (int n = -5; n <= 5; n++)
            {
                var q = e.Multiply(p, n); var mapped = map.Map(q);
                Assert.True(map.Target.IsOnCurve(mapped));
                Assert.Equal(q, map.MapBack(mapped)); Assert.Equal(q, map.Inverse().Map(mapped));
                Assert.Equal(mapped, map.Target.Multiply(map.Map(p), n));
            }
            Assert.True(e.TryGetIsomorphism(map.Target, out var found));
            Assert.Equal(p, found.MapBack(found.Map(p)));
        }
        Assert.False(e.TryGetIsomorphism(new EllipticCurveQ(0, 0, 0, 1, 1), out _));
    }

    [Fact]
    public void Lmfdb37a1HasMatchingCertifiedHeightsAndPeriods()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        AssertClose(e.CanonicalHeight(p), 0.051111408239968840235886099757);
        AssertClose(e.Regulator(new[] { p }), 0.051111408239968840235886099757);
        var periods = e.GetPeriods();
        AssertClose(periods.RealPeriod, 5.9869172924639192596640199589);
        AssertClose(periods.Area, 7.3381327407895767390707210033);
        Assert.Equal(0, periods.SecondPeriodRealPart.LowerBound);
        for (int n = 2; n <= 4; n++) AssertClose(e.CanonicalHeight(e.Multiply(p, n)), n * n * 0.051111408239968840235886099757);
        var changed = e.ChangeModel(new BigRational(2, 3), 4, -2, 3);
        AssertClose(changed.Target.CanonicalHeight(changed.Map(p)), 0.05111140823996896884);
        AssertClose(changed.Target.GetPeriods().RealPeriod, periods.RealPeriod.Approximation);
        Assert.Throws<ArithmeticException>(() => e.Regulator(new[] { p, e.Multiply(p, 2) }));
    }

    [Theory]
    [InlineData(0, 0, 1, -1, 0, 37, "I1", 1, 1)]
    [InlineData(0, 0, 0, -1, 0, 2, "III", 2, 5)]
    [InlineData(0, 0, 1, 0, 0, 3, "II", 1, 3)]
    public void LocalDataExamples(int a1, int a2, int a3, int a4, int a6, int prime, string kodaira, int cp, int f)
    {
        var e = new EllipticCurveQ(a1, a2, a3, a4, a6); var d = e.GetLocalData(prime);
        Assert.Equal(kodaira, d.KodairaSymbol); Assert.Equal(cp, d.TamagawaNumber); Assert.Equal(f, d.ConductorValuation);
        Assert.Equal(e.Conductor, e.LocalData.Aggregate(BigInteger.One, (n, l) => n * BigInteger.Pow(l.Prime, l.ConductorValuation)));
        Assert.Equal(ReductionType.Good, e.GetLocalData(101).ReductionType);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetLocalData(4));
    }

    [Fact]
    public void RealComputationsRespectLimitsAndCancellation()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        Assert.Throws<ArithmeticException>(() => e.CanonicalHeight(p, new RealComputationOptions { MaxIterations = 1 }));
        Assert.Throws<ArithmeticException>(() => e.GetPeriods(new RealComputationOptions { MaxRootWork = 0 }));
        Assert.ThrowsAny<OperationCanceledException>(() => e.CanonicalHeight(p, cancellationToken: new CancellationToken(true)));
        Assert.ThrowsAny<OperationCanceledException>(() => e.GetPeriods(cancellationToken: new CancellationToken(true)));
        Assert.Equal(0, e.CanonicalHeight(EllipticCurvePoint.Infinity).LowerBound);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    public void SaturationRecoversKnownGenerator(int prime)
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        var result = e.Saturate(new[] { e.Multiply(p, prime) }, new[] { prime });
        Assert.True(result.IsComplete, result.Reason); Assert.Equal(new BigInteger(prime), result.IndexGain);
        Assert.True(result.Generators[0].Equals(p) || result.Generators[0].Equals(e.Negate(p)));
    }

    [Fact]
    public void SaturationReportsUnprovedInputsAndLimits()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = new EllipticCurvePoint(0, 0);
        Assert.False(e.Saturate(new[] { p, e.Multiply(p, 2) }, new[] { 2 }).IndependenceCertified);
        var partial = e.Saturate(new[] { e.Multiply(p, 2) }, new[] { 2, 3 }, new SaturationOptions { MaxWork = 0 });
        Assert.False(partial.IsComplete); Assert.Equal(new[] { 2, 3 }, partial.UnresolvedPrimes);
        Assert.Equal(BigInteger.One, partial.IndexGain);
    }

    [Fact]
    public void LocalHeightsSumToCanonicalHeightIncludingDenominators()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = e.Multiply(new EllipticCurvePoint(0, 0), 5);
        var infinity = e.ArchimedeanHeight(p); var atTwo = e.LocalHeight(p, 2); var at37 = e.LocalHeight(p, 37); var height = e.CanonicalHeight(p);
        Assert.True(atTwo.LowerBound > 0);
        Assert.True(infinity.LowerBound + atTwo.LowerBound + at37.LowerBound <= height.UpperBound);
        Assert.True(infinity.UpperBound + atTwo.UpperBound + at37.UpperBound >= height.LowerBound);
        Assert.Equal(0, e.LocalHeight(p, 3).LowerBound); Assert.Equal(0, e.LocalHeight(p, 3).UpperBound);
        AssertClose(e.NaiveHeight(p), Math.Log(4));
        Assert.Throws<ArgumentException>(() => e.LocalHeight(EllipticCurvePoint.Infinity, 2));
        var torsion = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.All(torsion.TorsionPoints, q => Assert.True(torsion.CanonicalHeight(q).Contains(0)));
        var q = e.Negate(p); var pairing = e.HeightPairing(p, q);
        Assert.True(pairing.LowerBound <= -height.LowerBound && pairing.UpperBound >= -height.UpperBound);
    }

    [Theory]
    [InlineData(30, 512)]
    [InlineData(60, 1024)]
    public void HigherPrecisionEnclosuresContainIndependentReferences(int digits, int bits)
    {
        var options = new RealComputationOptions { DecimalDigits = digits, PrecisionBits = bits };
        var periodRow = File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "periods.csv")).First(x => x.StartsWith("0,0,1,-1,0,"));
        var heightRow = File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "heights.csv")).First(x => x.StartsWith("0,0,1,-1,0,"));
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var v = heightRow.Split(',');
        var p = new EllipticCurvePoint(ExtendedReferenceTests.Parse(v[5]), ExtendedReferenceTests.Parse(v[6]));
        var height = e.CanonicalHeight(p, options);
        ExtendedReferenceTests.ReferenceIn(height, v[7], digits);
        ExtendedReferenceTests.ReferenceIn(e.GetPeriods(options).RealPeriod, periodRow.Split(',')[7], digits);
        Assert.True(height.Width <= new BigRational(1, BigInteger.Pow(10, digits)));
    }

    [Fact]
    public void ExceptionalJInvariantsAndInvalidModelsAreHandled()
    {
        foreach (var e in new[] { new EllipticCurveQ(0, 0, 0, -1, 0), new EllipticCurveQ(0, 0, 0, 0, 1) })
        {
            var changed = e.ChangeModel(new BigRational(-2, 3), 7, 4, -5).Target;
            Assert.True(e.TryGetIsomorphism(changed, out var map)); Assert.Equal(changed, map.Target);
            Assert.True(e.GetMinimalModelIsomorphism().MapBack(EllipticCurvePoint.Infinity).IsInfinity);
        }
        Assert.False(new EllipticCurveQ(0, 0, 0, -1, 0).TryGetIsomorphism(new EllipticCurveQ(0, 0, 0, -4, 0), out _));
        Assert.False(new EllipticCurveQ(0, 0, 0, 0, 1).TryGetIsomorphism(new EllipticCurveQ(0, 0, 0, 0, 2), out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => new EllipticCurveQ(0, 0, 0, 0, 1).ChangeModel(0, 0, 0, 0));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).GetPeriods());
    }

    private static void AssertClose(RealEnclosure result, double expected)
    {
        Assert.True(result.Width <= new BigRational(1, BigInteger.Pow(10, 12)));
        Assert.InRange(result.Approximation, expected - 1e-12, expected + 1e-12);
    }
}

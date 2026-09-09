using System.Globalization;
using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class AnalyticRankTests
{
    private static string[] Rows => File.ReadAllLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "analytic.csv"));
    public static IEnumerable<object[]> References() => Rows.Select(row => new object[] { row });
    private static EllipticCurveQ Curve(string[] a) => new(new(BigInteger.Parse(a[0])), new(BigInteger.Parse(a[1])),
        new(BigInteger.Parse(a[2])), new(BigInteger.Parse(a[3])), new(BigInteger.Parse(a[4])));

    [Theory]
    [MemberData(nameof(References))]
    public void MatchesPariAnalyticRankAndLeadingDerivative(string row)
    {
        var a = row.Split(',');
        var result = Curve(a).EstimateAnalyticRank();
        int rank = int.Parse(a[7]);
        Assert.True(result.EstimatedRank == rank, $"Expected {rank}, got {result}: {result.Reason}");
        Assert.Equal(BigInteger.Parse(a[5]), result.Conductor);
        Assert.Equal(int.Parse(a[6]), result.RootNumber);
        double leading = double.Parse(a[8], CultureInfo.InvariantCulture);
        Assert.InRange(Math.Abs(result.Derivatives[rank] - leading), 0, 1e-9 * Math.Max(1, leading));
        for (int k = 0; k < rank; k++) Assert.InRange(Math.Abs(result.Derivatives[k]), 0, 1e-9);
        if (rank <= 1)
        {
            Assert.Equal(rank, result.ProvenRank);
            Assert.Equal(AnalyticRankStatus.Certified, result.Status);
            Assert.True(result.CertifiedLeadingLowerBound > BigRational.Zero);
            var reference = DecimalRational(a[8]);
            Assert.True(result.CertifiedLeadingLowerBound <= reference);
            Assert.True(result.CertifiedLeadingUpperBound >= reference);
        }
        else
        {
            Assert.Null(result.ProvenRank);
            Assert.Equal(AnalyticRankStatus.NumericalEstimate, result.Status);
        }
    }

    private static BigRational DecimalRational(string value)
    {
        var parts = value.Replace(" ", "").Split('E');
        int exponent = parts.Length == 2 ? int.Parse(parts[1]) : 0;
        value = parts[0];
        int dot = value.IndexOf('.');
        int scale = (dot < 0 ? 0 : value.Length - dot - 1) - exponent;
        var numerator = BigInteger.Parse(value.Replace(".", ""));
        return scale >= 0 ? new BigRational(numerator, BigInteger.Pow(10, scale)) : new BigRational(numerator * BigInteger.Pow(10, -scale));
    }

    [Fact]
    public void OrdinaryDerivativesThroughOrderEightMatchIndependentLFunctionEvaluation()
    {
        var results = Rows.Take(2).Select(row => Curve(row.Split(',')).EstimateAnalyticRank(new()
            { MaxDerivativeOrder = 8, CertifyLowRanks = false })).ToArray();
        foreach (string row in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "derivatives.csv")))
        {
            var a = row.Split(',');
            var result = results[int.Parse(a[0]) - 1];
            int k = int.Parse(a[1]);
            Assert.NotEqual(AnalyticRankStatus.Inconclusive, result.Status);
            double expected = double.Parse(a[2].Replace(" ", ""), CultureInfo.InvariantCulture);
            Assert.InRange(Math.Abs(result.Derivatives[k] - expected), 0, Math.Max(2e-12, result.EstimatedErrors[k]));
        }
    }

    [Fact]
    public void RigorousSpecialFunctionIntervalsContainPariValuesAcrossTheRange()
    {
        foreach (string row in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "special-functions.csv")))
        {
            var a = row.Split(',');
            var x = DyadicInterval.Fraction(BigInteger.Parse(a[0]), BigInteger.Parse(a[1]));
            var exp = DyadicInterval.ExpNegative(x);
            var expRef = DecimalRational(a[2]);
            Assert.True(exp.LowerRational <= expRef && exp.UpperRational >= expRef, $"exp(-x): {row}");
            var e1 = AnalyticCertificate.ExponentialIntegral(x, default);
            var e1Ref = DecimalRational(a[3]);
            Assert.True(e1.LowerRational <= e1Ref && e1.UpperRational >= e1Ref, $"E1(x): {row}");
        }
    }

    [Fact]
    public void CoefficientsMatchPariIncludingBadPrimesAndPrimePowers()
    {
        var curves = Rows.Select(row => Curve(row.Split(',')).GlobalMinimalModel).ToArray();
        var coefficients = curves.Select(e => AnalyticCoefficients.Compute(e, 100, 10000, default)).ToArray();
        foreach (string row in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "coefficients.csv")))
        {
            var a = row.Split(',').Select(int.Parse).ToArray();
            Assert.Equal(a[2], coefficients[a[0] - 1][a[1]]);
        }
    }

    [Fact]
    public void AllSmallDerivativesAreInconclusiveRatherThanAProvedLowerBound()
    {
        var result = new EllipticCurveQ(0, 1, 1, -2, 0).EstimateAnalyticRank(new() { MaxDerivativeOrder = 1 });
        Assert.Equal(AnalyticRankStatus.Inconclusive, result.Status);
        Assert.Null(result.EstimatedRank);
        Assert.Null(result.ProvenRank);
    }

    [Fact]
    public void SupportsRationalAndNonminimalModels()
    {
        var minimal = new EllipticCurveQ(0, 0, 1, -1, 0);
        // x=X/4, y=Y/8, then the inverse scaling.
        var integral = new EllipticCurveQ(0, 0, 8, -16, 0);
        var rational = new EllipticCurveQ(0, 0, new(1, 8), new(-1, 16), 0);
        var baseline = minimal.EstimateAnalyticRank();
        foreach (var e in new[] { integral, rational })
        {
            var result = e.EstimateAnalyticRank();
            Assert.Equal(1, result.ProvenRank);
            Assert.Equal(baseline.Derivatives, result.Derivatives);
            Assert.Equal(minimal.RootNumber, e.RootNumber);
        }
    }

    [Fact]
    public void LowRankNumericalEstimateDoesNotImplyProof()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        var result = e.EstimateAnalyticRank(new() { CertifyLowRanks = false });
        Assert.Equal(1, result.EstimatedRank);
        Assert.Null(result.ProvenRank);
        Assert.Equal(AnalyticRankStatus.NumericalEstimate, result.Status);
        // Too short an interval sum must retain its tail instead of claiming certainty.
        var shortCertificate = new EllipticCurveQ(0, 0, 0, -25, 0).EstimateAnalyticRank(new() { MaxCertificationTerms = 1 });
        Assert.Null(shortCertificate.ProvenRank);
        Assert.True(shortCertificate.CertifiedLeadingLowerBound < BigRational.Zero);
        Assert.True(shortCertificate.CertifiedLeadingUpperBound > BigRational.Zero);
    }

    [Fact]
    public void ResourceLimitsProduceAnUnknownRank()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        foreach (var options in new AnalyticRankOptions[]
        {
            new() { MaxTerms = 1 }, new() { MaxPointCountingWork = 1 }, new() { MaxIntegrationEvaluations = 1 }
        })
        {
            var result = e.EstimateAnalyticRank(options);
            Assert.Equal(AnalyticRankStatus.Inconclusive, result.Status);
            Assert.Null(result.EstimatedRank);
        }
    }

    [Fact]
    public void RejectsInvalidOptionsSingularCurvesAndCancellation()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        foreach (var options in new AnalyticRankOptions[]
        {
            new() { MaxDerivativeOrder = -1 }, new() { MaxDerivativeOrder = 9 }, new() { ZeroTolerance = double.NaN },
            new() { ZeroTolerance = 0 }, new() { MaxTerms = 0 }, new() { MaxTerms = int.MaxValue },
            new() { MaxPointCountingWork = 0 }, new() { MaxIntegrationEvaluations = 0 }, new() { MaxCertificationTerms = -1 }
        }) Assert.Throws<ArgumentOutOfRangeException>(() => e.EstimateAnalyticRank(options));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).EstimateAnalyticRank());
        var canceled = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => e.EstimateAnalyticRank(cancellationToken: canceled));
        Assert.Throws<OperationCanceledException>(() => e.GetRootNumber(canceled));
    }

    [Fact]
    public void IntervalRoundingHandlesNegativeValuesAndTranscendentals()
    {
        var x = DyadicInterval.Fraction(-1, 3);
        Assert.True(x.LowerRational <= new BigRational(-1, 3) && x.UpperRational >= new BigRational(-1, 3));
        var squared = x * x;
        Assert.True(squared.LowerRational <= new BigRational(1, 9) && squared.UpperRational >= new BigRational(1, 9));
        var pi = DecimalRational("3.14159265358979323846264338327950288419716939937510582097494");
        Assert.True(DyadicInterval.Pi().LowerRational <= pi && DyadicInterval.Pi().UpperRational >= pi);
        var exp = DyadicInterval.ExpNegative(DyadicInterval.Integer(1));
        var expRef = DecimalRational("0.367879441171442321595523770161460867445811131031767834507837");
        Assert.True(exp.LowerRational <= expRef && exp.UpperRational >= expRef);
        var e1 = AnalyticCertificate.ExponentialIntegral(DyadicInterval.Integer(1), default);
        var e1Ref = DecimalRational("0.219383934395520273677163775460121649031047293406908207577979");
        Assert.True(e1.LowerRational <= e1Ref && e1.UpperRational >= e1Ref);
    }
}

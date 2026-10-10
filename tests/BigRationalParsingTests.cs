using System.Globalization;
using System.Numerics;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class BigRationalParsingTests
{
    [Theory]
    [InlineData(" +12 ", "12", "1")]
    [InlineData("−12 / −3", "4", "1")]
    [InlineData(".5", "1", "2")]
    [InlineData("5.", "5", "1")]
    [InlineData("5.e1", "50", "1")]
    [InlineData("1E-3", "1", "1000")]
    [InlineData("1e−3", "1", "1000")]
    [InlineData("+.5e+2 / -1.25e1", "-4", "1")]
    [InlineData("-1.2E+3/.04", "-30000", "1")]
    [InlineData("1e2/0.4", "250", "1")]
    [InlineData("0.1", "1", "10")]
    [InlineData("\t-0.00 / +7e2\r\n", "0", "1")]
    [InlineData("1e00000000000000000000000000003", "1000", "1")]
    public void ParsesTheExistingExactInvariantGrammar(string text, string numerator, string denominator)
    {
        var expected = new BigRational(BigInteger.Parse(numerator, CultureInfo.InvariantCulture),
            BigInteger.Parse(denominator, CultureInfo.InvariantCulture));
        Assert.True(BigRational.TryParse(text, out var parsed));
        Assert.Equal(expected, parsed);
        Assert.Equal(expected, BigRational.Parse(text));
        Assert.True(RationalText.TryParse(text, out var editorValue));
        Assert.Equal(expected, editorValue);
    }

    [Theory]
    [InlineData(null)]
    [InlineData("")]
    [InlineData(" \t\r\n ")]
    [InlineData("+")]
    [InlineData(".")]
    [InlineData("-.")]
    [InlineData("1e")]
    [InlineData("1e-")]
    [InlineData("1e+-2")]
    [InlineData("1 2")]
    [InlineData("- 1")]
    [InlineData("1,25")]
    [InlineData("1,000")]
    [InlineData("١")]
    [InlineData("1/0")]
    [InlineData("0/0.000e3")]
    [InlineData("1/")]
    [InlineData("/1")]
    [InlineData("1/2/3")]
    [InlineData("(1/2)")]
    [InlineData("NaN")]
    [InlineData("Infinity")]
    [InlineData("1_000")]
    public void FailedParsingClearsTheResultAndParseReportsTheFailure(string text)
    {
        var value = new BigRational(123);
        Assert.False(BigRational.TryParse(text, out value));
        Assert.Equal(BigRational.Zero, value);
        Assert.False(RationalText.TryParse(text, out var editorValue));
        Assert.Equal(BigRational.Zero, editorValue);
        if (text == null) Assert.Throws<ArgumentNullException>(() => BigRational.Parse(text));
        else Assert.Throws<FormatException>(() => BigRational.Parse(text));
    }

    [Theory]
    [InlineData("ru-RU")]
    [InlineData("fr-FR")]
    [InlineData("en-US")]
    public void CultureDoesNotChangeTheNumberGrammarOrExactResult(string culture)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal(new BigRational(5, 2), BigRational.Parse("1.25/0.5"));
            Assert.False(BigRational.TryParse("1,25", out _));
            Assert.False(BigRational.TryParse("1,000", out _));
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }

    [Fact]
    public void LengthAndExponentLimitsIncludeTheirExactBoundaries()
    {
        var length = new RationalParsingOptions { MaxTextLength = 3 };
        Assert.Equal(new BigRational(1, 2), BigRational.Parse("1/2", length));
        Rejected(" 1/2", length);
        var exponent = new RationalParsingOptions { MaxExponent = 3 };
        Assert.Equal(new BigRational(1000), BigRational.Parse("1e3", exponent));
        Assert.Equal(new BigRational(1, 1000), BigRational.Parse("1e-3", exponent));
        Rejected("1e4", exponent);
        Rejected("1e-4", exponent);
        Rejected("0e4", exponent);
        exponent.MaxExponent = 0;
        Assert.Equal(BigRational.One, BigRational.Parse("1e-0", exponent));
        Rejected("1e1", exponent);
    }

    [Fact]
    public void FinalValueLimitUsesNormalizedNumeratorAndDenominator()
    {
        var limits = new RationalParsingOptions { MaxValueBits = 3 };
        Assert.Equal(new BigRational(7, 4), BigRational.Parse("7/4", limits));
        Assert.Equal(new BigRational(5, 4), BigRational.Parse("1.25", limits));
        Assert.Equal(new BigRational(3, 2), BigRational.Parse("1.5000", limits));
        Rejected("8", limits);
        Rejected("1/8", limits);
        Rejected("9/4", limits);
    }

    [Fact]
    public void IntermediateLimitBoundsPowersAndComponentsBeforeFinalReduction()
    {
        var limits = new RationalParsingOptions { MaxIntermediateBits = 4 };
        Assert.Equal(new BigRational(10), BigRational.Parse("1e1", limits));
        Assert.Equal(new BigRational(1, 10), BigRational.Parse("1e-1", limits));
        Assert.Equal(new BigRational(3, 2), BigRational.Parse("1.5", limits));
        Rejected("16", limits);
        Rejected("1e2", limits);
        Rejected("1e-2", limits);
        Rejected("1.25", limits);
        Rejected("3/0.1", limits);
    }

    [Fact]
    public void ZeroAndTrailingZerosDoNotRequireUnnecessaryLargePowers()
    {
        var limits = new RationalParsingOptions
        { MaxValueBits = 1, MaxIntermediateBits = 1, MaxExponent = int.MaxValue };
        Assert.Equal(BigRational.One, BigRational.Parse("10000e-4", limits));
        Assert.Equal(BigRational.One, BigRational.Parse("0.000010000e5", limits));
        Assert.Equal(BigRational.Zero, BigRational.Parse("0.000e2147483647", limits));
        Rejected("1e2147483647", limits);
        Rejected("1e-2147483647", limits);
        Rejected("1e2147483648", limits);
        Rejected("0e-2147483648", limits);
    }

    [Fact]
    public void ZeroFractionsStillValidateTheirDenominatorBeforeNormalization()
    {
        var limits = new RationalParsingOptions { MaxValueBits = 1, MaxIntermediateBits = 64 };
        Assert.Equal(BigRational.Zero, BigRational.Parse("0/123456789", limits));
        Rejected("0/0", limits);
        Rejected("0/not-a-number", limits);
        Rejected("0/1e5000", limits);
    }

    [Fact]
    public void LargeFractionComponentsMayCancelWithinTheExistingEditorBudgets()
    {
        var component = (BigInteger.One << 32768).ToString(CultureInfo.InvariantCulture);
        var text = component + "/" + component;
        Assert.InRange(text.Length, 1, RationalText.MaxTextLength);
        Assert.Equal(BigRational.One, BigRational.Parse(text));
        Assert.True(RationalText.TryParse(text, out var editorValue));
        Assert.Equal(BigRational.One, editorValue);
    }

    [Fact]
    public void ConfigurableLimitsSupportLargeCanonicalStringRoundTrips()
    {
        var original = new BigRational((BigInteger.One << 50000) + 17, (BigInteger.One << 48000) + 11);
        var text = original.ToString();
        Assert.False(BigRational.TryParse(text, out var defaultValue));
        Assert.Equal(BigRational.Zero, defaultValue);
        var limits = new RationalParsingOptions
        { MaxTextLength = 40000, MaxValueBits = 60000, MaxIntermediateBits = 120000 };
        Assert.Equal(original, BigRational.Parse(text, limits));
        Assert.True(BigRational.TryParse(text, limits, out var restored));
        Assert.Equal(original, restored);
    }

    [Theory]
    [InlineData(nameof(RationalParsingOptions.MaxTextLength))]
    [InlineData(nameof(RationalParsingOptions.MaxExponent))]
    [InlineData(nameof(RationalParsingOptions.MaxValueBits))]
    [InlineData(nameof(RationalParsingOptions.MaxIntermediateBits))]
    public void InvalidOptionsThrowInsteadOfBecomingAParseFailure(string property)
    {
        var limits = new RationalParsingOptions();
        typeof(RationalParsingOptions).GetProperty(property)!.SetValue(limits,
            property == nameof(RationalParsingOptions.MaxExponent) ? -1 : 0);
        var error = Assert.Throws<ArgumentOutOfRangeException>(() => BigRational.TryParse("1", limits, out _));
        Assert.Equal(property, error.ParamName);
        Assert.Throws<ArgumentOutOfRangeException>(() => BigRational.Parse("1", limits));
    }

    [Fact]
    public void LibraryTextParsingDoesNotCreateDuplicateMathematicalMenuActions()
    {
        Assert.DoesNotContain(CalculationCatalog.All, action => action.Member?.Name is "Parse" or "TryParse"
            || action.Member?.Name.EndsWith("Checked", StringComparison.Ordinal) == true);
        Assert.Single(CalculationCatalog.All, action => action.Id == "EllipticCurveQ.Multiply(EllipticCurvePoint,BigInteger)");
    }

    private static void Rejected(string text, RationalParsingOptions limits)
    {
        Assert.False(BigRational.TryParse(text, limits, out var value));
        Assert.Equal(BigRational.Zero, value);
        Assert.Throws<FormatException>(() => BigRational.Parse(text, limits));
    }
}

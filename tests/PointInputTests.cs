using System.Numerics;
using EllipticCurves.Explorer.Computations;
using Xunit;

namespace EllipticCurves.Tests;

public class PointInputTests
{
    private static EllipticCurvePoint[] Parse(string text) =>
        (EllipticCurvePoint[])CalculationInput.ParseScalar(typeof(IReadOnlyList<EllipticCurvePoint>), text);

    [Theory]
    [InlineData("-523548280341848, -2806322695726774350150")]
    [InlineData("(-523548280341848, -2806322695726774350150)")]
    [InlineData("-523548280341848; -2806322695726774350150")]
    public void LargeCoordinatesRetainEveryDigit(string text)
    {
        var point = Assert.Single(Parse(text));
        Assert.Equal(new BigRational(BigInteger.Parse("-523548280341848")), point.X);
        Assert.Equal(new BigRational(BigInteger.Parse("-2806322695726774350150")), point.Y);
    }

    [Theory]
    [InlineData("1439725582641991678/9, -1727498601505665033867479240/27")]
    [InlineData("(1439725582641991678/9, -1727498601505665033867479240/27)")]
    public void FractionalCoordinatesRemainExact(string text)
    {
        var point = Assert.Single(Parse(text));
        Assert.Equal(new BigRational(BigInteger.Parse("1439725582641991678"), 9), point.X);
        Assert.Equal(new BigRational(BigInteger.Parse("-1727498601505665033867479240"), 27), point.Y);
    }

    [Theory]
    [InlineData("1.5,-2.5")]
    [InlineData("( 1.5, -2.5 )")]
    [InlineData("1.5; -2.5")]
    [InlineData("(1.5; -2.5)")]
    [InlineData("1.5e0, -25e-1")]
    public void DecimalCoordinatesAcceptBothPointSeparators(string text)
    {
        var point = Assert.Single(Parse(text));
        Assert.Equal(new BigRational(3, 2), point.X);
        Assert.Equal(new BigRational(-5, 2), point.Y);
    }

    [Fact]
    public void MixedListsAcceptWhitespaceAndInfinity()
    {
        const string text = " \r\n0, 2\r\n\t\r\n (1, 0) \n2; 0\n o \n";
        var expected = new[] { new EllipticCurvePoint(0, 2), new(1, 0), new(2, 0), EllipticCurvePoint.Infinity };
        Assert.Equal(expected, Parse(text));
        Assert.Equal(expected, (EllipticCurvePoint[])CalculationInput.ParseScalar(typeof(IEnumerable<EllipticCurvePoint>), text));
        Assert.Empty(Parse(" \r\n\t\n "));
    }

    [Theory]
    [InlineData("1,2,3")]
    [InlineData("(1,5,2,5)")]
    [InlineData("1,5; -2,5")]
    [InlineData("(1,5; -2,5)")]
    [InlineData("(1, 2")]
    [InlineData("1, 2)")]
    [InlineData("((1, 2))")]
    [InlineData("1,")]
    [InlineData(",2")]
    [InlineData("1/0, 2")]
    [InlineData("1; 2; 3")]
    public void MalformedCoordinatesAreRejected(string text) =>
        Assert.Throws<FormatException>(() => Parse(text));
}

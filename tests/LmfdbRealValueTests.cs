using System.Globalization;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class LmfdbRealValueTests
{
    [Theory]
    [InlineData(".-1")]
    [InlineData(".+1")]
    [InlineData("1.e 2")]
    [InlineData("1e+ 2")]
    [InlineData("1.2.3")]
    [InlineData("NaN")]
    [InlineData("Infinity")]
    [InlineData("")]
    [InlineData(null)]
    [InlineData("1e2147483648")]
    [InlineData("1e10001")]
    [InlineData("1e-10001")]
    public void RejectsMalformedDecimals(string text)
    {
        Assert.Throws<FormatException>(() => new LmfdbRealValue(text, null));
    }

    [Theory]
    [InlineData(".5", 1, 2)]
    [InlineData("-.5", -1, 2)]
    [InlineData("1.", 1, 1)]
    [InlineData("  +1.25e+2  ", 125, 1)]
    [InlineData("-12.5E-2", -1, 8)]
    public void ValidDecimalFormsPreserveTheirExactValue(string text, int numerator, int denominator)
    {
        var stored = new LmfdbRealValue(text, 97);
        Assert.Equal(new BigRational(numerator, denominator), stored.AsRational());
        Assert.Equal(text, stored.DecimalValue);
        Assert.Equal(double.Parse(text, CultureInfo.InvariantCulture), stored.Approximation);
        Assert.Equal(97, stored.StoredPrecisionBits);
    }
}

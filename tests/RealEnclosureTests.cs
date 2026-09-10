using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class RealEnclosureTests
{
    [Theory]
    [InlineData(-1)]
    [InlineData(1)]
    public void ApproximationPreservesFiniteValuesNearOverflow(int sign)
    {
        var value = new BigRational(sign * (5 * BigInteger.Pow(2, 1200) + 1), 7 * BigInteger.Pow(2, 176) + 1);
        var interval = new RealEnclosure(value, value);
        Assert.Equal(sign * (10.0 / 7) * Math.Pow(2, 1023), interval.Approximation);
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(1)]
    public void ApproximationRoundsVerySmallRationalsBeforeUnderflow(int sign)
    {
        var value = new BigRational(sign * (3 * BigInteger.Pow(2, 1200) + 1), 2 * BigInteger.Pow(2, 2275) + 1);
        Assert.Equal(sign * double.Epsilon, new RealEnclosure(value, value).Approximation);
    }

    [Fact]
    public void ApproximationUsesRoundToNearestEvenAtTies()
    {
        var scale = BigInteger.One << 53;
        Assert.Equal(1.0, RealEnclosure.ToDouble(new BigRational(scale + 1, scale)));
        Assert.Equal(Math.BitIncrement(Math.BitIncrement(1.0)), RealEnclosure.ToDouble(new BigRational(scale + 3, scale)));
        var subnormalDenominator = BigInteger.One << 1075;
        Assert.Equal(0L, BitConverter.DoubleToInt64Bits(RealEnclosure.ToDouble(new BigRational(1, subnormalDenominator))));
        Assert.Equal(long.MinValue, BitConverter.DoubleToInt64Bits(RealEnclosure.ToDouble(new BigRational(-1, subnormalDenominator))));
        Assert.Equal(2 * double.Epsilon, RealEnclosure.ToDouble(new BigRational(3, subnormalDenominator)));
        Assert.Equal(double.Epsilon, RealEnclosure.ToDouble(new BigRational((BigInteger.One << 100) + 1, BigInteger.One << 1175)));
    }

    [Fact]
    public void ApproximationRoundsAtTheNormalAndOverflowBoundaries()
    {
        var overflowMidpoint = ((BigInteger.One << 54) - 1) << 970;
        Assert.Equal(double.MaxValue, RealEnclosure.ToDouble(new BigRational(overflowMidpoint - 1)));
        Assert.Equal(double.PositiveInfinity, RealEnclosure.ToDouble(new BigRational(overflowMidpoint)));
        Assert.Equal(-double.MaxValue, RealEnclosure.ToDouble(new BigRational(1 - overflowMidpoint)));
        Assert.Equal(double.NegativeInfinity, RealEnclosure.ToDouble(new BigRational(-overflowMidpoint)));

        var normalMidpointNumerator = (BigInteger.One << 54) - 2;
        var denominator = BigInteger.One << 1076;
        Assert.Equal(BitConverter.Int64BitsToDouble(0x000fffffffffffff),
            RealEnclosure.ToDouble(new BigRational(normalMidpointNumerator - 1, denominator)));
        Assert.Equal(BitConverter.Int64BitsToDouble(0x0010000000000000),
            RealEnclosure.ToDouble(new BigRational(normalMidpointNumerator, denominator)));
    }

    [Fact]
    public void ApproximationRoundTripsFiniteDoubleBitPatterns()
    {
        var samples = new List<long> { 1, 2, 0x000fffffffffffff, 0x0010000000000000, 0x3fefffffffffffff,
            0x3ff0000000000000, 0x3ff0000000000001, 0x7fefffffffffffff };
        var random = new Random(20260910);
        for (int i = 0; i < 1000; i++) samples.Add(random.NextInt64(1, 0x7ff0000000000000));
        foreach (long bits in samples)
        {
            int exponent = (int)((bits >> 52) & 0x7ff);
            BigInteger mantissa = bits & 0x000fffffffffffff;
            if (exponent != 0) mantissa += BigInteger.One << 52;
            int power = exponent == 0 ? -1074 : exponent - 1023 - 52;
            var exact = power >= 0 ? new BigRational(mantissa << power) : new BigRational(mantissa, BigInteger.One << -power);
            Assert.Equal(bits, BitConverter.DoubleToInt64Bits(RealEnclosure.ToDouble(exact)));
            Assert.Equal(bits | long.MinValue, BitConverter.DoubleToInt64Bits(RealEnclosure.ToDouble(-exact)));
        }
        Assert.Equal(double.PositiveInfinity, RealEnclosure.ToDouble(new BigRational(BigInteger.One << 1024)));
        Assert.Equal(double.NegativeInfinity, RealEnclosure.ToDouble(new BigRational(-(BigInteger.One << 1024))));
        Assert.Equal(0.0, RealEnclosure.ToDouble(new BigRational(1, BigInteger.One << 4000)));
    }
}

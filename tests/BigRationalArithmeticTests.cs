using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class BigRationalArithmeticTests
{
    [Theory]
    [InlineData(14, 15, 25, 21, 10, 9)]
    [InlineData(-14, 15, 25, 21, -10, 9)]
    [InlineData(14, -15, -25, 21, 10, 9)]
    [InlineData(3, 1, 2, 9, 2, 3)]
    [InlineData(2, 9, 3, 1, 2, 3)]
    [InlineData(-7, 1, -3, 1, 21, 1)]
    public void MultiplicationCancelsCrossFactorsAndPreservesSigns(
        long a, long b, long c, long d, long numerator, long denominator)
    {
        AssertCanonicalEqual(new BigRational(numerator, denominator), new BigRational(a, b) * new BigRational(c, d));
    }

    [Theory]
    [InlineData(14, 15, 21, 25, 10, 9)]
    [InlineData(-14, 15, 21, 25, -10, 9)]
    [InlineData(14, 15, -21, 25, -10, 9)]
    [InlineData(-14, 15, -21, 25, 10, 9)]
    [InlineData(6, 35, 9, 14, 4, 15)]
    [InlineData(7, 1, -3, 1, -7, 3)]
    [InlineData(3, 1, 9, 2, 2, 3)]
    [InlineData(2, 9, 3, 1, 2, 27)]
    public void DivisionCancelsNumeratorAndDenominatorFactorsAndNormalizesSigns(
        long a, long b, long c, long d, long numerator, long denominator)
    {
        AssertCanonicalEqual(new BigRational(numerator, denominator), new BigRational(a, b) / new BigRational(c, d));
    }

    [Theory]
    [InlineData(1, 6, 1, 10, 4, 15, 1, 15)]
    [InlineData(1, 6, 1, 15, 7, 30, 1, 10)]
    [InlineData(2, 3, 4, 5, 22, 15, -2, 15)]
    [InlineData(1, 7, 6, 7, 1, 1, -5, 7)]
    [InlineData(3, 1, 2, 7, 23, 7, 19, 7)]
    [InlineData(2, 7, 3, 1, 23, 7, -19, 7)]
    [InlineData(-5, 1, 3, 1, -2, 1, -8, 1)]
    [InlineData(3, -14, -3, 14, -3, 7, 0, 1)]
    public void AdditionAndSubtractionHandleSharedCoprimeAndIntegerDenominators(
        long a, long b, long c, long d, long sumNum, long sumDen, long differenceNum, long differenceDen)
    {
        var left = new BigRational(a, b);
        var right = new BigRational(c, d);
        AssertCanonicalEqual(new BigRational(sumNum, sumDen), left + right);
        AssertCanonicalEqual(new BigRational(differenceNum, differenceDen), left - right);
    }

    [Fact]
    public void DefaultAndZeroRetainCanonicalSemanticsAndDivisionByZeroThrows()
    {
        BigRational zero = default;
        var value = new BigRational(-11, 17);
        foreach (var z in new[] { zero, BigRational.Zero, new BigRational(0, -123) })
        {
            AssertCanonicalEqual(value, value + z);
            AssertCanonicalEqual(value, z + value);
            AssertCanonicalEqual(value, value - z);
            AssertCanonicalEqual(-value, z - value);
            AssertCanonicalEqual(BigRational.Zero, z + z);
            AssertCanonicalEqual(BigRational.Zero, z - z);
            AssertCanonicalEqual(BigRational.Zero, z * value);
            AssertCanonicalEqual(BigRational.Zero, value * z);
            AssertCanonicalEqual(BigRational.Zero, z / value);
            Assert.Throws<DivideByZeroException>(() => value / z);
            Assert.Throws<DivideByZeroException>(() => z / z);
        }
    }

    [Fact]
    public void LargeCrossFactorsAndExactCancellationRemainCanonical()
    {
        var scale = BigInteger.One << 4096;
        var value = new BigRational(scale, scale + 1);
        var reciprocal = new BigRational(scale + 1, scale);
        AssertCanonicalEqual(BigRational.One, value * reciprocal);
        AssertCanonicalEqual(BigRational.One, value / value);
        AssertCanonicalEqual(BigRational.Zero, value - value);
        AssertCanonicalEqual(BigRational.Zero, value + (-value));
        AssertCanonicalEqual(BigRational.One, new BigRational(1, scale + 1) + value);
    }

    [Fact]
    public void DeterministicLargeOperandsMatchIndependentUnreducedReferenceFormulas()
    {
        var random = new Random(271828);
        foreach (int bits in new[] { 8, 64, 256, 1024, 4096 })
        {
            for (int sample = 0; sample < 24; sample++)
            {
                var a = SignedInteger(random, bits);
                var b = PositiveInteger(random, bits);
                var c = SignedInteger(random, bits);
                var d = PositiveInteger(random, bits);
                switch (sample % 6)
                {
                    case 0: d = b; break;
                    case 1: d = BigInteger.Abs(a); break;
                    case 2: c = b; d = BigInteger.Abs(a); break;
                    case 3: b = 1; break;
                    case 4: d = 1; break;
                }
                if (sample % 11 == 0) a = 0;
                if (sample % 13 == 0) c = 0;
                if (sample % 2 == 0) b = -b;
                if (sample % 3 == 0) d = -d;
                var left = new BigRational(a, b);
                var right = new BigRational(c, d);

                AssertCanonicalEqual(new BigRational(left.Num * right.Den + right.Num * left.Den, left.Den * right.Den), left + right);
                AssertCanonicalEqual(new BigRational(left.Num * right.Den - right.Num * left.Den, left.Den * right.Den), left - right);
                AssertCanonicalEqual(new BigRational(left.Num * right.Num, left.Den * right.Den), left * right);
                if (!right.IsZero)
                    AssertCanonicalEqual(new BigRational(left.Num * right.Den, left.Den * right.Num), left / right);
                else Assert.Throws<DivideByZeroException>(() => left / right);
            }
        }
    }

    private static BigInteger PositiveInteger(Random random, int bits)
    {
        var bytes = new byte[(bits + 7) / 8 + 1];
        random.NextBytes(bytes);
        bytes[bytes.Length - 1] = 0;
        int topBit = (bits - 1) % 8;
        bytes[bytes.Length - 2] &= (byte)((1 << (topBit + 1)) - 1);
        bytes[bytes.Length - 2] |= (byte)(1 << topBit);
        return new BigInteger(bytes);
    }

    private static BigInteger SignedInteger(Random random, int bits)
    {
        var value = PositiveInteger(random, bits);
        return random.Next(2) == 0 ? value : -value;
    }

    private static void AssertCanonicalEqual(BigRational expected, BigRational actual)
    {
        Assert.Equal(expected, actual);
        Assert.Equal(expected.Num, actual.Num);
        Assert.Equal(expected.Den, actual.Den);
        Assert.Equal(expected.GetHashCode(), actual.GetHashCode());
        Assert.Equal(expected.ToString(), actual.ToString());
        Assert.True(actual.Den.Sign > 0);
        Assert.Equal(BigInteger.One, BigInteger.GreatestCommonDivisor(BigInteger.Abs(actual.Num), actual.Den));
        if (actual.IsZero) Assert.Equal(BigInteger.One, actual.Den);
    }
}

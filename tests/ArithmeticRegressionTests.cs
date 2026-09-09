using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class ArithmeticRegressionTests
{
    [Fact]
    public void DefaultRationalIsUsableAsZero()
    {
        BigRational zero = default;
        Assert.Equal(BigInteger.One, zero.Den);
        Assert.Equal(BigRational.Zero, zero);
        Assert.Equal(BigRational.Zero.GetHashCode(), zero.GetHashCode());
        Assert.Equal("0", zero.ToString());
        Assert.Equal(new BigRational(3, 2), zero + new BigRational(3, 2));
        Assert.True(zero < BigRational.One);
        Assert.True(BigRational.IsSquare(zero, out var root));
        Assert.Equal(BigRational.Zero, root);
        Assert.Equal(new EllipticCurvePoint(0, 0), default(EllipticCurvePoint));
        Assert.Throws<DivideByZeroException>(() => new BigRational(0, 0));
    }

    [Fact]
    public void CurveEqualityAcceptsNull()
    {
        var e = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.False(e.Equals((EllipticCurveQ)null));
        Assert.False(e.Equals((object)null));
        Assert.True(e.Equals(new EllipticCurveQ(0, 0, 0, -1, 0)));
    }

    [Fact]
    public void IntegralPointsRequireBothCoordinatesToBeIntegers()
    {
        var fractional = new EllipticCurveQ(0, 0, 0, 0, new(1, 4));
        Assert.Contains(new EllipticCurvePoint(0, new(1, 2)), fractional.RationalPoints(0, 1));
        Assert.Equal(new[] { EllipticCurvePoint.Infinity }, fractional.IntegralPoints(0));
        var integral = new EllipticCurveQ(0, 0, 0, -1, 1);
        Assert.Contains(new EllipticCurvePoint(0, 1), integral.IntegralPoints(0));
    }

    [Theory]
    [InlineData(-1, 1)]
    [InlineData(int.MinValue, 1)]
    [InlineData(1, 0)]
    [InlineData(1, -1)]
    public void RationalPointSearchRejectsInvalidBounds(int numerator, int denominator)
    {
        var e = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.RationalPoints(numerator, denominator).ToArray());
    }

    [Fact]
    public void SixthRootsDoNotStartNewtonIterationBelowTheRoot()
    {
        var power = BigInteger.Pow(5, 6);
        Assert.True(InternalMath.TryIntegerKthRoot(power, 6, out var root));
        Assert.Equal(new BigInteger(5), root);
        Assert.Equal(new BigInteger(4), InternalMath.IntegerKthRoot(power - 1, 6));
        Assert.Equal(new BigInteger(5), InternalMath.IntegerKthRoot(power + 1, 6));
        // For j=0, isomorphism needs a sixth root rather than a fourth root.
        var e = new EllipticCurveQ(0, 0, 0, 0, new(power));
        var other = new EllipticCurveQ(0, 0, 0, 0, 1);
        Assert.True(e.IsIsomorphic(other, out var u));
        Assert.Equal(new BigRational(5), u);
    }

    [Fact]
    public void IsomorphismScaleRelatesTheActualInputModels()
    {
        foreach (var e in new[] { new EllipticCurveQ(0, 0, 1, -1, 0), new EllipticCurveQ(0, 0, 0, -1, 0),
            new EllipticCurveQ(0, 0, 0, 0, 1) })
        {
            var other = NativeArithmeticTests.ChangeCoordinates(e, new(5, 2), new(1, 3), new(-2, 3), new(3, 7));
            Assert.True(e.IsIsomorphic(other, out var u));
            Assert.Equal(new BigRational(5, 2), u);
            Assert.Equal(e.C4, BigRational.Pow(u, 4) * other.C4);
            Assert.Equal(e.C6, BigRational.Pow(u, 6) * other.C6);
            Assert.Equal(e.Discriminant, BigRational.Pow(u, 12) * other.Discriminant);
            Assert.True(other.IsIsomorphic(e, out var inverse));
            Assert.Equal(BigRational.One, u * inverse);
        }
        Assert.False(new EllipticCurveQ(0, 0, 0, 0, 1).IsIsomorphic(new(0, 0, 0, 0, 2)));
    }

    [Fact]
    public void ReturnedTorsionCollectionCannotCorruptTheCurveCache()
    {
        var e = new EllipticCurveQ(0, 0, 0, -1, 0);
        var points = e.TorsionPoints;
        if (points is ICollection<EllipticCurvePoint> collection)
        {
            if (collection.IsReadOnly) Assert.Throws<NotSupportedException>(() => collection.Clear());
            else collection.Clear(); // A mutable snapshot may be changed, but never the internal cache.
        }
        Assert.Equal(4, e.TorsionPoints.Count());
        Assert.Equal("Z/2Z x Z/2Z", e.TorsionStructure);
        Assert.Equal(2, e.TorsionOrder(new(0, 0)));
    }

    [Fact]
    public void TorsionDivisorHelperDoesNotTreatAPseudoprimeAsPrime()
    {
        var n = BigInteger.Parse("341550071728321");
        var factors = InternalMath.FactorAbs(n);
        Assert.Equal(2, factors.Count);
        Assert.Equal(1, factors[10670053]);
        Assert.Equal(1, factors[32010157]);
        Assert.Contains(new BigInteger(10670053), InternalMath.EnumerateDivisorsAbs(n));
    }

    [Theory]
    [InlineData(-1, 0, 4, "Z/2Z x Z/2Z")]
    [InlineData(0, 1, 6, "Z/6Z")]
    [InlineData(0, 4, 3, "Z/3Z")]
    [InlineData(-1, 1, 1, "Z/1Z")]
    public void TorsionEnumerationPreservesKnownGroupsAndOriginalCoordinates(
        int a4, int a6, int count, string structure)
    {
        var original = new EllipticCurveQ(0, 0, 0, a4, a6);
        var shifted = NativeArithmeticTests.ChangeCoordinates(original, 1, 0, new(1, 2), new(1, 2));
        foreach (var e in new[] { original, shifted })
        {
            var points = e.TorsionPoints.ToArray();
            Assert.Equal(count, points.Length);
            Assert.Equal(structure, e.TorsionStructure);
            Assert.All(points, p =>
            {
                Assert.True(e.IsOnCurve(p));
                var order = e.TorsionOrder(p);
                Assert.NotNull(order);
                Assert.Equal(0, count % order.Value);
                Assert.True(e.Multiply(p, order.Value).IsInfinity);
            });
        }
    }
}

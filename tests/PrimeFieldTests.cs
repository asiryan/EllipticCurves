using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class PrimeFieldTests
{
    [Theory]
    [InlineData(2), InlineData(3), InlineData(5), InlineData(7), InlineData(13), InlineData(17), InlineData(29)]
    public void EnumerationAndGroupLawAgreeWithExhaustiveEquations(int prime)
    {
        var random = new Random(101 + prime);
        for (int sample = 0; sample < 12; sample++)
        {
            var a = Enumerable.Range(0, 5).Select(_ => random.Next(prime)).ToArray();
            var rational = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
            if (rational.Discriminant.Num % prime == 0) continue;
            var e = new EllipticCurveFp(prime, a[0], a[1], a[2], a[3], a[4]);
            var expected = new HashSet<(int, int)>();
            for (int x = 0; x < prime; x++) for (int y = 0; y < prime; y++)
                if ((y * y + a[0] * x * y + a[2] * y - x * x * x - a[1] * x * x - a[3] * x - a[4]) % prime == 0) expected.Add((x, y));
            var points = e.Points().ToArray(); Assert.Equal(expected.Count + 1, points.Length);
            Assert.Equal(new BigInteger(points.Length), e.CountPoints()); Assert.Equal(points.Length, points.Distinct().Count());
            Assert.True(expected.SetEquals(points.Where(p => !p.IsInfinity).Select(p => ((int)p.X, (int)p.Y))));
            Assert.Equal(NativeNumberTheory.Mod(rational.JInvariant.Num * BigInteger.ModPow(rational.JInvariant.Den, prime - 2, prime), prime), e.JInvariant);
            foreach (var p in points)
            {
                Assert.True(e.Multiply(p, points.Length).IsInfinity); Assert.True(e.Add(p, e.Negate(p)).IsInfinity);
                var order = e.GetPointOrder(p); Assert.True(e.Multiply(p, order).IsInfinity);
                for (int n = 1; n < order; n++) Assert.False(e.Multiply(p, n).IsInfinity);
                Assert.Equal(e.Negate(e.Double(p)), e.Multiply(p, -2));
                var q = points[random.Next(points.Length)]; var r = points[random.Next(points.Length)];
                Assert.True(e.IsOnCurve(e.Add(p, q))); Assert.Equal(e.Add(p, q), e.Add(q, p));
                Assert.Equal(e.Add(e.Add(p, q), r), e.Add(p, e.Add(q, r)));
            }
        }
    }

    [Fact]
    public void ReductionRespectsTheGroupLawAndModelChanges()
    {
        var original = new EllipticCurveQ(0, 0, 1, -1, 0); var generator = new EllipticCurvePoint(0, 0);
        var change = original.ChangeModel(new BigRational(2, 3), 5, -1, 7); var e = change.Target;
        foreach (int prime in new[] { 2, 3, 5, 7, 11 })
        {
            var finite = e.ReduceModuloPrime(prime); var p = e.ReducePointModuloPrime(change.Map(generator), prime);
            Assert.Equal(e.CountPoints(prime), finite.CountPoints());
            for (int n = -8; n <= 8; n++)
                Assert.Equal(finite.Multiply(p, n), e.ReducePointModuloPrime(change.Map(original.Multiply(generator, n)), prime));
        }
        Assert.True(original.ReducePointModuloPrime(original.Multiply(generator, 5), 2).IsInfinity);
    }

    [Fact]
    public void InvalidFieldsPointsAndExcessiveSearchesAreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new EllipticCurveFp(9, 0, 0, 0, -1, 0));
        Assert.Throws<ArgumentException>(() => new EllipticCurveFp(5, 0, 0, 0, 0, 0));
        var e = new EllipticCurveFp(7, 0, 0, 0, -1, 0); var f = new EllipticCurveFp(11, 0, 0, 0, -1, 0);
        Assert.Throws<ArgumentException>(() => e.Add(e.CreatePoint(0, 0), f.CreatePoint(0, 0)));
        Assert.Throws<ArgumentException>(() => e.CreatePoint(0, 1)); Assert.False(e.IsOnCurve(default));
        Assert.Equal(e.CreatePoint(0, 0), e.CreatePoint(7, -7));
        Assert.Throws<ArithmeticException>(() => e.CountPoints(6)); Assert.Throws<ArithmeticException>(() => e.Points(6).ToArray());
        Assert.Throws<OperationCanceledException>(() => e.CountPoints(cancellationToken: new CancellationToken(true)));
        Assert.Throws<ArgumentException>(() => new EllipticCurveQ(0, 0, 1, -1, 0).ReduceModuloPrime(37));
    }
}

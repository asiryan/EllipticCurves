using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class PrimePointSearchTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    public void AllSmallGeneralEquationsMatchDirectCoordinateSearch(int prime)
    {
        int equations = (int)Math.Pow(prime, 5);
        for (int encoded = 0; encoded < equations; encoded++)
        {
            var a = Decode(encoded, prime);
            var expected = BrutePoints(prime, a);
            var rational = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
            // Also compare singular equations: their counts give the bad-prime
            // Fourier coefficient and must not be subjected to a smoothness test.
            Assert.Equal(prime - expected.Count, AnalyticCoefficients.Trace(rational, prime, default));
            if (prime != 2)
                Assert.Equal(new BigInteger(expected.Count + 1), PrimeFieldPointSearch.CountWithoutTable(
                    prime, a[0], a[1], a[2], a[3], a[4], default));
            if (rational.Discriminant.Num % prime == 0) continue;
            var curve = new EllipticCurveFp(prime, a[0], a[1], a[2], a[3], a[4]);
            var points = curve.Points(prime).ToArray();
            Assert.True(points[0].IsInfinity);
            Assert.Equal(expected.Count + 1, points.Length);
            Assert.Equal(points.Length, points.Distinct().Count());
            Assert.True(expected.SetEquals(points.Skip(1).Select(p => ((int)p.X, (int)p.Y))));
            Assert.Equal(new BigInteger(points.Length), curve.CountPoints(prime));
        }
    }

    [Theory]
    [InlineData(7)]
    [InlineData(17)]
    [InlineData(101)]
    public void SeededGeneralModelsAndCoefficientRepresentativesMatchDirectSearch(int prime)
    {
        var random = new Random(67211 + prime);
        for (int sample = 0; sample < 24; sample++)
        {
            var a = Enumerable.Range(0, 5).Select(_ => random.Next(-2 * prime, 2 * prime)).ToArray();
            var expected = BrutePoints(prime, a);
            var rational = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
            Assert.Equal(new BigInteger(expected.Count + 1), PrimeFieldPointSearch.CountWithoutTable(
                prime, a[0], a[1], a[2], a[3], a[4], default));
            if (rational.Discriminant.Num % prime == 0) continue;
            var curve = new EllipticCurveFp(prime, a[0], a[1], a[2], a[3], a[4]);
            var actual = curve.Points(prime).ToArray();
            Assert.Equal(expected.Count + 1, actual.Length);
            Assert.True(expected.SetEquals(actual.Skip(1).Select(p => ((int)p.X, (int)p.Y))));
            Assert.Equal(new BigInteger(actual.Length), curve.CountPoints(prime));
        }
    }

    [Fact]
    public void CachedCountsAndPointOrdersStillEnforceLimitsAndCancellation()
    {
        var curve = new EllipticCurveFp(101, 0, 0, 0, -1, 0);
        var points = curve.Points(101).ToArray();
        Assert.Equal(new BigInteger(points.Length), curve.CountPoints(101));
        foreach (var p in points)
        {
            var expected = 1;
            var multiple = p;
            while (!multiple.IsInfinity) { multiple = curve.Add(multiple, p); expected++; }
            Assert.Equal(new BigInteger(expected), curve.GetPointOrder(p, 101));
        }
        var affine = points.First(p => !p.IsInfinity);
        Assert.Throws<ArithmeticException>(() => curve.CountPoints(100));
        Assert.Throws<ArithmeticException>(() => curve.GetPointOrder(affine, 100));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.CountPoints(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetPointOrder(EllipticCurvePointFp.Infinity, -1));
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(101, new CancellationToken(true)));
        Assert.Throws<OperationCanceledException>(() => curve.GetPointOrder(affine, 101, new CancellationToken(true)));
        Assert.Throws<ArgumentException>(() => curve.GetPointOrder(default, 101));
        Assert.Equal(BigInteger.One, curve.GetPointOrder(EllipticCurvePointFp.Infinity, 0));
    }

    [Fact]
    public async Task ConcurrentCountsPublishOnlyTheCompletedValue()
    {
        var curve = new EllipticCurveFp(100003, 1, 2, 3, 4, 5);
        var counts = await Task.WhenAll(Enumerable.Range(0, 6).Select(_ => Task.Run(() => curve.CountPoints(100003))));
        Assert.All(counts, count => Assert.Equal(counts[0], count));
        Assert.Equal(counts[0], new BigInteger(curve.Points(100003).Count()));
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(100003, new CancellationToken(true)));
    }

    [Fact]
    public void LargePrimeUsesBoundedStorageAndKeepsEnumerationLazy()
    {
        const long prime = 2147483647; // A proven Mersenne prime; larger than every table cap.
        var curve = new EllipticCurveFp(prime, 0, 0, 0, -1, 0);
        using var enumerator = curve.Points(prime).GetEnumerator();
        Assert.True(enumerator.MoveNext());
        Assert.True(enumerator.Current.IsInfinity);
        var sample = curve.Points(prime).Take(25).ToArray();
        Assert.Equal(sample.Length, sample.Distinct().Count());
        Assert.All(sample, p => Assert.True(curve.IsOnCurve(p)));
        using var timeout = new CancellationTokenSource(TimeSpan.FromMilliseconds(20));
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(prime, timeout.Token));
        Assert.Throws<ArithmeticException>(() => curve.CountPoints(prime - 1));
    }

    [Fact]
    public void CancellationBetweenEnumerationResultsAndRetryAreSupported()
    {
        var curve = new EllipticCurveFp(1009, 1, 2, 3, 4, 5);
        using var cancellation = new CancellationTokenSource();
        using var enumerator = curve.Points(1009, cancellation.Token).GetEnumerator();
        Assert.True(enumerator.MoveNext());
        Assert.True(enumerator.Current.IsInfinity);
        Assert.True(enumerator.MoveNext());
        cancellation.Cancel();
        Assert.Throws<OperationCanceledException>(() => enumerator.MoveNext());
        Assert.Equal(new BigInteger(curve.Points(1009).Count()), curve.CountPoints(1009));
    }

    private static int[] Decode(int encoded, int prime)
    {
        var a = new int[5];
        for (int i = 0; i < a.Length; i++) { a[i] = encoded % prime; encoded /= prime; }
        return a;
    }

    private static HashSet<(int x, int y)> BrutePoints(int prime, int[] a)
    {
        var points = new HashSet<(int x, int y)>();
        for (long x = 0; x < prime; x++) for (long y = 0; y < prime; y++)
            if ((y * y + a[0] * x * y + a[2] * y - x * x * x - a[1] * x * x - a[3] * x - a[4]) % prime == 0)
                points.Add(((int)x, (int)y));
        return points;
    }
}

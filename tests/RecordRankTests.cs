using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;
using Xunit;
using Xunit.Abstractions;

namespace EllipticCurves.Tests;

public class RecordRankTests(ITestOutputHelper output)
{
    private const int PrimeBound = 1009;

    [Fact]
    public void PublishedRecordHas31IndependentRationalPoints()
    {
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var timer = Stopwatch.StartNew();
        var (curve, points) = LoadRecord();
        Assert.Equal(31, points.Length);
        Assert.Equal(31, points.Distinct().Count());
        Assert.False(curve.IsSingular);
        Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));

        // The file's rank, torsion, conductor and factorization are not proof inputs.
        int torsionWitness = ProveNoRationalTwoTorsion(curve, timeout.Token);
        output.WriteLine("All 31 points satisfy the curve equation exactly.");
        output.WriteLine($"No rational 2-torsion: good prime {torsionWitness}, cubic has no roots.");
        foreach (int bound in new[] { 101, 251, 503, PrimeBound })
        {
            var certificate = Certify(curve, points, bound, timeout.Token);
            output.WriteLine($"Good odd primes <= {bound}: binary image dimension = {certificate.ImageDimension}; rank >= {certificate.LowerBound}.");
            if (bound == PrimeBound)
            {
                Assert.Equal(31, certificate.ImageDimension);
                Assert.Equal(31, certificate.LowerBound);
            }
        }
        output.WriteLine($"Elapsed verification time: {timer.Elapsed.TotalMilliseconds:F1} ms.");
        output.WriteLine("Unconditional lower bound only: the exact rank and a full Mordell-Weil basis are not certified.");
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void DependentReplacementCannotCertify31(bool useSum)
    {
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var (curve, points) = LoadRecord();
        ProveNoRationalTwoTorsion(curve, timeout.Token);
        points[30] = useSum ? curve.Add(points[0], points[1]) : points[0];
        var certificate = Certify(curve, points, PrimeBound, timeout.Token);
        Assert.Equal(30, certificate.ImageDimension);
        Assert.Equal(30, certificate.LowerBound);
    }

    [Fact]
    public void CorruptedCoordinateIsRejected()
    {
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var (curve, points) = LoadRecord();
        ProveNoRationalTwoTorsion(curve, timeout.Token);
        var corrupted = new EllipticCurvePoint(points[0].X, points[0].Y + 1);
        Assert.False(curve.IsOnCurve(corrupted));
        var certificate = Certify(curve, Array.Empty<EllipticCurvePoint>(), PrimeBound, timeout.Token);
        Assert.Throws<InvalidOperationException>(() => certificate.Add(corrupted));
    }

    private static RationalPointRank Certify(EllipticCurveQ curve, IEnumerable<EllipticCurvePoint> points,
        int primeBound, CancellationToken token)
    {
        var budget = new DescentBudget(new RankComputationOptions
        {
            ReductionPrimeBound = primeBound,
            SearchBound = 0,
            MaxDescentWork = 0,
            MaxPointSearchWork = 0
        }, token);
        // Each caller proves torsionDimension = 0 first. The input model is integral;
        // good reduction suffices, so no minimalization or factorization is needed.
        var certificate = new RationalPointRank(curve, 0, budget);
        foreach (var point in points) certificate.Add(point);
        Assert.Equal(0L, budget.Work);
        Assert.Equal(0L, budget.PointWork);
        return certificate;
    }

    private static int ProveNoRationalTwoTorsion(EllipticCurveQ curve, CancellationToken token)
    {
        Assert.All(new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 },
            coefficient => Assert.Equal(BigInteger.One, coefficient.Den));
        // X=4x, Y=8y+4a1*x+4a3: Y^2=X^3+b2*X^2+8*b4*X+16*b6.
        // At a good odd prime, rational 2-torsion injects into E(F_p)[2].
        // An irreducible cubic mod p therefore certifies E(Q)[2] = 0.
        for (int p = 3; p <= PrimeBound; p += 2)
        {
            token.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(p, token) || curve.Discriminant.Num % p == 0) continue;
            long a = (long)NativeNumberTheory.Mod(curve.B2.Num, p);
            long b = (long)NativeNumberTheory.Mod(8 * curve.B4.Num, p);
            long c = (long)NativeNumberTheory.Mod(16 * curve.B6.Num, p);
            bool hasRoot = false;
            for (long x = 0; x < p; x++)
                if ((((x + a) * x + b) * x + c) % p == 0) { hasRoot = true; break; }
            if (!hasRoot) return p;
        }
        throw new InvalidOperationException("No certificate excluding rational 2-torsion was found.");
    }

    private static (EllipticCurveQ curve, EllipticCurvePoint[] points) LoadRecord()
    {
        using var document = JsonDocument.Parse(File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "Fixtures", "icarm-302.json")));
        var root = document.RootElement;
        var a = root.GetProperty("ainvs").EnumerateArray().Select(value => Parse(value.GetString()!)).ToArray();
        var curve = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        var points = root.GetProperty("points").EnumerateArray().Select(point =>
            new EllipticCurvePoint(Parse(point[0].GetString()!), Parse(point[1].GetString()!))).ToArray();
        return (curve, points);
    }

    private static BigRational Parse(string value)
    {
        var parts = value.Split('/');
        if (parts.Length is < 1 or > 2) throw new FormatException("Expected an integer or a rational fraction.");
        return new BigRational(BigInteger.Parse(parts[0], CultureInfo.InvariantCulture),
            parts.Length == 2 ? BigInteger.Parse(parts[1], CultureInfo.InvariantCulture) : BigInteger.One);
    }
}

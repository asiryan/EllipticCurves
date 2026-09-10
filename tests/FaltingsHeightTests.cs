using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class FaltingsHeightTests
{
    public static IEnumerable<object[]> Rows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "faltings-heights.csv")).Select(s => new object[] { s });

    [Theory, MemberData(nameof(Rows))]
    public void HeightsEncloseIndependentPariReferences(string row)
    {
        var v = row.Split(',');
        var a = v.Take(5).Select(ExtendedReferenceTests.Parse).ToArray();
        var e = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        ExtendedReferenceTests.ReferenceIn(e.FaltingsHeight(), v[5]);
        ExtendedReferenceTests.ReferenceIn(e.StableFaltingsHeight(), v[6]);
    }

    [Theory]
    [InlineData("37a1"), InlineData("11a1"), InlineData("48a3"), InlineData("389a1"), InlineData("5077a1")]
    public void MinimalModelNormalizationMatchesLmfdb(string label)
    {
        var stored = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture(label));
        var e = stored.GlobalMinimalModel;
        var tolerance = new BigRational(1, BigInteger.Pow(10, 20));
        var h = e.FaltingsHeight(); var stable = e.StableFaltingsHeight();
        var reference = stored.FaltingsHeight.AsRational();
        Assert.True(h.LowerBound <= reference + tolerance && h.UpperBound >= reference - tolerance);
        reference = stored.StableFaltingsHeight.AsRational();
        Assert.True(stable.LowerBound <= reference + tolerance && stable.UpperBound >= reference - tolerance);
    }

    [Theory]
    [InlineData(12, 256), InlineData(30, 512), InlineData(60, 1024)]
    public void PrecisionAndRationalModelChangesPreserveBothHeights(int digits, int bits)
    {
        var row = Rows().Select(r => ((string)r[0]).Split(',')).First(v => v[0] == "0" && v[1] == "0" && v[2] == "0" && v[3] == "-1" && v[4] == "0");
        var e = new EllipticCurveQ(0, 0, 0, -1, 0);
        var changed = e.ChangeModel(new BigRational(-7, 3), 17, -5, 9).Target;
        var options = new RealComputationOptions { DecimalDigits = digits, PrecisionBits = bits };
        ExtendedReferenceTests.ReferenceIn(changed.FaltingsHeight(options), row[5], digits);
        ExtendedReferenceTests.ReferenceIn(changed.StableFaltingsHeight(options), row[6], digits);
    }

    [Fact]
    public void StableHeightIsTwistInvariantAndSemistableHeightsAgree()
    {
        foreach (var e in new[] { new EllipticCurveQ(0, 0, 1, -1, 0), new EllipticCurveQ(0, 0, 0, -1, 0), new EllipticCurveQ(0, 0, 0, 0, 1) })
        {
            var h = e.StableFaltingsHeight();
            foreach (int d in new[] { -7, -1, 2, 3, 5 })
            {
                var twist = e.QuadraticTwist(d).StableFaltingsHeight();
                Assert.True(h.LowerBound <= twist.UpperBound && h.UpperBound >= twist.LowerBound);
            }
        }
        var semistable = new EllipticCurveQ(0, 0, 1, -1, 0);
        var ordinary = semistable.FaltingsHeight(); var stable = semistable.StableFaltingsHeight();
        Assert.True(ordinary.LowerBound <= stable.UpperBound && ordinary.UpperBound >= stable.LowerBound);
        var additive = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.True(additive.StableFaltingsHeight().UpperBound < additive.FaltingsHeight().LowerBound);
    }

    [Fact]
    public void InvalidInputsCancellationAndPrecisionFailureAreExplicit()
    {
        var singular = new EllipticCurveQ(0, 0, 0, 0, 0);
        Assert.Throws<InvalidOperationException>(() => singular.FaltingsHeight());
        Assert.Throws<InvalidOperationException>(() => singular.StableFaltingsHeight());
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        Assert.Throws<OperationCanceledException>(() => e.FaltingsHeight(cancellationToken: new CancellationToken(true)));
        Assert.Throws<OperationCanceledException>(() => e.StableFaltingsHeight(cancellationToken: new CancellationToken(true)));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.FaltingsHeight(new RealComputationOptions { DecimalDigits = 0 }));
        Assert.Throws<ArithmeticException>(() => e.StableFaltingsHeight(new RealComputationOptions { MaxRootWork = 0 }));
        Assert.Throws<ArithmeticException>(() => e.FaltingsHeight(new RealComputationOptions { DecimalDigits = 100, PrecisionBits = 64 }));
    }
}

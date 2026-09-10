using System.Globalization;
using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class ExtendedReferenceTests
{
    public static IEnumerable<object[]> LocalRows() => Rows("local-data.csv");
    public static IEnumerable<object[]> HeightRows() => Rows("heights.csv");
    public static IEnumerable<object[]> PeriodRows() => Rows("periods.csv");
    private static IEnumerable<object[]> Rows(string file) => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", file)).Distinct().Select(x => new object[] { x });
    private static EllipticCurveQ Curve(string[] v) => new(Parse(v[0]), Parse(v[1]), Parse(v[2]), Parse(v[3]), Parse(v[4]));
    internal static BigRational Parse(string s)
    {
        if (s.Contains('/')) { var q = s.Split('/'); return new BigRational(BigInteger.Parse(q[0]), BigInteger.Parse(q[1])); }
        var parts = s.Split('e', 'E'); int exponent = parts.Length == 2 ? int.Parse(parts[1], CultureInfo.InvariantCulture) : 0;
        int dot = parts[0].IndexOf('.'); if (dot >= 0) exponent -= parts[0].Length - dot - 1;
        var n = BigInteger.Parse(parts[0].Replace(".", ""), CultureInfo.InvariantCulture);
        return exponent >= 0 ? new BigRational(n * BigInteger.Pow(10, exponent)) : new BigRational(n, BigInteger.Pow(10, -exponent));
    }
    internal static void ReferenceIn(RealEnclosure enclosure, string reference, int digits = 12)
    {
        var value = Parse(reference); var referenceError = new BigRational(1, BigInteger.Pow(10, 65));
        Assert.True(enclosure.LowerBound <= value + referenceError && enclosure.UpperBound >= value - referenceError,
            $"Reference {reference} is outside [{enclosure.LowerBound}, {enclosure.UpperBound}].");
        Assert.True(enclosure.Width <= new BigRational(1, BigInteger.Pow(10, digits)));
    }
    [Theory, MemberData(nameof(LocalRows))]
    public void LocalInvariantsMatchPari(string row)
    {
        var v = row.Split(','); var e = Curve(v); var d = e.GetLocalData(BigInteger.Parse(v[5]));
        int code = int.Parse(v[8]), k = Math.Abs(code);
        string symbol = (k == 1 ? "I0" : k == 2 ? "II" : k == 3 ? "III" : k == 4 ? "IV" : "I" + (k - 4)) + (code < 0 ? "*" : "");
        Assert.Equal(int.Parse(v[6]), d.DiscriminantValuation);
        Assert.Equal(int.Parse(v[7]), d.ConductorValuation);
        Assert.Equal(symbol, d.KodairaSymbol);
        Assert.Equal(int.Parse(v[9]), d.TamagawaNumber);
        Assert.Equal(int.Parse(v[10]), d.RootNumber);
        Assert.Equal(code >= 5 ? (d.RootNumber == -1 ? ReductionType.SplitMultiplicative : ReductionType.NonSplitMultiplicative) : ReductionType.Additive, d.ReductionType);
    }
    [Theory, MemberData(nameof(HeightRows))]
    public void CanonicalHeightsEnclosePariValues(string row)
    {
        var v = row.Split(','); var e = Curve(v); var p = new EllipticCurvePoint(Parse(v[5]), Parse(v[6]));
        ReferenceIn(e.CanonicalHeight(p), v[7]);
    }
    [Theory, MemberData(nameof(PeriodRows))]
    public void PeriodsEnclosePariValues(string row)
    {
        var v = row.Split(','); var result = Curve(v).GetPeriods();
        ReferenceIn(result.PrimitiveRealPeriod, v[5]);
        ReferenceIn(result.SecondPeriodImaginaryPart, v[6]);
        ReferenceIn(result.RealPeriod, v[7]);
        ReferenceIn(result.Area, v[8]);
    }
}

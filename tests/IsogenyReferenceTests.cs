using EllipticCurves;
using Xunit;
using static EllipticCurves.Tests.ExtendedReferenceTests;

namespace EllipticCurves.Tests;

public class IsogenyReferenceTests
{
    public static IEnumerable<object[]> Rows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "isogenies.csv")).Select(row => new object[] { row });
    [Theory, MemberData(nameof(Rows))]
    public void VeluQuotientInvariantsMatchPari(string row)
    {
        var v = row.Split(','); var e = new EllipticCurveQ(Parse(v[0]), Parse(v[1]), Parse(v[2]), Parse(v[3]), Parse(v[4]));
        var isogeny = e.CreateIsogeny(new[] { new EllipticCurvePoint(Parse(v[5]), Parse(v[6])) });
        Assert.Equal(int.Parse(v[7]), isogeny.Degree); Assert.Equal(Parse(v[8]), isogeny.Target.C4); Assert.Equal(Parse(v[9]), isogeny.Target.C6);
    }
}

using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

// This regression intentionally uses many CPU workers. Keep it separate from
// the suite's concurrent tests with short cancellation budgets.
[CollectionDefinition("Parallel factorization", DisableParallelization = true)]
public class ParallelFactorizationCollection { }

[Collection("Parallel factorization")]
public class ParallelFactorizationTests
{
    [Fact]
    public void SeventyFiveDigitResidualConductorMatchesProvenPariResult()
    {
        var curve = new EllipticCurveQ(1, 0, 0,
            new BigRational(BigInteger.Parse("-20820207864197471248300179976626")),
            new BigRational(BigInteger.Parse("36732936589138673862895758597955508398047757956")));
        using var timeout = new CancellationTokenSource(TimeSpan.FromMinutes(5));
        var conductor = curve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = 12 }, out var factors, timeout.Token);
        Assert.Equal(BigInteger.Parse("21785392458764315483988614758901932764423833726404496768609056369259382575196330"), conductor);
        Assert.Equal(new[] { "2", "3", "5", "7", "13", "19", "86044359449701746413144681",
            "4881201479084572277272794646126863424004819546339" }, factors.Keys.Select(p => p.ToString()));
        Assert.All(factors.Values, exponent => Assert.Equal(1, exponent));
        Assert.Equal(conductor, factors.Aggregate(BigInteger.One, (n, factor) => n * BigInteger.Pow(factor.Key, factor.Value)));
    }
}

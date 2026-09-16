using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class NativeNumberTheoryTests
{
    public static IEnumerable<object[]> LargeFactorizations() => File.ReadLines(
        Path.Combine(AppContext.BaseDirectory, "Fixtures", "factorization.csv"))
        .Skip(1).Where(line => !string.IsNullOrWhiteSpace(line))
        .Select(line => line.Split(',').Cast<object>().ToArray());

    [Theory]
    [MemberData(nameof(LargeFactorizations))]
    [InlineData("10001234603", "1000007654329")]
    [InlineData("1000002469163", "100000015308653")]
    [InlineData("100000003703729", "10000000022962973")]
    [InlineData("10000000004938337", "1000000000030617347")]
    [InlineData("1000000000006172849", "100000000000038271633")]
    public void FactorsSemiprimesWithCertifiedFactors(string left, string right)
    {
        var p = BigInteger.Parse(left);
        var q = BigInteger.Parse(right);
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var factors = NativeNumberTheory.Factor(-p * q, timeout.Token);
        Assert.Equal(2, factors.Count);
        Assert.Equal(1, factors[p]);
        Assert.Equal(1, factors[q]);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    public void ExtractsLargePerfectPowers(int exponent)
    {
        var p = BigInteger.Pow(2, 89) - 1;
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(5));
        var factors = NativeNumberTheory.Factor(-72 * BigInteger.Pow(p, exponent), timeout.Token);
        Assert.Equal(3, factors.Count);
        Assert.Equal(3, factors[2]);
        Assert.Equal(2, factors[3]);
        Assert.Equal(exponent, factors[p]);
    }

    [Fact]
    public void PreservesRepeatedFactorsAfterSplitting()
    {
        BigInteger p = 1000003, q = 1000033;
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(10));
        var factors = NativeNumberTheory.Factor(BigInteger.Pow(p, 3) * q * q, timeout.Token);
        Assert.Equal(2, factors.Count);
        Assert.Equal(3, factors[p]);
        Assert.Equal(2, factors[q]);
    }

    [Theory]
    [InlineData("100000003703729", "10000000022962973", 1)]
    [InlineData("1575838430456954508271967", "81274068710384465721193186106423", 1)]
    [InlineData("1575838430456954508271967", "81274068710384465721193186106423", 4)]
    [InlineData("3093889989073154286827653", "5930788502963551593663361", 1)]
    [InlineData("3093889989073154286827653", "5930788502963551593663361", 4)]
    [InlineData("1575838430456954508271967", "81274068710384465721193186106423", 8)]
    [InlineData("3093889989073154286827653", "5930788502963551593663361", 24)]
    public void QuadraticSieveSplitsWithoutKnownFactors(string left, string right, int workers)
    {
        var p = BigInteger.Parse(left);
        var q = BigInteger.Parse(right);
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var divisor = NativeQuadraticSieve.FindDivisor(p * q, timeout.Token, workers);
        Assert.True(divisor == p || divisor == q);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(4)]
    [InlineData(24)]
    public void FactorizationHonorsCancellationAndUnitConventions(int workers)
    {
        Assert.Empty(NativeNumberTheory.Factor(1, default));
        Assert.Empty(NativeNumberTheory.Factor(-1, default));
        Assert.Throws<ArgumentException>(() => NativeNumberTheory.Factor(0, default));
        Assert.Throws<OperationCanceledException>(() => NativeNumberTheory.Factor(1, new CancellationToken(true)));
        var n = BigInteger.Parse("128074800873422933261289680790071262754898704850689544041");
        using var timeout = new CancellationTokenSource(TimeSpan.FromMilliseconds(50));
        Assert.Throws<OperationCanceledException>(() => NativeQuadraticSieve.FindDivisor(n, timeout.Token, workers));
    }

    [Fact]
    public void WorkerLimitsAllowAvailableCpusWithoutOversubscription()
    {
        Assert.Equal(24, NativeQuadraticSieve.WorkerCount(75, 0, 24));
        Assert.Equal(12, NativeQuadraticSieve.WorkerCount(75, 12, 24));
        Assert.Equal(8, NativeQuadraticSieve.WorkerCount(57, 8, 24));
        Assert.Equal(4, NativeQuadraticSieve.WorkerCount(57, 0, 24));
        Assert.Equal(4, NativeQuadraticSieve.WorkerCount(75, int.MaxValue, 4));
        Assert.Equal(1, NativeQuadraticSieve.WorkerCount(75, 0, 1));
        Assert.Equal(1, NativeQuadraticSieve.WorkerCount(75, 1, 24));
        Assert.Equal(1, NativeQuadraticSieve.WorkerCount(30, 24, 24));
        Assert.Throws<ArgumentOutOfRangeException>(() => NativeQuadraticSieve.WorkerCount(75, -1, 24));
        Assert.Throws<ArgumentOutOfRangeException>(() => NativeNumberTheory.Factor(1, default, -1));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(8)]
    public void FactorizationWithExplicitWorkerLimitsCertifiesLargeFactors(int workers)
    {
        var p = BigInteger.Parse("1575838430456954508271967");
        var q = BigInteger.Parse("81274068710384465721193186106423");
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(30));
        var factors = NativeNumberTheory.Factor(p * q, timeout.Token, workers);
        Assert.Equal(2, factors.Count);
        Assert.Equal(1, factors[p]);
        Assert.Equal(1, factors[q]);
    }

    [Fact]
    public void EcmFirstStageFindsASmoothOrderFactor()
    {
        var n = 1051 * BigInteger.Parse("1000000000000000000000000000007");
        Assert.Equal(new BigInteger(1051), NativeEcmFactorization.FindDivisor(n, 1, 5, 5, default));
    }

    [Theory]
    [InlineData(1009, 5, 100)]
    [InlineData(10009, 200, 2000)]
    [InlineData(10061, 200, 2000)]
    public void EcmSecondStageFindsFactorsMissedByTheFirst(int prime, int b1, int b2)
    {
        var n = prime * BigInteger.Parse("1000000000000000000000000000007");
        Assert.Equal(BigInteger.One, NativeEcmFactorization.FindDivisor(n, 1, b1, b1, default));
        Assert.Equal(new BigInteger(prime), NativeEcmFactorization.FindDivisor(n, 1, b1, b2, default));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(101)]
    [InlineData(1000003)]
    public void EcmDoesNotReturnTheWholeInput(int prime)
        => Assert.Equal(BigInteger.One, NativeEcmFactorization.FindDivisor(prime, 8, 50, 500, default));

    [Fact]
    public void EcmHonorsCancellation()
    {
        var n = BigInteger.Pow(2, 127) - 1;
        Assert.Throws<OperationCanceledException>(() =>
            NativeEcmFactorization.FindDivisor(n, 16, 500, 5000, new CancellationToken(true)));
        using var timeout = new CancellationTokenSource(TimeSpan.FromMilliseconds(20));
        Assert.Throws<OperationCanceledException>(() =>
            NativeEcmFactorization.FindDivisor(n, 1000, 10000, 100000, timeout.Token));
    }

    [Fact]
    public void DoesNotAcceptTheOldMillerRabinPseudoprime()
    {
        var factors = NativeNumberTheory.Factor(BigInteger.Parse("341550071728321"), default);
        Assert.Equal(2, factors.Count);
        Assert.Equal(1, factors[10670053]);
        Assert.Equal(1, factors[32010157]);
    }

    [Fact]
    public void ProvesPrimalityBeyondUnsigned64Bits()
    {
        var prime = BigInteger.Pow(2, 89) - 1;
        var factors = NativeNumberTheory.Factor(-prime * 72, default);
        Assert.Equal(3, factors.Count);
        Assert.Equal(3, factors[2]);
        Assert.Equal(2, factors[3]);
        Assert.Equal(1, factors[prime]);
    }
}

using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class NativeNumberTheoryTests
{
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

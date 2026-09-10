using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class CoefficientAndCmTests
{
    [Theory]
    [InlineData("37a1"), InlineData("11a1"), InlineData("48a3"), InlineData("389a1"), InlineData("5077a1")]
    public void PublicCoefficientsAndMetadataAgreeWithStoredData(string label)
    {
        var data = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture(label)); var e = data.GlobalMinimalModel;
        var coefficients = e.GetFourierCoefficients(data.FourierCoefficients.Count - 1);
        Assert.Equal(data.FourierCoefficients, coefficients.Select(x => new BigInteger(x)));
        int index = 0;
        for (int prime = 2; index < data.PrimeFourierCoefficients.Count; prime++)
        {
            if (!NativeNumberTheory.IsPrime(prime, default)) continue;
            Assert.Equal(data.PrimeFourierCoefficients[index++], e.GetFrobeniusTrace(prime));
            if (e.Discriminant.Num % prime != 0) Assert.Equal(prime + 1 - e.GetFrobeniusTrace(prime), e.CountPoints(prime));
        }
        Assert.Equal(e.CmDiscriminant, data.CmDiscriminant);
        Assert.Equal(e.TorsionPoints.Count(), data.TorsionOrder);
        Assert.Equal(data.IsogenyClassSize, data.IsogenyMatrix.Count);
        Assert.NotNull(data.FaltingsHeight); Assert.NotNull(data.StableFaltingsHeight);
        Assert.NotNull(data.LeadingLValue); Assert.NotNull(data.AnalyticShaOrder); Assert.NotNull(data.IntegralPointXCoordinates);
        Assert.True(data.ModularDegree > 0); Assert.True(data.ManinConstant > 0);
    }

    [Theory]
    [InlineData(0L, -3), InlineData(1728L, -4), InlineData(-3375L, -7), InlineData(8000L, -8)]
    [InlineData(-32768L, -11), InlineData(54000L, -12), InlineData(287496L, -16), InlineData(-884736L, -19)]
    [InlineData(-12288000L, -27), InlineData(16581375L, -28), InlineData(-884736000L, -43)]
    [InlineData(-147197952000L, -67), InlineData(-262537412640768000L, -163)]
    public void RationalCmClassificationIsExactAndModelInvariant(long j, int discriminant)
    {
        var e = EllipticCurveQ.FromJInvariant(j);
        Assert.True(e.HasComplexMultiplication); Assert.Equal(discriminant, e.CmDiscriminant);
        Assert.Equal(discriminant, e.QuadraticTwist(-7).ChangeModel(new BigRational(-2, 3), 5, -1, 4).Target.CmDiscriminant);
        Assert.False(EllipticCurveQ.FromJInvariant(new BigRational(j * new BigInteger(2) + 1, 2)).HasComplexMultiplication);
    }

    [Fact]
    public void CoefficientsValidateInputsAndPreserveModelInvariance()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0); var changed = e.ChangeModel(new BigRational(2, 3), 5, -1, 7).Target;
        Assert.Equal(e.GetFourierCoefficients(100), changed.GetFourierCoefficients(100));
        Assert.Equal(new long[] { 0 }, e.GetFourierCoefficients(0)); Assert.Equal(1, e.GetFourierCoefficient(1));
        Assert.Equal(e.GetFourierCoefficients(50)[49], e.GetFourierCoefficient(49));
        Assert.Throws<ArgumentException>(() => e.CountPoints(37));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetFrobeniusTrace(9));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetFourierCoefficient(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => e.GetFourierCoefficients(-1));
        Assert.Throws<ArithmeticException>(() => e.GetFrobeniusTrace(101, 100));
        Assert.Throws<ArithmeticException>(() => e.GetFourierCoefficients(100, 100));
        Assert.Throws<OperationCanceledException>(() => e.GetFourierCoefficients(1, cancellationToken: new CancellationToken(true)));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).GetFourierCoefficients(0));
        Assert.Throws<InvalidOperationException>(() => new EllipticCurveQ(0, 0, 0, 0, 0).CmDiscriminant);
    }
}

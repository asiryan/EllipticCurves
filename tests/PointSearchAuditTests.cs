using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class PointSearchAuditTests
{
    [Theory]
    [InlineData(1048573)] // Immediately below the square-root table cap.
    [InlineData(1048583)] // Immediately above it.
    [InlineData(7340033)] // p-1 = 7*2^20: exercises the long Tonelli-Shanks chain.
    public void PrimeEnumerationAtTableBoundaryMatchesIndependentCharacters(int prime)
    {
        const int bound = 512;
        var curve = new EllipticCurveFp(prime, -1, 2, -3, 4, 5);
        var points = curve.Points(prime).TakeWhile(point => point.IsInfinity || point.X < bound).ToArray();
        Assert.True(points[0].IsInfinity);
        Assert.Equal(points.Length, points.Distinct().Count());
        Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));
        var byX = points.Where(point => !point.IsInfinity).ToLookup(point => point.X);
        for (int x = 0; x < bound; x++)
        {
            BigInteger b = curve.A1 * x + curve.A3;
            BigInteger rhs = ((x + curve.A2) * x + curve.A4) * x + curve.A6;
            var discriminant = Mod(b * b + 4 * rhs, prime);
            int expected = discriminant.IsZero ? 1
                : BigInteger.ModPow(discriminant, (prime - 1) / 2, prime).IsOne ? 2 : 0;
            Assert.Equal(expected, byX[new BigInteger(x)].Count());
        }
    }

    [Theory]
    [InlineData(16777199)] // Below the residue-bitset cap.
    [InlineData(16777259)] // Above it: exact Jacobi fallback.
    public void CardinalityAcrossResidueTableBoundaryMatchesExactSymmetry(int prime)
    {
        Assert.Equal(3, prime % 4);
        // For y^2=x^3-x, x and -x have opposite nonzero quadratic characters
        // when p=3 mod 4. Their contributions cancel, proving #E(F_p)=p+1.
        var curve = new EllipticCurveFp(prime, 0, 0, 0, -1, 0);
        Assert.Equal(new BigInteger(prime) + 1, curve.CountPoints(prime));
    }

    [Theory]
    [InlineData(3), InlineData(7), InlineData(101), InlineData(1009)]
    public void HugeSignedCoefficientsAgreeWithDirectEulerCounting(int prime)
    {
        BigInteger[] reduced = { 1, 2, 3, 4, 5 };
        var huge = reduced.Select((value, index) => value + (index % 2 == 0 ? 1 : -1) *
            prime * (BigInteger.One << (4096 + 37 * index))).ToArray();
        BigInteger expected = 1;
        for (int x = 0; x < prime; x++)
        {
            var b = reduced[0] * x + reduced[2];
            var rhs = ((x + reduced[1]) * x + reduced[3]) * x + reduced[4];
            var discriminant = Mod(b * b + 4 * rhs, prime);
            expected += discriminant.IsZero ? 1
                : BigInteger.ModPow(discriminant, (prime - 1) / 2, prime).IsOne ? 2 : 0;
        }
        Assert.Equal(expected, PrimeFieldPointSearch.Count(prime, huge[0], huge[1], huge[2], huge[3], huge[4], default));
        Assert.Equal(expected, PrimeFieldPointSearch.CountWithoutTable(prime, huge[0], huge[1], huge[2], huge[3], huge[4], default));
        var rational = new EllipticCurveQ(new(huge[0]), new(huge[1]), new(huge[2]), new(huge[3]), new(huge[4]));
        Assert.Equal((long)(prime + 1 - expected), AnalyticCoefficients.Trace(rational, prime, default));
        if (rational.Discriminant.Num % prime != 0)
            Assert.Equal(expected, new EllipticCurveFp(prime, huge[0], huge[1], huge[2], huge[3], huge[4]).CountPoints(prime));
    }

    [Fact]
    public void ExtensionFallbackSquaresAndNonsquaresAgreeWithPrimeFieldNorm()
    {
        const int prime = 1031;
        var field = FiniteFieldTests.Field(prime, new[] { 1, 0, 1 });
        Assert.True(field.Order > FiniteFieldQuadraticSolver.MaximumTableOrder);
        var solver = new FiniteFieldQuadraticSolver(field, default);
        var random = new Random(20261120);
        int squares = 0, nonsquares = 0;
        for (int sample = 0; sample < 200; sample++)
        {
            int a = sample == 0 ? 0 : random.Next(prime);
            int b = sample == 0 ? 0 : random.Next(prime);
            var value = field.CreateElement(a, b);
            bool expected = IsSquareByNorm(value, prime);
            Assert.Equal(expected, solver.IsSquare(value, default));
            Assert.Equal(expected, solver.TrySquareRoot(value, out var root, default));
            if (expected) { squares++; Assert.Equal(value, root * root); }
            else nonsquares++;
        }
        Assert.True(squares > 1 && nonsquares > 0);
    }

    [Fact]
    public void PublicExtensionFallbackEnumerationMatchesNormCharacters()
    {
        const int prime = 1031, bound = 64;
        var field = FiniteFieldTests.Field(prime, new[] { 1, 0, 1 });
        var curve = new EllipticCurveFq(field, field.Generator, field.One, field.One, field.Generator, field.One);
        var points = curve.Points((long)field.Order).TakeWhile(point => point.IsInfinity ||
            point.X.Coordinates.Length <= 1 && (point.X.IsZero || point.X.Coordinates[0] < bound)).ToArray();
        Assert.True(points[0].IsInfinity);
        Assert.Equal(points.Length, points.Distinct().Count());
        Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));
        var byX = points.Where(point => !point.IsInfinity).ToLookup(point => point.X);
        for (int code = 0; code < bound; code++)
        {
            var x = field.CreateElement(code);
            var b = curve.A1 * x + curve.A3;
            var rhs = ((x + curve.A2) * x + curve.A4) * x + curve.A6;
            var discriminant = b * b + 4 * rhs;
            int expected = discriminant.IsZero ? 1 : IsSquareByNorm(discriminant, prime) ? 2 : 0;
            Assert.Equal(expected, byX[x].Count());
        }
    }

    [Fact]
    public void PrimeSubfieldJDoesNotTriggerFalseCoefficientDescent()
    {
        var field = FiniteFieldTests.Field(3, new[] { 1, 0, 1 });
        var curve = new EllipticCurveFq(field, field.Zero, field.Zero, field.Zero, field.Generator, field.Zero);
        var baseChange = new EllipticCurveFq(field, 0, 0, 0, 1, 0);
        Assert.Equal(field.Zero, curve.JInvariant);
        Assert.Equal(curve.JInvariant, baseChange.JInvariant);
        // For x=a+b*t, x^3+t*x=(a-b)(1+t). Its two nonzero values
        // are nonsquares in F_9, so only three affine points and O remain.
        Assert.Equal(new BigInteger(4), curve.CountPoints());
        Assert.Equal(new BigInteger(16), baseChange.CountPoints());
        var expected = BrutePoints(curve);
        Assert.Equal(4, expected.Count);
        Assert.True(expected.SetEquals(curve.Points()));
    }

    [Fact]
    public void EveryGeneralWeierstrassModelOverF4MatchesCoordinatePairSearch()
    {
        var field = FiniteFieldTests.Field(2, new[] { 1, 1, 1 });
        var elements = field.Elements().ToArray();
        int nonsingular = 0;
        for (int code = 0; code < 1024; code++)
        {
            int remaining = code;
            var a = new FiniteFieldElement[5];
            for (int i = 0; i < a.Length; i++) { a[i] = elements[remaining % 4]; remaining /= 4; }
            EllipticCurveFq curve;
            try { curve = new EllipticCurveFq(field, a[0], a[1], a[2], a[3], a[4]); }
            catch (ArgumentException) { continue; }
            nonsingular++;
            var expected = BrutePoints(curve);
            var actual = curve.Points(4).ToArray();
            Assert.Equal(expected.Count, actual.Length);
            Assert.True(expected.SetEquals(actual));
            Assert.Equal(new BigInteger(expected.Count), curve.CountPoints(4));
        }
        Assert.Equal(768, nonsingular);
    }

    [Theory]
    [InlineData(3), InlineData(5), InlineData(7)]
    public void DegreeOneNonzeroGeneratorPresentationAgreesWithPrimeField(int prime)
    {
        var field = FiniteFieldTests.Field(prime, new[] { 2, 1 });
        Assert.Equal(field.CreateElement(-2), field.Generator);
        var curve = new EllipticCurveFq(field, 0, 0, 1, -1, 0);
        var baseline = new EllipticCurveFp(prime, 0, 0, 1, -1, 0);
        Assert.Equal(baseline.CountPoints(), curve.CountPoints());
        var expected = baseline.Points().Select(point => point.IsInfinity ? (-1, -1) : ((int)point.X, (int)point.Y)).ToHashSet();
        int Constant(FiniteFieldElement value) => value.IsZero ? 0 : (int)value.Coordinates[0];
        var actual = curve.Points().Select(point => point.IsInfinity ? (-1, -1) : (Constant(point.X), Constant(point.Y))).ToArray();
        Assert.Equal(expected.Count, actual.Length);
        Assert.True(expected.SetEquals(actual));
    }

    private static bool IsSquareByNorm(FiniteFieldElement value, int prime)
    {
        // In F_p[t]/(t^2+1), Norm(a+b*t)=a^2+b^2. Its Legendre symbol
        // equals the quadratic character over F_(p^2), independently of the solver.
        var a = value.Coordinates.Length > 0 ? value.Coordinates[0] : BigInteger.Zero;
        var b = value.Coordinates.Length > 1 ? value.Coordinates[1] : BigInteger.Zero;
        var norm = Mod(a * a + b * b, prime);
        return norm.IsZero || BigInteger.ModPow(norm, (prime - 1) / 2, prime).IsOne;
    }

    private static HashSet<EllipticCurvePointFq> BrutePoints(EllipticCurveFq curve)
    {
        var result = new HashSet<EllipticCurvePointFq> { EllipticCurvePointFq.Infinity };
        var elements = curve.Field.Elements().ToArray();
        var squares = elements.Select(value => value * value).ToArray();
        foreach (var x in elements)
        {
            var xSquared = x * x;
            var rhs = xSquared * x + curve.A2 * xSquared + curve.A4 * x + curve.A6;
            var linear = curve.A1 * x + curve.A3;
            for (int y = 0; y < elements.Length; y++)
                if (squares[y] + linear * elements[y] == rhs)
                    result.Add(curve.CreatePoint(x, elements[y]));
        }
        return result;
    }

    private static BigInteger Mod(BigInteger value, BigInteger prime)
    {
        var residue = value % prime;
        return residue.Sign < 0 ? residue + prime : residue;
    }
}

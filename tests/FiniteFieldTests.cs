using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class FiniteFieldTests
{
    public static IEnumerable<object[]> Fields()
    {
        yield return new object[] { 2, new int[] { 1, 1, 1 } };
        yield return new object[] { 2, new int[] { 1, 1, 0, 1 } };
        yield return new object[] { 2, new int[] { 1, 1, 0, 0, 1 } };
        yield return new object[] { 3, new int[] { 1, 0, 1 } };
        yield return new object[] { 3, new int[] { 1, 2, 0, 1 } };
        yield return new object[] { 5, new int[] { 2, 0, 1 } };
        yield return new object[] { 7, new int[] { 1, 0, 1 } };
        yield return new object[] { 5, new int[] { 1, 1, 0, 1 } };
    }
    internal static FiniteField Field(int prime, int[] polynomial) => new(prime, polynomial.Select(x => new BigInteger(x)).ToArray());
    internal static FiniteFieldElement Decode(FiniteField field, int code)
    {
        var coefficients = new List<BigInteger>();
        while (code > 0) { coefficients.Add(code % field.Characteristic); code /= (int)field.Characteristic; }
        return field.CreateElement(coefficients.ToArray());
    }
    public static IEnumerable<object[]> ReferenceRows() => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "finite-fields.csv")).Select(s => new object[] { s });

    [Theory, MemberData(nameof(ReferenceRows))]
    public void ArithmeticAgreesWithPari(string row)
    {
        var v = row.Split(','); var field = Field(int.Parse(v[0]), v[1].Split(':').Select(int.Parse).ToArray());
        FiniteFieldElement Element(int index) => Decode(field, int.Parse(v[index]));
        var a = Element(2); var b = Element(3);
        Assert.Equal(Element(4), a + b); Assert.Equal(Element(5), a * b);
        if (v[6] != "-1") Assert.Equal(Element(6), a.Inverse());
    }

    [Theory, MemberData(nameof(Fields))]
    public void FieldLawsFrobeniusAndInverses(int prime, int[] polynomial)
    {
        var field = Field(prime, polynomial); var elements = field.Elements().ToArray();
        Assert.Equal(field.Order, elements.Length); Assert.Equal(elements.Length, elements.Distinct().Count());
        var random = new Random(prime + polynomial.Length);
        foreach (var a in elements)
        {
            Assert.Equal(a, a.Pow(field.Order)); Assert.Equal(field.One, a.Pow(0));
            Assert.Equal(field.Zero, a + -a); Assert.Equal(field.Zero, prime * a);
            if (!a.IsZero) { Assert.Equal(field.One, a * a.Inverse()); Assert.Equal(a.Inverse() * a.Inverse(), a.Pow(-2)); }
            foreach (var b in elements.Take(12))
            {
                var c = elements[random.Next(elements.Length)];
                Assert.Equal(a, a + b - b); Assert.Equal(a * b, b * a);
                Assert.Equal(a * (b + c), a * b + a * c);
                Assert.Equal((a * b) * c, a * (b * c));
                Assert.Equal((a + b).Pow(prime), a.Pow(prime) + b.Pow(prime));
                if (!b.IsZero) Assert.Equal(a, a / b * b);
            }
        }
    }

    // Independent exhaustive trial division, deliberately not Rabin's criterion.
    private static bool HasProperFactor(int[] f, int prime)
    {
        for (int degree = 1; degree <= (f.Length - 1) / 2; degree++)
        for (int code = 0; code < (int)Math.Pow(prime, degree); code++)
        {
            var g = new int[degree + 1]; g[degree] = 1; int value = code;
            for (int i = 0; i < degree; i++) { g[i] = value % prime; value /= prime; }
            var remainder = (int[])f.Clone();
            for (int i = f.Length - 1; i >= degree; i--)
            {
                int factor = remainder[i];
                for (int j = 0; j <= degree; j++) remainder[i - degree + j] = (remainder[i - degree + j] - factor * g[j] % prime + prime) % prime;
            }
            if (remainder.All(x => x == 0)) return true;
        }
        return false;
    }

    [Theory]
    [InlineData(2, 6), InlineData(3, 4), InlineData(5, 3)]
    public void IrreducibilityAgreesWithExhaustiveTrialFactorization(int prime, int maxDegree)
    {
        for (int degree = 1; degree <= maxDegree; degree++)
        for (int code = 0; code < (int)Math.Pow(prime, degree); code++)
        {
            var f = new int[degree + 1]; f[degree] = 1; int value = code;
            for (int i = 0; i < degree; i++) { f[i] = value % prime; value /= prime; }
            if (HasProperFactor(f, prime)) Assert.Throws<ArgumentException>(() => Field(prime, f));
            else Assert.Equal(degree, Field(prime, f).Degree);
        }
    }

    [Fact]
    public void PresentationsAreCanonicalImmutableAndNotImplicitlyInterchangeable()
    {
        var coefficients = new BigInteger[] { 2, 0, 2, 0 }; var field = new FiniteField(3, coefficients);
        var same = Field(3, new[] { 1, 0, 1 }); var other = Field(3, new[] { 2, 1, 1 });
        coefficients[0] = 0;
        Assert.Equal(same, field); Assert.Equal(same.GetHashCode(), field.GetHashCode()); Assert.NotEqual(field, other);
        var input = new BigInteger[] { 4, -1, 1 }; var a = field.CreateElement(input); input[0] = 7;
        Assert.Equal(field.CreateElement(0, 2), a);
        Assert.Equal(field.Generator, same.Generator); Assert.Equal(field.One, field.Generator * same.Generator + 2);
        Assert.Throws<NotSupportedException>(() => ((IList<BigInteger>)field.Modulus)[0] = 0);
        Assert.Throws<NotSupportedException>(() => ((IList<BigInteger>)a.Coefficients)[0] = 0);
        Assert.Throws<ArgumentException>(() => field.Add(field.One, other.One));
        Assert.Throws<ArgumentException>(() => field.Multiply(field.Zero, other.Zero));
        Assert.Throws<ArgumentException>(() => field.Divide(field.One, default));
        Assert.Throws<ArgumentException>(() => new EllipticCurveFq(field, other.One, field.Zero, field.One, field.Zero, field.Zero));
        Assert.False(default(FiniteFieldElement).IsZero);
        Assert.NotEqual(field.Zero, default(FiniteFieldElement));
        Assert.Throws<InvalidOperationException>(() => default(FiniteFieldElement).Pow(0));
    }

    [Fact]
    public void LargeCharacteristicLinearFieldsAndInputLimits()
    {
        var prime = (BigInteger.One << 61) - 1;
        var field = new FiniteField(prime, new BigInteger[] { 1, 0, 1 });
        Assert.True(field.Order > ulong.MaxValue);
        Assert.Equal(-field.One, field.Generator.Pow(2));
        var a = field.Generator + 123; Assert.Equal(field.One, a * a.Inverse());
        Assert.Throws<ArithmeticException>(() => field.Elements().ToArray());
        var linear = Field(5, new[] { 2, 1 });
        Assert.Equal(linear.CreateElement(3), linear.Generator);
        Assert.Equal(linear.CreateElement(4), linear.CreateElement(1, 2, 3));
        Assert.Throws<ArgumentException>(() => Field(2, new[] { 1, 0, 1 }));
        Assert.Throws<ArgumentException>(() => Field(2, new[] { 1, 1, 1, 1, 1, 1, 1 })); // product of two cubics, with no roots
        Assert.Throws<ArgumentException>(() => Field(5, new[] { 1, 5 }));
        Assert.Throws<ArgumentException>(() => Field(5, Array.Empty<int>()));
        Assert.Throws<ArgumentOutOfRangeException>(() => Field(9, new[] { 1, 0, 1 }));
        Assert.Throws<ArgumentNullException>(() => new FiniteField(3, null));
        Assert.Throws<DivideByZeroException>(() => field.Zero.Inverse());
        Assert.Throws<DivideByZeroException>(() => field.One / field.Zero);
        Assert.Throws<DivideByZeroException>(() => field.Zero.Pow(-1));
        Assert.Throws<ArgumentNullException>(() => field.CreateElement(null));
        Assert.Throws<ArgumentOutOfRangeException>(() => linear.Elements(-1).ToArray());
        Assert.Throws<ArithmeticException>(() => linear.Elements(4).ToArray());
        var cancelled = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => new FiniteField(2, new BigInteger[] { 1, 1, 1 }, cancelled));
        Assert.Throws<OperationCanceledException>(() => a.Pow(0, cancelled));
        Assert.Throws<OperationCanceledException>(() => field.Multiply(a, a, cancelled));
        Assert.Throws<OperationCanceledException>(() => field.Elements(cancellationToken: cancelled).ToArray());
    }
}

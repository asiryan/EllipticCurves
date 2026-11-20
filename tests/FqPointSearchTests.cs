using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class FqPointSearchTests
{
    public static IEnumerable<object[]> SmallFields() => FiniteFieldTests.Fields()
        .Concat(new[] { new object[] { 2, new int[] { 0, 1 } }, new object[] { 3, new int[] { 0, 1 } } });

    public static IEnumerable<object[]> BinaryFields() => SmallFields().Where(row => (int)row[0] == 2)
        .Concat(new[] { new object[] { 2, new int[] { 1, 1, 0, 1, 1, 0, 0, 0, 1 } } });

    [Theory, MemberData(nameof(BinaryFields))]
    public void BinaryPolynomialInverseSquaresMatchGeneralFieldInversion(int prime, int[] polynomial)
    {
        var field = FiniteFieldTests.Field(prime, polynomial);
        var solver = new FiniteFieldQuadraticSolver(field, default);
        foreach (var value in field.Elements().Where(value => !value.IsZero))
        {
            var inverse = field.Inverse(value);
            var actual = solver.InverseSquare(value, default);
            Assert.Equal(inverse * inverse, actual);
            Assert.Equal(field.One, value * value * actual);
        }
        Assert.Throws<DivideByZeroException>(() => solver.InverseSquare(field.Zero, default));
        Assert.Throws<OperationCanceledException>(() => solver.InverseSquare(field.One, new CancellationToken(true)));
    }

    [Fact]
    public void BinarySolversWorkAcrossMachineWordBoundaries()
    {
        var polynomial = new BigInteger[129];
        foreach (int bit in new[] { 0, 1, 2, 7, 128 }) polynomial[bit] = 1;
        var field = new FiniteField(2, polynomial);
        var solver = new FiniteFieldQuadraticSolver(field, default);
        var coefficients = new BigInteger[128];
        foreach (int bit in new[] { 0, 65, 127 }) coefficients[bit] = 1;
        var value = field.CreateElement(coefficients);
        var square = value * value;
        Assert.Equal(field.One, square * solver.InverseSquare(value, default));
        Assert.True(solver.TrySquareRoot(square, default, out var root));
        Assert.Equal(square, root * root);
        var image = square + value;
        Assert.True(solver.TryArtinSchreierRoot(image, default, out root));
        Assert.Equal(image, root * root + root);
    }

    [Theory, MemberData(nameof(SmallFields))]
    public void GeneralEquationsMatchIndependentCoordinatePairSearch(int prime, int[] polynomial)
    {
        var field = FiniteFieldTests.Field(prime, polynomial);
        var elements = field.Elements().ToArray();
        var squares = elements.Select(value => value * value).ToArray();
        var random = new Random(61207 + prime + field.Degree);
        for (int sample = 0; sample < 10; sample++)
        {
            var a = Enumerable.Range(0, 5).Select(_ => elements[random.Next(elements.Length)]).ToArray();
            if (sample < 3) { a[0] = field.Zero; a[2] = sample == 0 ? field.Zero : field.One; }
            EllipticCurveFq curve;
            try { curve = new EllipticCurveFq(field, a[0], a[1], a[2], a[3], a[4]); }
            catch (ArgumentException) { continue; }
            var expected = new HashSet<EllipticCurvePointFq> { EllipticCurvePointFq.Infinity };
            foreach (var x in elements)
            {
                var xSquared = x * x;
                var rhs = xSquared * x + a[1] * xSquared + a[3] * x + a[4];
                var linear = a[0] * x + a[2];
                for (int y = 0; y < elements.Length; y++)
                    if (squares[y] + linear * elements[y] == rhs)
                        expected.Add(curve.CreatePoint(x, elements[y]));
            }
            var actual = curve.Points((long)field.Order).ToArray();
            Assert.True(actual[0].IsInfinity);
            Assert.Equal(expected.Count, actual.Length);
            Assert.Equal(actual.Length, actual.Distinct().Count());
            Assert.True(expected.SetEquals(actual));
            Assert.Equal(new BigInteger(expected.Count), curve.CountPoints((long)field.Order));
        }
    }

    [Theory, MemberData(nameof(SmallFields))]
    public void QuadraticSolversMatchExhaustiveImages(int prime, int[] polynomial)
    {
        var field = FiniteFieldTests.Field(prime, polynomial);
        var elements = field.Elements().ToArray();
        var squares = elements.Select(value => value * value).ToHashSet();
        var artinSchreier = prime == 2 ? elements.Select(value => value * value + value).ToHashSet() : null;
        var solver = new FiniteFieldQuadraticSolver(field, default);
        foreach (var value in elements)
        {
            Assert.Equal(squares.Contains(value), solver.IsSquare(value, default));
            bool exists = solver.TrySquareRoot(value, default, out var root);
            Assert.Equal(squares.Contains(value), exists);
            if (exists) Assert.Equal(value, root * root);
            if (prime != 2) continue;
            Assert.Equal(artinSchreier.Contains(value), solver.HasArtinSchreierRoot(value, default));
            exists = solver.TryArtinSchreierRoot(value, default, out root);
            Assert.Equal(artinSchreier.Contains(value), exists);
            if (exists) Assert.Equal(value, root * root + root);
        }
    }

    [Theory]
    [InlineData(1000000007), InlineData(1000000009)]
    public void LargePrimeSquareRootsUseExactBoundedMemoryFallback(int prime)
    {
        var field = FiniteFieldTests.Field(prime, new[] { 0, 1 });
        Assert.True(field.Order > FiniteFieldQuadraticSolver.MaximumTableOrder);
        var solver = new FiniteFieldQuadraticSolver(field, default);
        foreach (int code in new[] { 0, 1, 2, 31, 197 })
        {
            var value = field.CreateElement(code); var square = value * value;
            Assert.True(solver.TrySquareRoot(square, default, out var root));
            Assert.Equal(square, root * root);
        }
        int nonSquare = 2;
        while (BigInteger.ModPow(nonSquare, (prime - 1) / 2, prime) == 1) nonSquare++;
        Assert.False(solver.TrySquareRoot(field.CreateElement(nonSquare), default, out _));
    }

    [Theory]
    [InlineData(1031)]
    [InlineData(1000000007)]
    public void EvenDegreeFallbackFindsNonresiduesOutsidePrimeSubfield(int prime)
    {
        var field = FiniteFieldTests.Field(prime, new[] { 1, 0, 1 });
        Assert.True(field.Order > FiniteFieldQuadraticSolver.MaximumTableOrder);
        var solver = new FiniteFieldQuadraticSolver(field, default);
        // Searching the p constants first would never find a nonresidue in
        // even degree, and is infeasible for the billion-element subfield.
        using var timeout = new CancellationTokenSource(TimeSpan.FromSeconds(10));
        foreach (var value in new[] { field.One, field.Generator, field.Generator + 13 })
        {
            var square = value * value;
            Assert.True(solver.TrySquareRoot(square, timeout.Token, out var root));
            Assert.Equal(square, root * root);
        }
    }

    [Fact]
    public void ModerateExtensionEnumerationNeedsOnlyOneWorkUnitPerAbscissa()
    {
        var field = FiniteFieldTests.Field(101, new[] { 2, 0, 1 });
        var curve = new EllipticCurveFq(field, field.Zero, field.Zero, field.One, field.Generator, field.One);
        long work = (long)field.Order;
        Assert.Throws<ArithmeticException>(() => curve.Points(work - 1).First());
        var points = curve.Points(work).ToArray();
        Assert.Equal(points.Length, points.Distinct().Count());
        Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));
        Assert.Equal(new BigInteger(points.Length), curve.CountPoints(work));
    }

    [Theory]
    [InlineData(new int[] { 1, 1, 0, 1, 1, 0, 0, 0, 1 })]
    [InlineData(new int[] { 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 })]
    public void CharacteristicTwoLargerEvenAndOddDegreesEnumerateExactly(int[] polynomial)
    {
        var field = FiniteFieldTests.Field(2, polynomial);
        foreach (var a1 in new[] { field.Zero, field.One, field.Generator })
        {
            var curve = new EllipticCurveFq(field, a1, field.Zero, field.One, field.Generator, field.One);
            var points = curve.Points((long)field.Order).ToArray();
            Assert.Equal(points.Length, points.Distinct().Count());
            Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));
            Assert.Equal(new BigInteger(points.Length), curve.CountPoints((long)field.Order));
        }
    }

    [Fact]
    public void CachedCountsStillValidateBudgetsAndCancellation()
    {
        var field = FiniteFieldTests.Field(3, new[] { 1, 0, 1 });
        var curve = new EllipticCurveFq(field, 0, 0, 1, -1, 0);
        var cancelled = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(cancellationToken: cancelled));
        var count = curve.CountPoints();
        Assert.Equal(count, curve.CountPoints(9));
        Assert.Throws<ArithmeticException>(() => curve.CountPoints(8));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.CountPoints(-1));
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(cancellationToken: cancelled));
        Assert.Throws<OperationCanceledException>(() => curve.Points(cancellationToken: cancelled).First());
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.Points(-1).First());
        using var source = new CancellationTokenSource();
        using var enumerator = curve.Points(cancellationToken: source.Token).GetEnumerator();
        Assert.True(enumerator.MoveNext());
        Assert.True(enumerator.Current.IsInfinity);
        source.Cancel();
        Assert.Throws<OperationCanceledException>(() => enumerator.MoveNext());
    }

    [Fact]
    public void CancellationAfterQuadraticPreparationLeavesCountRetryUsable()
    {
        var field = FiniteFieldTests.Field(101, new[] { 2, 0, 1 });
        var curve = new EllipticCurveFq(field, field.Zero, field.Zero, field.One, field.Generator, field.One);
        using var source = new CancellationTokenSource();
        using var enumerator = curve.Points(cancellationToken: source.Token).GetEnumerator();
        Assert.True(enumerator.MoveNext());
        Assert.True(enumerator.Current.IsInfinity);
        // The first affine point guarantees that quadratic preparation finished.
        // Cancel at this boundary instead of racing a timer against the machine.
        Assert.True(enumerator.MoveNext());
        Assert.False(enumerator.Current.IsInfinity);
        source.Cancel();
        Assert.Throws<OperationCanceledException>(() => enumerator.MoveNext());
        Assert.Throws<OperationCanceledException>(() => curve.CountPoints(cancellationToken: source.Token));
        var result = curve.CountPoints();
        Assert.Equal(new BigInteger(curve.Points().Count()), result);
    }
}

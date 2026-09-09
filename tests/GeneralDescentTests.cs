using System.Numerics;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class GeneralDescentTests
{
    private static DescentBudget Budget() => new(new RankComputationOptions { MaxDescentWork = 20000000 }, default);

    public static IEnumerable<object[]> PariCases => File.ReadLines(Path.Combine(AppContext.BaseDirectory, "Fixtures", "general-descent.csv"))
        .Skip(8).Where(x => !string.IsNullOrWhiteSpace(x)).Select(x => new object[] { x });

    [Theory]
    [MemberData(nameof(PariCases))]
    public void GeneralSelmerDimensionsMatchIndependentPari(string row)
    {
        var a = row.Split(',').Select(int.Parse).ToArray();
        var e = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        var result = e.GetRankBounds(new RankComputationOptions { PreferGeneralTwoDescent = true });
        Assert.True(result.TwoSelmerDimension == a[7], $"{row}: {result}, Selmer={result.TwoSelmerDimension}, work={result.DescentWork}, {result.Reason}");
        Assert.True(result.LowerBound <= a[6], $"{row}: {result}");
        Assert.True(result.UpperBound >= a[5], $"{row}: {result}");
    }

    [Theory]
    [InlineData(0, -1, 1, -10, -20, 0)]
    [InlineData(0, 0, 1, -1, 0, 1)]
    [InlineData(0, 1, 1, -2, 0, 2)]
    [InlineData(0, 0, 1, -7, 6, 3)]
    [InlineData(1, -1, 0, -79, 289, 4)]
    public void ProvesRanksWithoutRationalTwoTorsion(int a1, int a2, int a3, int a4, int a6, int rank)
    {
        var e = new EllipticCurveQ(a1, a2, a3, a4, a6);
        var options = new RankComputationOptions { MaxDescentWork = 10000000 };
        var result = e.GetRankBounds(options);
        Assert.True(result.ExactRank == rank, $"{e}: {result}, {result.Reason}, work={result.DescentWork}");
        Assert.Equal(rank, result.TwoSelmerDimension);
        Assert.True(result.UsedGeneralTwoDescent);
        Assert.False(result.UsedTwoIsogenyDescent);
        var changed = NativeArithmeticTests.ChangeCoordinates(e, new BigRational(-2, 3), 4, -2, 3);
        var other = changed.GetRankBounds(options);
        Assert.Equal(result.ExactRank, other.ExactRank);
        Assert.Equal(result.TwoSelmerDimension, other.TwoSelmerDimension);
    }

    [Theory]
    [InlineData(0, 0, 1, -1, 0, 1)]
    [InlineData(0, 1, 1, -2, 0, 2)]
    [InlineData(0, 0, 1, -7, 6, 3)]
    public void ReductionCharactersProveHigherLowerBoundsBeforeDescent(int a1, int a2, int a3, int a4, int a6, int rank)
    {
        var result = new EllipticCurveQ(a1, a2, a3, a4, a6).GetRankBounds(new RankComputationOptions { MaxDescentWork = 0 });
        Assert.Equal(rank, result.LowerBound);
        Assert.Null(result.UpperBound);
        Assert.Null(result.TwoSelmerDimension);
        Assert.Contains("MaxDescentWork", result.Reason);
    }

    [Theory]
    [InlineData(0, 0, 1, 9, 9)]
    [InlineData(0, 1, 1, 6, 5)]
    public void NontrivialShaDoesNotBecomeAnExactRank(int a1, int a2, int a3, int a4, int a6)
    {
        // PARI independently proves rank 0 and computes 2-Selmer dimension 2.
        var result = new EllipticCurveQ(a1, a2, a3, a4, a6).GetRankBounds();
        Assert.Equal(0, result.LowerBound);
        Assert.Equal(2, result.UpperBound);
        Assert.Equal(2, result.TwoSelmerDimension);
        Assert.Null(result.ExactRank);
    }

    [Theory]
    [InlineData(1, 0)]
    [InlineData(5, 1)]
    [InlineData(34, 2)]
    public void GeneralDescentAccountsForFullRationalTwoTorsion(int n, int rank)
    {
        var result = new EllipticCurveQ(0, 0, 0, -n*n, 0).GetRankBounds(new RankComputationOptions { PreferGeneralTwoDescent = true, MaxDescentWork = 10000000 });
        Assert.Equal(rank + 2, result.TwoSelmerDimension);
        Assert.Equal(rank, result.ExactRank);
    }

    [Fact]
    public void IncompleteEnumerationNeverUsesThePartialClassCountAsUpperBound()
    {
        var e = new EllipticCurveQ(0, 0, 1, -7, 6);
        var complete = e.GetRankBounds();
        Assert.Equal(3, complete.ExactRank);
        foreach (long work in new[] { 1L, complete.DescentWork / 2, complete.DescentWork - 1 })
        {
            var partial = e.GetRankBounds(new RankComputationOptions { MaxDescentWork = work });
            Assert.Null(partial.UpperBound);
            Assert.Null(partial.TwoSelmerDimension);
            Assert.Null(partial.ExactRank);
            Assert.Equal(3, partial.LowerBound);
        }
        var fewClasses = e.GetRankBounds(new RankComputationOptions { MaxSquareClasses = 2 });
        Assert.Null(fewClasses.UpperBound);
        Assert.Null(fewClasses.TwoSelmerDimension);
        Assert.Contains("MaxSquareClasses", fewClasses.Reason);
    }

    [Fact]
    public void CancellationInterruptsAnActiveDescent()
    {
        using var cancellation = new CancellationTokenSource(TimeSpan.FromMilliseconds(50));
        Assert.Throws<OperationCanceledException>(() => new EllipticCurveQ(0, 0, 0, -100000007, -100000037)
            .GetRankBounds(new RankComputationOptions { SearchBound = 0, MaxDescentWork = long.MaxValue }, cancellation.Token));
    }

    [Fact]
    public void DependentMultiplesDoNotCreateExtraRank()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        var proof = new RationalPointRank(e, 0, Budget());
        var p = new EllipticCurvePoint(0, 0);
        for (int n = -6; n <= 6; n++) proof.Add(e.Multiply(p, n));
        Assert.Equal(1, proof.LowerBound);
        // Even multiples have trivial reduction characters, but infinite order is
        // still certified by exact additions and Mazur's torsion bound.
        var even = new RationalPointRank(e, 0, Budget());
        even.Add(e.Multiply(p, 2));
        Assert.Equal(0, even.ImageDimension);
        Assert.Equal(1, even.LowerBound);
    }

    [Fact]
    public void TorsionAndRepeatedPointsDoNotCreatePositiveRank()
    {
        var e = new EllipticCurveQ(0, -17, 0, 72, 0).GlobalMinimalModel;
        var proof = new RationalPointRank(e, 2, Budget());
        foreach (var p in e.TorsionPoints) { proof.Add(p); proof.Add(p); }
        Assert.Equal(0, proof.LowerBound);
    }

    [Fact]
    public void MissingPointSearchDoesNotProveRankZero()
    {
        var result = new EllipticCurveQ(0, 1, 1, -2, 0).GetRankBounds(new RankComputationOptions { SearchBound = 0 });
        Assert.Equal(0, result.LowerBound);
        Assert.Equal(2, result.UpperBound);
        Assert.Null(result.ExactRank);
        Assert.Equal(0, result.PointSearchWork);
    }

    [Fact]
    public void PointBudgetExhaustionKeepsCompletedSelmerBound()
    {
        var result = new EllipticCurveQ(0, 1, 1, -2, 0).GetRankBounds(new RankComputationOptions { MaxPointSearchWork = 0 });
        Assert.Equal(0, result.LowerBound);
        Assert.Equal(2, result.UpperBound);
        Assert.Contains("MaxPointSearchWork", result.Reason);
    }

    [Fact]
    public void ValidatesOptionsAndHonorsCancellation()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        foreach (var options in new RankComputationOptions[] {
            new() { SearchBound = -1 }, new() { MaxSquareClasses = 1 }, new() { MaxDescentWork = -1 },
            new() { MaxPointSearchWork = -1 }, new() { ReductionPrimeBound = 2 }, new() { ReductionPrimeBound = 10001 } })
            Assert.Throws<ArgumentOutOfRangeException>(() => e.GetRankBounds(options));
        Assert.Throws<ArgumentNullException>(() => e.GetRankBounds((RankComputationOptions)null!));
        Assert.Throws<OperationCanceledException>(() => e.GetRankBounds(new RankComputationOptions(), new CancellationToken(true)));
    }

    [Fact]
    public void ExactRootIsolationHandlesRepeatedInteriorRoots()
    {
        // (x+2)(x-1)^2(x-3) = x^4-3x^3-3x^2+11x-6.
        var roots = DescentPolynomial.RealRoots(new BigRational[] { -6, 11, -3, -3, 1 }, new BigRational(1, 100000), Budget());
        Assert.Equal(3, roots.Count);
        foreach (int expected in new[] { -2, 1, 3 }) Assert.Contains(roots, r => r.Lower <= expected && r.Upper >= expected);
        Assert.True(DescentPolynomial.HasRationalRoot(new BigRational[] { -6, 11, -3, -3, 1 }, Budget()));
        Assert.False(DescentPolynomial.HasRationalRoot(new BigRational[] { 4, 0, -4, 0, 1 }, Budget())); // (x^2-2)^2
        Assert.True(DescentPolynomial.HasRationalRoot(new BigRational[] { -1, 0, 0, 0, 81 }, Budget()));
    }

    [Fact]
    public void QuarticEquivalenceAndCovariantMapSurviveCoordinateChanges()
    {
        var q = new BinaryQuartic(2, 3, -5, 7, 9);
        var translated = new BinaryQuartic(q.A, 4*q.A+q.B, 6*q.A+3*q.B+q.C,
            4*q.A+3*q.B+2*q.C+q.D, q.A+q.B+q.C+q.D+q.E);
        Assert.True(q.Equivalent(translated, Budget()));
        Assert.True(q.Equivalent(q.Reverse(), Budget()));
        Assert.True(q.Equivalent(q, Budget()));
        var image = q.MapPoint(0, 1, 3);
        var jacobian = new EllipticCurveQ(0, 0, 0, new BigRational(-27*q.I), new BigRational(-27*q.J));
        Assert.True(jacobian.IsOnCurve(new EllipticCurvePoint(image.x, image.y)));
    }

    [Fact]
    public void EquivalenceHandlesReducibleResolventsAndZeroSexticSeminvariants()
    {
        // Counterexample in Cremona-Fisher (2009), section 5, to the older
        // equivalence criterion in Cremona's book. Both have I=592, J=-27776.
        var first = new BinaryQuartic(2, 0, -8, -8, 22);
        var second = new BinaryQuartic(3, 0, 22, -16, 3);
        Assert.False(first.Equivalent(second, Budget()));
        Assert.False(second.Equivalent(first, Budget()));
        // For J=0, negating a quartic need not preserve its covering class.
        var q = new BinaryQuartic(1, 0, -30, 0, 25);
        Assert.False(q.Equivalent(q.Scale(-1), Budget()));
        Assert.True(q.Equivalent(q, Budget()));
        Assert.True(q.Equivalent(q.Reverse(), Budget()));
        Assert.True(q.Scale(-1).Equivalent(new BinaryQuartic(1,0,-270,-2400,-5975), Budget()));
        // Trivial coverings with a branch point at infinity also compare correctly.
        var trivial = new BinaryQuartic(0, 1, 0, -1, 0);
        Assert.True(trivial.Equivalent(trivial.Reverse(), Budget()));
    }

    [Fact]
    public void LocalSolubilityIncludesInfinityAndGenuineLocalObstructions()
    {
        Assert.False(QuarticLocalSolubility.AtPrime(new BinaryQuartic(3, 0, 0, 0, 3), 2, Budget()));
        // No affine solution mod 3, but the infinity chart has the rational point (1:0:1).
        Assert.True(QuarticLocalSolubility.AtPrime(new BinaryQuartic(1, 0, -1, 0, 2), 3, Budget()));
        // Reichardt-Lind covering: everywhere locally soluble despite having no Q point.
        Assert.True(QuarticLocalSolubility.Everywhere(new BinaryQuartic(2, 0, 0, 0, -34), new BigInteger[] { 2, 17 }, Budget()));
        Assert.False(new BinaryQuartic(-1, 0, 0, 0, -1).HasRealPoint(Budget()));
        Assert.True(new BinaryQuartic(-1, 0, 10, 0, -1).HasRealPoint(Budget()));
    }

    [Fact]
    public void LocalSolubilitySurvivesDeepNonintegralChangesOfCoordinates()
    {
        for (int k = 1; k <= 9; k++)
        {
            var fourth = BigInteger.One << (4*k);
            Assert.False(QuarticLocalSolubility.AtPrime(new BinaryQuartic(3*fourth, 0, 0, 0, 3), 2, Budget()));
            Assert.True(QuarticLocalSolubility.AtPrime(new BinaryQuartic(2*fourth, 0, 0, 0, -34), 2, Budget()));
            Assert.True(QuarticLocalSolubility.AtPrime(new BinaryQuartic(2, 0, 0, 0, -34*fourth), 2, Budget()));
        }
    }

    [Fact]
    public void LocalAnswersRespectIndependentExhaustiveCongruenceObstructions()
    {
        var random = new Random(8675309);
        foreach (int p in new[] { 2, 3, 5, 7 })
        {
            int modulus = p == 2 ? 1024 : p*p*p*p;
            var squares = new bool[modulus];
            for (long y = 0; y < modulus; y++) squares[y*y % modulus] = true;
            int insoluble = 0;
            for (int trial = 0; trial <= 60; trial++)
            {
                long a = random.Next(1, 20), b = random.Next(-20, 21), c = random.Next(-20, 21), d = random.Next(-20, 21), e = random.Next(1, 20);
                if (trial == 0)
                {
                    // A nonsquare times (x^2 - nonsquare)^2 modulo odd p
                    // guarantees at least one obstruction at every prime tested.
                    long nonsquare = p == 3 || p == 5 ? 2 : 3;
                    a = nonsquare; b = d = 0; c = -2*nonsquare*nonsquare; e = nonsquare*nonsquare*nonsquare+p;
                    if (p == 2) { a = e = 3; c = 0; }
                }
                var q = new BinaryQuartic(a,b,c,d,e);
                if (4*q.I*q.I*q.I-q.J*q.J == 0) continue;
                bool residue = false;
                for (long x = 0; x < modulus && !residue; x++)
                {
                    long value = ((((a*x+b)*x+c)*x+d)*x+e) % modulus;
                    residue = squares[(value+modulus)%modulus];
                    if (x%p == 0)
                    {
                        value = ((((e*x+d)*x+c)*x+b)*x+a) % modulus;
                        residue |= squares[(value+modulus)%modulus];
                    }
                }
                bool local = QuarticLocalSolubility.AtPrime(q, p, Budget());
                if (!residue) { insoluble++; Assert.False(local); }
                if (local) Assert.True(residue);
            }
            Assert.True(insoluble > 0, $"No congruence obstruction tested at {p}.");
        }
    }

    [Fact]
    public void PrescribedRationalPointsAreNeverRejectedLocally()
    {
        var random = new Random(10101);
        foreach (int p in new[] { 2,3,5,7,11 })
        for (int k = 1; k <= 5; k++)
        {
            BigInteger a = random.Next(2,20), b = random.Next(-20,20), c = random.Next(-20,20), d = random.Next(-20,20);
            BigInteger y = BigInteger.Pow(p,k), e = y*y-a-b-c-d, scale = BigInteger.Pow(p,k);
            // The starting quartic has (x,y)=(1,p^k). Replace u by p^k*u;
            // its point now has denominator p^k and can require deep lifting.
            var q = new BinaryQuartic(a*BigInteger.Pow(scale,4),b*BigInteger.Pow(scale,3),c*scale*scale,d*scale,e);
            if (4*q.I*q.I*q.I-q.J*q.J == 0) continue;
            Assert.True(QuarticLocalSolubility.AtPrime(q,p,Budget()));
        }
    }
}

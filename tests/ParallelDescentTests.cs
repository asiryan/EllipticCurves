using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class ParallelDescentTests
{
    public static IEnumerable<object[]> KnownRanks() =>
        new[]
        {
            new[] { 0, -1, 1, -10, -20, 0 }, new[] { 0, 0, 1, -1, 0, 1 },
            new[] { 0, 1, 1, -2, 0, 2 }, new[] { 0, 0, 1, -7, 6, 3 },
            new[] { 1, -1, 0, -79, 289, 4 }
        }.SelectMany(curve => new[] { 2, 4 }.Select(workers => new object[] { curve, workers }));

    [Theory, MemberData(nameof(KnownRanks))]
    public void ParallelRanksAndCoordinateChangesMatchIndependentExamples(int[] a, int workers)
    {
        var curve = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        var options = new RankComputationOptions { MaxDescentWork = 10000000, MaxDegreeOfParallelism = workers };
        var result = curve.GetRankBounds(options);
        Assert.Equal(a[5], result.ExactRank);
        Assert.Equal(a[5], result.TwoSelmerDimension);
        var changed = NativeArithmeticTests.ChangeCoordinates(curve, new BigRational(-2, 3), 4, -2, 3);
        var other = changed.GetRankBounds(options);
        Assert.Equal(result.ExactRank, other.ExactRank);
        Assert.Equal(result.TwoSelmerDimension, other.TwoSelmerDimension);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(4)]
    public void ScreenshotCurveCompletesWithAnAdequateSharedBudget(int workers)
    {
        var curve = new EllipticCurveQ(0, 0, 0, new BigRational(-49, 50), new BigRational(1, 2));
        var result = curve.GetRankBounds(new RankComputationOptions
            { MaxDegreeOfParallelism = workers, MaxDescentWork = 20000000 });
        Assert.True(result.IsExact, result.Reason);
        Assert.Equal(1, result.ExactRank);
        Assert.Equal(1, result.TwoSelmerDimension);
        Assert.InRange(result.DescentWork, 1, 20000000);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    public void IncompleteParallelSearchNeverClaimsAnUpperBound(int workers)
    {
        var curve = new EllipticCurveQ(0, 0, 0, new BigRational(-49, 50), new BigRational(1, 2));
        foreach (long limit in new[] { 0L, 1, 128, 1000, 10000 })
        {
            var result = curve.GetRankBounds(new RankComputationOptions
                { MaxDegreeOfParallelism = workers, MaxDescentWork = limit });
            // The rank-one witness is found during descent on this curve, not
            // necessarily before the work allowance expires.
            Assert.InRange(result.LowerBound, 0, 1);
            Assert.Null(result.UpperBound);
            Assert.Null(result.TwoSelmerDimension);
            Assert.Null(result.ExactRank);
            Assert.Contains("MaxDescentWork", result.Reason);
            Assert.Equal(limit, result.DescentWork);
        }
        var rankThree = new EllipticCurveQ(0, 0, 1, -7, 6);
        var fewClasses = rankThree.GetRankBounds(new RankComputationOptions
            { MaxDegreeOfParallelism = workers, MaxSquareClasses = 2 });
        Assert.Equal(3, fewClasses.LowerBound);
        Assert.Null(fewClasses.UpperBound);
        Assert.Null(fewClasses.TwoSelmerDimension);
        Assert.Contains("MaxSquareClasses", fewClasses.Reason);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    public void MissingWitnessesAndNontrivialShaRemainUnproved(int workers)
    {
        foreach (bool noSearch in new[] { false, true })
        {
            var options = new RankComputationOptions { MaxDegreeOfParallelism = workers };
            if (noSearch) options.SearchBound = 0;
            else options.MaxPointSearchWork = 0;
            var result = new EllipticCurveQ(0, 1, 1, -2, 0).GetRankBounds(options);
            Assert.Equal(0, result.LowerBound);
            Assert.Equal(2, result.UpperBound);
            Assert.Null(result.ExactRank);
            Assert.Equal(0, result.PointSearchWork);
        }
        foreach (var curve in new[] { new EllipticCurveQ(0, 0, 1, 9, 9), new EllipticCurveQ(0, 1, 1, 6, 5) })
        {
            var result = curve.GetRankBounds(new RankComputationOptions { MaxDegreeOfParallelism = workers });
            Assert.Equal(0, result.LowerBound);
            Assert.Equal(2, result.UpperBound);
            Assert.Equal(2, result.TwoSelmerDimension);
            Assert.Null(result.ExactRank);
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    public void RepeatedOverlappingRegionsDoNotDuplicateCoveringClasses(int workers)
    {
        // These curves exercise three real resolvent roots and overlapping regions,
        // reducible coverings, and the distinction between rank and Selmer dimension.
        foreach (var (n, rank) in new[] { (1, 0), (5, 1), (34, 2) })
        for (int repeat = 0; repeat < 3; repeat++)
        {
            var result = new EllipticCurveQ(0, 0, 0, -n * n, 0).GetRankBounds(new RankComputationOptions
                { MaxDegreeOfParallelism = workers, PreferGeneralTwoDescent = true, MaxDescentWork = 10000000 });
            Assert.Equal(rank, result.ExactRank);
            Assert.Equal(rank + 2, result.TwoSelmerDimension);
        }
    }

    [Fact]
    public void WorkerSettingLeavesTheTwoIsogenyMethodAndItsBudgetUnchanged()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -25, 0);
        var sequential = curve.GetRankBounds();
        var configured = curve.GetRankBounds(new RankComputationOptions { MaxDegreeOfParallelism = 4 });
        Assert.True(configured.UsedTwoIsogenyDescent);
        Assert.Equal(sequential.ExactRank, configured.ExactRank);
        Assert.Equal(sequential.DescentWork, configured.DescentWork);
        Assert.Equal(sequential.PointSearchWork, configured.PointSearchWork);
    }

    [Fact]
    public async Task ConcurrentCallsOnTheSameCurveKeepProofsAndBudgetsSeparate()
    {
        var curve = new EllipticCurveQ(0, 1, 1, -2, 0);
        var sharedOptions = new RankComputationOptions { MaxDegreeOfParallelism = 2 };
        var results = await Task.WhenAll(Enumerable.Range(0, 6).Select(index => Task.Run(() =>
            index % 2 == 0 ? curve.GetRankBounds(sharedOptions)
                : curve.GetRankBounds(new RankComputationOptions { MaxDegreeOfParallelism = 2, MaxDescentWork = 1 }))));
        for (int index = 0; index < results.Length; index++)
        {
            Assert.Equal(2, results[index].LowerBound);
            if (index % 2 == 0) Assert.Equal(2, results[index].ExactRank);
            else { Assert.Null(results[index].UpperBound); Assert.Equal(1, results[index].DescentWork); }
        }
        Assert.Equal(5000000, sharedOptions.MaxDescentWork);
    }

    [Fact]
    public void WorkerBudgetViewsShareExactlyOneAllowanceUnderContention()
    {
        var budget = new DescentBudget(new RankComputationOptions
            { MaxDegreeOfParallelism = 4, MaxDescentWork = 100003, MaxPointSearchWork = 10007 }, default);
        long steps = 0, pointSteps = 0;
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            var worker = budget.WithToken(default);
            try { while (true) { worker.Step(); Interlocked.Increment(ref steps); } }
            catch (DescentLimitException) { }
            while (worker.PointStep()) Interlocked.Increment(ref pointSteps);
        });
        Assert.Equal(100003, steps);
        Assert.Equal(steps, budget.Work);
        Assert.Equal(10007, pointSteps);
        Assert.Equal(pointSteps, budget.PointWork);
        Assert.True(budget.PointSearchExhausted);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void WorkerFailureStopsAndJoinsPeersBeforeRethrowing(bool workLimit)
    {
        using var started = new CountdownEvent(2);
        var budget = new DescentBudget(new RankComputationOptions
            { MaxDegreeOfParallelism = 2, MaxDescentWork = long.MaxValue }, default);
        Exception failure = workLimit ? new DescentLimitException("test allowance exhausted") : new ArithmeticException("test failure");
        int active = 0;
        var observed = Record.Exception(() => DescentParallelism.ForEach(budget, _ => new[] { 0, 1 }, (item, worker) =>
        {
            Interlocked.Increment(ref active);
            try
            {
                started.Signal();
                Assert.True(started.Wait(TimeSpan.FromSeconds(15)), "Two bodies did not execute concurrently.");
                if (item == 0) throw failure;
                while (true) worker.Step();
            }
            finally { Interlocked.Decrement(ref active); }
        }));
        Assert.Same(failure, observed);
        Assert.Equal(0, active);
    }

    [Fact]
    public void SimultaneousArithmeticFailureIsNotHiddenByWorkExhaustion()
    {
        using var started = new CountdownEvent(2);
        var budget = new DescentBudget(new RankComputationOptions { MaxDegreeOfParallelism = 2 }, default);
        var failure = new ArithmeticException("must not become a partial rank");
        var observed = Record.Exception(() => DescentParallelism.ForEach(budget, _ => new[] { 0, 1 }, (item, worker) =>
        {
            started.Signal();
            Assert.True(started.Wait(TimeSpan.FromSeconds(15)));
            if (item == 0) throw new DescentLimitException("limit reached first");
            Assert.True(worker.Token.WaitHandle.WaitOne(TimeSpan.FromSeconds(15)));
            throw failure;
        }));
        Assert.Same(failure, observed);
    }

    [Fact]
    public void EnumerationFailureAlsoStopsAndJoinsActiveWorkers()
    {
        using var started = new ManualResetEventSlim();
        var budget = new DescentBudget(new RankComputationOptions
            { MaxDegreeOfParallelism = 2, MaxDescentWork = long.MaxValue }, default);
        var failure = new DescentLimitException("enumeration allowance exhausted");
        int active = 0;
        IEnumerable<int> Items(DescentBudget _)
        {
            yield return 0;
            Assert.True(started.Wait(TimeSpan.FromSeconds(15)));
            throw failure;
        }
        var observed = Record.Exception(() => DescentParallelism.ForEach(budget, Items, (_, worker) =>
        {
            Interlocked.Increment(ref active);
            try { started.Set(); while (true) worker.Step(); }
            finally { Interlocked.Decrement(ref active); }
        }));
        Assert.Same(failure, observed);
        Assert.Equal(0, active);
    }

    [Fact]
    public void CancellationDuringParallelWorkJoinsWorkersAndPreservesCallerToken()
    {
        using var cancellation = new CancellationTokenSource();
        using var started = new CountdownEvent(2);
        var budget = new DescentBudget(new RankComputationOptions
            { MaxDegreeOfParallelism = 2, MaxDescentWork = long.MaxValue }, cancellation.Token);
        int active = 0;
        var error = Assert.Throws<OperationCanceledException>(() => DescentParallelism.ForEach(budget, _ => new[] { 0, 1 }, (item, worker) =>
        {
            Interlocked.Increment(ref active);
            try
            {
                started.Signal();
                Assert.True(started.Wait(TimeSpan.FromSeconds(15)));
                if (item == 0) cancellation.Cancel();
                while (true) worker.Step();
            }
            finally { Interlocked.Decrement(ref active); }
        }));
        Assert.Equal(cancellation.Token, error.CancellationToken);
        Assert.Equal(0, active);
        Assert.Throws<OperationCanceledException>(() => new EllipticCurveQ(0, 0, 1, -1, 0)
            .GetRankBounds(new RankComputationOptions { MaxDegreeOfParallelism = 4 }, cancellation.Token));
    }
}

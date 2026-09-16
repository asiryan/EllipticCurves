using System.Numerics;
using Xunit;

namespace EllipticCurves.Tests;

public class ArithmeticCacheTests
{
    [Fact]
    public void MinimalModelAndItsPointMapsShareReadOnlyDiscriminantFactors()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0)
            .ChangeModel(new BigRational(2, 3), 5, -3, 7).Target;
        Assert.Equal(new BigInteger(32), curve.GetConductor(new FactorizationOptions(), out var conductorFactors));
        var minimal = curve.GetGlobalMinimalModel();
        var factors = curve.GetMinimalDiscriminantFactorization(default);

        Assert.Same(minimal, curve.GetGlobalMinimalModel());
        Assert.Same(minimal, minimal.GetGlobalMinimalModel());
        Assert.Same(factors, minimal.GetMinimalDiscriminantFactorization(default));
        for (int i = 0; i < 2; i++)
        {
            var map = curve.GetMinimalModelIsomorphism();
            Assert.Same(minimal, map.Target.GetGlobalMinimalModel());
            Assert.Same(factors, map.Target.GetMinimalDiscriminantFactorization(default));
            var point = map.MapBack(new EllipticCurvePoint(0, 0));
            Assert.True(curve.CanonicalHeight(point).Contains(0));
        }

        Assert.Equal(6, Assert.Single(factors).Value);
        Assert.Equal(5, Assert.Single(conductorFactors).Value);
        Assert.Equal(BigInteger.Abs(minimal.Discriminant.Num),
            factors.Aggregate(BigInteger.One, (n, factor) => n * BigInteger.Pow(factor.Key, factor.Value)));
        Assert.Throws<NotSupportedException>(() => ((IDictionary<BigInteger, int>)factors)[2] = 99);
        Assert.Equal(6, Assert.Single(curve.GetLocalData()).DiscriminantValuation);
        Assert.Equal(1, curve.GetRootNumber());

        var separate = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.NotSame(factors, separate.GetMinimalDiscriminantFactorization(default));
        var other = new EllipticCurveQ(0, 0, 1, -1, 0).GetMinimalDiscriminantFactorization(default);
        Assert.Equal(new BigInteger(37), Assert.Single(other).Key);
        Assert.Equal(1, Assert.Single(other).Value);
    }

    [Fact]
    public void WarmCurveCachesStillHonorCancellationAndValidateWorkerOptions()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        curve.GetMinimalDiscriminantFactorization(default);
        var token = new CancellationToken(true);
        Assert.Throws<OperationCanceledException>(() => curve.GetGlobalMinimalModel(token));
        Assert.Throws<OperationCanceledException>(() => curve.GetMinimalModelIsomorphism(token));
        Assert.Throws<OperationCanceledException>(() => curve.GetMinimalDiscriminantFactorization(token));
        Assert.Throws<OperationCanceledException>(() => curve.GetConductor(token));
        Assert.Throws<OperationCanceledException>(() => curve.GetRootNumber(token));
        Assert.Throws<OperationCanceledException>(() => curve.GetLocalData(token));
        Assert.Throws<OperationCanceledException>(() => curve.CanonicalHeight(new EllipticCurvePoint(0, 0), cancellationToken: token));
        Assert.Throws<ArgumentOutOfRangeException>(() => curve.GetConductor(
            new FactorizationOptions { MaxDegreeOfParallelism = -1 }, default));
    }

    [Fact]
    public void FailedOrCancelledPreparationCanBeRetriedAndSuccessIsReused()
    {
        var cache = new CachedComputation<object>();
        Assert.Throws<ArithmeticException>(() => cache.Get(() => throw new ArithmeticException(), default));
        using var cancelled = new CancellationTokenSource();
        Assert.Throws<OperationCanceledException>(() => cache.Get(() =>
        {
            cancelled.Cancel();
            return new object();
        }, cancelled.Token));

        var expected = new object();
        Assert.Same(expected, cache.Get(() => expected, default));
        Assert.Same(expected, cache.Get(() => throw new InvalidOperationException("Already computed"), default));
        Assert.Throws<OperationCanceledException>(() => cache.Get(() => expected, cancelled.Token));
    }

    [Fact]
    public async Task ConcurrentCallersShareOnePreparationAndAWaitingCallerCanCancel()
    {
        var cache = new CachedComputation<object>();
        var expected = new object();
        var entered = new TaskCompletionSource<bool>(TaskCreationOptions.RunContinuationsAsynchronously);
        using var release = new ManualResetEventSlim();
        using var cancelled = new CancellationTokenSource();
        int preparations = 0;
        var owner = Task.Run(() => cache.Get(() =>
        {
            Interlocked.Increment(ref preparations);
            entered.SetResult(true);
            Assert.True(release.Wait(TimeSpan.FromSeconds(30)));
            return expected;
        }, default));
        var waiters = new List<Task<object>>();
        try
        {
            await entered.Task.WaitAsync(TimeSpan.FromSeconds(10));
            // The owner stays blocked throughout cancellation of this caller.
            cancelled.CancelAfter(TimeSpan.FromMilliseconds(50));
            var error = Assert.ThrowsAny<OperationCanceledException>(() => cache.Get(() =>
                throw new InvalidOperationException("A waiting caller must not prepare a second value"), cancelled.Token));
            Assert.Equal(cancelled.Token, error.CancellationToken);
            Assert.False(owner.IsCompleted);
            for (int i = 0; i < 4; i++)
                waiters.Add(Task.Run(() => cache.Get(() =>
                {
                    Interlocked.Increment(ref preparations);
                    return new object();
                }, default)));
        }
        finally { release.Set(); }
        Assert.Same(expected, await owner.WaitAsync(TimeSpan.FromSeconds(10)));
        foreach (var result in await Task.WhenAll(waiters).WaitAsync(TimeSpan.FromSeconds(10))) Assert.Same(expected, result);
        Assert.Equal(1, preparations);
    }

    [Fact]
    public async Task WaitingCallerRetriesAfterThePreparingCallerIsCancelled()
    {
        var cache = new CachedComputation<object>();
        var expected = new object();
        var entered = new TaskCompletionSource<bool>(TaskCreationOptions.RunContinuationsAsynchronously);
        using var release = new ManualResetEventSlim();
        using var cancelled = new CancellationTokenSource();
        var owner = Task.Run(() => cache.Get(() =>
        {
            entered.SetResult(true);
            Assert.True(release.Wait(TimeSpan.FromSeconds(30)));
            return new object();
        }, cancelled.Token));
        Task<object> waiter;
        try
        {
            await entered.Task.WaitAsync(TimeSpan.FromSeconds(10));
            waiter = Task.Run(() => cache.Get(() => expected, default));
            cancelled.Cancel();
        }
        finally { release.Set(); }
        await Assert.ThrowsAnyAsync<OperationCanceledException>(async () => await owner.WaitAsync(TimeSpan.FromSeconds(10)));
        Assert.Same(expected, await waiter.WaitAsync(TimeSpan.FromSeconds(10)));
        Assert.Same(expected, cache.Get(() => throw new InvalidOperationException("Already computed"), default));
    }
}

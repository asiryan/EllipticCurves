using System.Text.Json;
using System.Diagnostics;
using System.Numerics;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public class ElkiesSearchTests
{
    [Theory]
    [InlineData(0, 1, 15)]
    [InlineData(1, 1, 15)]
    [InlineData(-1, 1, 15)]
    [InlineData(2, 3, 17)]
    [InlineData(-9529, 5471, 17)]
    public void PublishedSectionsGiveExactPointsAndHonestLowerBounds(int a, int b, int expectedLower)
    {
        var (curve, points) = ElkiesSearchFamily.Create(a, b);
        Assert.False(curve.IsSingular);
        Assert.Equal(17, points.Length);
        Assert.All(points, point => Assert.True(curve.IsOnCurve(point)));
        var certificate = curve.GetRankLowerBound(points);
        Assert.Equal(expectedLower, certificate.LowerBound);
        Assert.Equal(expectedLower == 17, certificate.IndependenceCertified);
        // The 28-rank specialization only supplies the 17 family sections here.
        Assert.Equal(0, certificate.TwoTorsionDimensionUpperBound);
    }

    [Fact]
    public void CertificateHandlesTorsionDependenceAndRationalModels()
    {
        var torsion = new EllipticCurveQ(0, 0, 0, -1, 0);
        Assert.Equal(0, torsion.GetRankLowerBound(new[] { new EllipticCurvePoint(0, 0), new(1, 0), new(-1, 0) }).LowerBound);
        var (curve, points) = ElkiesSearchFamily.Create(2, 3);
        var map = curve.ChangeModel(6, new BigRational(1, 7), new BigRational(1, 5), new BigRational(1, 3));
        var mapped = points.Select(map.Map).ToArray();
        Assert.Equal(17, map.Target.GetRankLowerBound(mapped).LowerBound);
        points[16] = curve.Add(points[0], points[1]);
        var dependent = curve.GetRankLowerBound(points);
        Assert.Equal(16, dependent.LowerBound);
        Assert.False(dependent.IndependenceCertified);
        points[16] = new(points[0].X, points[0].Y + 1);
        Assert.Throws<ArgumentException>(() => curve.GetRankLowerBound(points));
        Assert.Throws<OperationCanceledException>(() => curve.GetRankLowerBound(Array.Empty<EllipticCurvePoint>(), cancellationToken: new(true)));
    }

    [Fact]
    public void SearchCanResumeAtAnExactCheckpointAndRejectTamperedPoints()
    {
        var options = new CurveSearchOptions(-20, 20, 8, 31, 1009, 2, 0, 1);
        var checkpoints = new List<CurveSearchState>();
        var completed = CurveSearchEngine.Run(new(options), update => checkpoints.Add(CurveSearchState.Parse(update.Result!)));
        Assert.True(completed.Complete);
        Assert.Equal(Enumerable.Range(1, 8).Sum(b => Enumerable.Range(-20, 41).Count(a => System.Numerics.BigInteger.GreatestCommonDivisor(a, b) == 1)), completed.Tested);
        var checkpoint = checkpoints.First(s => s.NextSlot > 0 && !s.Complete);
        var fake = checkpoint.Results.Select(c => c with { LowerBound = 0, PointsText = "invalid" }).ToArray();
        var resumed = CurveSearchEngine.Run(checkpoint with { Candidates = fake, Options = options with { Workers = 2 } });
        Assert.Equal(completed.Results, resumed.Results);
        Assert.Equal(completed.Tested, resumed.Tested);
        Assert.All(resumed.Results, c => Assert.InRange(c.LowerBound, 1, 17));
        Assert.Equal(completed, CurveSearchState.Parse(JsonSerializer.Serialize(completed)) with { Candidates = completed.Candidates });
        Assert.Throws<System.IO.InvalidDataException>(() => (completed with { NextSlot = options.Slots + 1 }).Validate());
    }

    [Fact]
    public void ResidueTableScoresMatchDirectCountsIncludingDenominatorsDivisibleByPrimes()
    {
        var state = CurveSearchEngine.Run(new(new(1, 1, 7, 31, 101, 7, 0, 2)));
        Assert.Equal(7, state.Results.Length);
        foreach (var candidate in state.Results)
        {
            var (curve, _) = ElkiesSearchFamily.Create(candidate.Numerator, candidate.Denominator);
            double expected = 0;
            foreach (int p in new[] { 5, 7, 11, 13, 17, 19, 23, 29, 31 })
            {
                if (curve.Discriminant.Num % p == 0) continue;
                int count = 1;
                for (int x = 0; x < p; x++)
                {
                    var rhs = ((BigInteger)x * x * x + curve.A4.Num * x + curve.A6.Num) % p;
                    rhs = (rhs + p) % p;
                    count += rhs.IsZero ? 1 : BigInteger.ModPow(rhs, (p - 1) / 2, p) == 1 ? 2 : 0;
                }
                expected += Math.Log((double)count / p);
            }
            Assert.Equal(expected, candidate.Score, 12);
        }
    }

    private static CalculationRunner Worker() => new(() =>
    {
        var start = new ProcessStartInfo("dotnet");
        start.ArgumentList.Add(Path.Combine(AppContext.BaseDirectory, "worker", "WorkerHost.dll"));
        return start;
    });

    [Fact]
    public async Task DefaultSearchRunsInWorkerAndSavedReportsAreRechecked()
    {
        using var workbench = new WorkbenchViewModel();
        using var model = new CurveSearchViewModel(workbench, Worker());
        Assert.True(model.CanStart);
        var run = model.StartAsync();
        Assert.True(model.IsBusy);
        Assert.False(workbench.CanRun);
        await run.WaitAsync(TimeSpan.FromSeconds(30));
        Assert.True(model.Snapshot.Complete, model.Error);
        Assert.Equal("", model.Error);
        Assert.True(workbench.CanRun);
        Assert.NotEmpty(model.Candidates);
        var saved = model.Snapshot;
        string path = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ecsearch");
        try
        {
            saved.Save(path);
            model.Load(CurveSearchState.Load(path));
            Assert.Contains("SAVED REPORT", model.Details);
            await model.StartAsync().WaitAsync(TimeSpan.FromSeconds(30));
            Assert.DoesNotContain("SAVED REPORT", model.Details);
            Assert.Equal(saved.Results, model.Snapshot.Results);
        }
        finally { File.Delete(path); }
        model.NewSearch();
        model.Settings[0].Text = "not a number";
        Assert.False(model.CanStart);
        Assert.NotEmpty(model.InputError);
    }

    [Fact]
    public async Task PauseKillsWorkerAndRetainsCompletedCheckpoint()
    {
        using var workbench = new WorkbenchViewModel();
        using var model = new CurveSearchViewModel(workbench, Worker());
        model.Settings[0].Text = "-1000";
        model.Settings[1].Text = "1000";
        model.Settings[2].Text = "1000";
        bool paused = false;
        model.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName != nameof(model.Percent) || !model.HasProgress || paused) return;
            paused = true;
            workbench.Cancel();
        };
        await model.StartAsync().WaitAsync(TimeSpan.FromSeconds(30));
        Assert.True(paused);
        Assert.True(model.HasProgress);
        Assert.False(model.Snapshot.Complete);
        Assert.True(model.CanStart);
        Assert.True(workbench.CanRun);
        Assert.Contains("Paused", model.Status);
        Assert.NotEmpty(model.Candidates);
        Assert.False(model.CanEdit);
    }
}

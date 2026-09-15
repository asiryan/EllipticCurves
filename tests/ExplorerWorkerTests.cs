using System.Collections.Concurrent;
using System.Diagnostics;
using System.Text.Json;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerWorkerTests
{
    private static ProcessStartInfo StartInfo(string extra = null)
    {
        var start = new ProcessStartInfo("dotnet");
        start.ArgumentList.Add(Path.Combine(AppContext.BaseDirectory, "worker", "WorkerHost.dll"));
        if (extra != null) start.ArgumentList.Add(extra);
        return start;
    }
    private static CalculationRequest Request(string id = "Q.TorsionStructure") => new(id, "y^2 = x^3 - x", new());

    [Fact]
    public async Task WorkerProtocolReturnsProgressAndFinalResult()
    {
        var progress = new ConcurrentQueue<CalculationUpdate>();
        var outcome = await new CalculationRunner(() => StartInfo()).RunAsync(Request(), progress.Enqueue);
        Assert.Equal("Completed", outcome.Status);
        Assert.Contains("2", outcome.Text);
        Assert.Contains(progress, p => p.Kind == "progress");
        Assert.Contains(progress, p => p.Kind == "completed");
    }

    [Fact]
    public async Task ParallelRankRunsThroughTheRealWorkerProtocol()
    {
        var operation = CalculationCatalog.All.Single(o => o.Member?.Name == nameof(EllipticCurveQ.GetRankBounds)
            && o.Parameters.Any(p => p.Key == "options.MaxDegreeOfParallelism"));
        var request = new CalculationRequest(operation.Id, "y^2 = x^3 - 49/50*x + 1/2",
            new() { ["options.MaxDescentWork"] = "20000000" });
        var outcome = await new CalculationRunner(() => StartInfo()).RunAsync(request, _ => { })
            .WaitAsync(TimeSpan.FromSeconds(60));
        Assert.Equal("Completed", outcome.Status);
        Assert.Contains("Exact Rank: 1", outcome.Text);
        Assert.Contains("Two Selmer Dimension: 1", outcome.Text);
    }

    [Fact]
    public async Task WorkerPropagatesLibraryErrors()
    {
        var op = CalculationCatalog.All.Single(o => o.Member?.Name == "GetConductor");
        var request = Request(op.Id) with { Equation = "y^2 = x^3" };
        var result = await new CalculationRunner(() => StartInfo()).RunAsync(request, _ => { });
        Assert.Equal("Failed", result.Status);
        Assert.Contains("singular", result.Message.ToLowerInvariant());
    }

    [Fact]
    public async Task LargeConductorCompletesThroughTheRealWorkerProtocol()
    {
        var operation = CalculationCatalog.All.Single(o => o.Member?.Name == "GetConductor");
        var request = new CalculationRequest(operation.Id,
            "y^2 = x^3 + x^2 - 221556180740323405132844117936*x + 35386140191724122461245294467670188433973860",
            new(), TimeoutSeconds: 15);
        var result = await new CalculationRunner(() => StartInfo()).RunAsync(request, _ => { })
            .WaitAsync(TimeSpan.FromSeconds(30));
        Assert.Equal("Completed", result.Status);
        Assert.Contains("1103561624055499058867562340698878392772504928025988266715523317532246643920", result.Text);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public async Task StopAndTimeoutKillNonCooperativeProcess(bool timeout)
    {
        var ready = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        using var cancellation = new CancellationTokenSource();
        var request = Request() with { TimeoutSeconds = timeout ? 2 : 0 };
        var watch = Stopwatch.StartNew();
        var task = new CalculationRunner(() => StartInfo("--unresponsive")).RunAsync(request, _ => ready.TrySetResult(), cancellation.Token);
        await ready.Task.WaitAsync(TimeSpan.FromSeconds(10));
        if (!timeout) cancellation.Cancel();
        var outcome = await task.WaitAsync(TimeSpan.FromSeconds(10));
        Assert.Equal(timeout ? "Timed out" : "Cancelled", outcome.Status);
        Assert.True(watch.Elapsed < TimeSpan.FromSeconds(12));
    }

    [Fact]
    public async Task ClosingHostInputDoesNotOrphanWorker()
    {
        var start = StartInfo();
        start.UseShellExecute = false;
        start.CreateNoWindow = true;
        start.RedirectStandardInput = start.RedirectStandardOutput = true;
        using var process = Process.Start(start);
        try
        {
            // Huge integer scalar multiplication on a nontorsion point keeps the worker busy.
            var multiply = CalculationCatalog.All.Single(o => o.Member?.Name == "Multiply" && o.Context == CalculationContext.RationalCurve);
            var request = new CalculationRequest(multiply.Id, "y^2 + y = x^3 - x",
                new() { ["P.infinity"] = "False", ["P.x"] = "0", ["P.y"] = "0", ["n"] = new string('9', 1000) });
            await process.StandardInput.WriteLineAsync(JsonSerializer.Serialize(request));
            await process.StandardInput.FlushAsync();
            Assert.NotNull(await process.StandardOutput.ReadLineAsync().WaitAsync(TimeSpan.FromSeconds(10)));
            process.StandardInput.Close();
            await process.WaitForExitAsync().WaitAsync(TimeSpan.FromSeconds(10));
            Assert.Equal(2, process.ExitCode);
        }
        finally { if (!process.HasExited) process.Kill(true); }
    }

    [Fact]
    public void DeleteMovesSelectionAndClearsTheLastResult()
    {
        using var workbench = new WorkbenchViewModel();
        var newest = new CalculationJobViewModel(Request(), "Newest") { Status = "Completed" };
        var middle = new CalculationJobViewModel(Request(), "Middle") { Status = "Failed" };
        var oldest = new CalculationJobViewModel(Request(), "Oldest") { Status = "Cancelled" };
        foreach (var job in new[] { newest, middle, oldest }) workbench.Jobs.Add(job);
        workbench.Selected = middle;
        var notifications = new List<string>();
        workbench.PropertyChanged += (_, args) => notifications.Add(args.PropertyName);

        workbench.Delete(middle);
        Assert.Equal(new[] { newest, oldest }, workbench.Jobs);
        Assert.Same(oldest, workbench.Selected);
        Assert.Contains(nameof(workbench.Summary), notifications);
        workbench.Delete(oldest);
        Assert.Same(newest, workbench.Selected);
        workbench.Delete(newest);
        Assert.Empty(workbench.Jobs);
        Assert.Null(workbench.Selected);
        Assert.False(workbench.HasSelection);
        Assert.False(workbench.HasResults);
        Assert.False(workbench.CanDelete(workbench.Selected));
        Assert.Contains(nameof(workbench.HasResults), notifications);
        Assert.Equal("Choose a calculation in Tools", workbench.Summary);
        workbench.Delete(null); // An empty history is a harmless no-op.
    }

    [Fact]
    public async Task ClearHistoryRemovesAllResultsAndAllowsNewCalculations()
    {
        using var workbench = new WorkbenchViewModel(new CalculationRunner(() => StartInfo()));
        Assert.False(workbench.CanClearHistory);
        var jobs = new[] { "Completed", "Failed", "Cancelled", "Timed out" }
            .Select(status => new CalculationJobViewModel(Request(), status) { Status = status }).ToArray();
        foreach (var job in jobs) workbench.Jobs.Add(job);
        workbench.Selected = jobs[1];
        var notifications = new List<string>();
        workbench.PropertyChanged += (_, args) => notifications.Add(args.PropertyName);
        Assert.True(workbench.CanClearHistory);

        workbench.ClearHistory();
        Assert.Empty(workbench.Jobs);
        Assert.Null(workbench.Selected);
        Assert.False(workbench.HasResults);
        Assert.False(workbench.HasSelection);
        Assert.False(workbench.CanClearHistory);
        Assert.Equal("Choose a calculation in Tools", workbench.Summary);
        Assert.Contains(nameof(workbench.CanClearHistory), notifications);
        Assert.Contains(nameof(workbench.HasResults), notifications);
        Assert.Contains(nameof(workbench.HasSelection), notifications);
        workbench.ClearHistory(); // Clearing an empty history is harmless.

        await workbench.RunAsync(Request()).WaitAsync(TimeSpan.FromSeconds(10));
        Assert.Equal("Completed", Assert.Single(workbench.Jobs).Status);
        Assert.True(workbench.CanClearHistory);
    }

    [Fact]
    public async Task ClearHistoryDoesNotRemoveOrCancelAnActiveCalculation()
    {
        using var workbench = new WorkbenchViewModel(new CalculationRunner(() => StartInfo("--unresponsive")));
        var older = new CalculationJobViewModel(Request(), "Older") { Status = "Completed" };
        workbench.Jobs.Add(older);
        var calculation = workbench.RunAsync(Request());
        var active = workbench.Active;
        try
        {
            Assert.False(workbench.CanClearHistory);
            workbench.ClearHistory();
            Assert.Equal(new[] { active, older }, workbench.Jobs);
            Assert.Same(active, workbench.Selected);
            Assert.True(workbench.IsBusy);
            Assert.Equal("Running", active.Status);
        }
        finally
        {
            workbench.Cancel();
            await calculation.WaitAsync(TimeSpan.FromSeconds(10));
        }
        Assert.True(workbench.CanClearHistory);
        workbench.ClearHistory();
        Assert.Empty(workbench.Jobs);
        Assert.Null(workbench.Selected);
    }

    [Fact]
    public async Task DeleteProtectsTheActiveCalculationButAllowsOlderResults()
    {
        using var workbench = new WorkbenchViewModel(new CalculationRunner(() => StartInfo("--unresponsive")));
        var older = new CalculationJobViewModel(Request(), "Older") { Status = "Completed" };
        workbench.Jobs.Add(older);
        var calculation = workbench.RunAsync(Request());
        var active = workbench.Active;
        try
        {
            Assert.False(workbench.CanDelete(active));
            workbench.Delete(active);
            Assert.Contains(active, workbench.Jobs);
            Assert.True(workbench.CanDelete(older));
            workbench.Delete(older);
            Assert.Same(active, workbench.Selected);
            Assert.Same(active, Assert.Single(workbench.Jobs));
            Assert.True(workbench.IsBusy);
        }
        finally
        {
            workbench.Cancel();
            await calculation.WaitAsync(TimeSpan.FromSeconds(10));
        }
        Assert.True(workbench.CanDelete(active));
        workbench.Delete(active);
        Assert.Empty(workbench.Jobs);
    }

    [Fact]
    public async Task WorkbenchKeepsCapturedInputsAndCanRunAgainAfterCancellation()
    {
        var calls = 0;
        var runner = new CalculationRunner(() => StartInfo(++calls == 1 ? "--unresponsive" : null));
        using var workbench = new WorkbenchViewModel(runner);
        var request = Request();
        var first = workbench.RunAsync(request);
        Assert.True(workbench.IsBusy);
        Assert.False(workbench.CanRun);
        workbench.Cancel();
        await first.WaitAsync(TimeSpan.FromSeconds(10));
        Assert.Equal("Cancelled", workbench.Selected.Status);
        await workbench.RunAsync(request with { Equation = "y^2 = x^3 + 1" });
        Assert.Equal(2, workbench.Jobs.Count);
        Assert.Equal("y^2 = x^3 - x", workbench.Jobs[1].Equation);
        Assert.Equal("Completed", workbench.Selected.Status);
        Assert.True(workbench.CanRun);
        Assert.Contains("Captured plot: y^2 = x^3 + 1", workbench.Selected.Report);
    }
}

using System.Collections.Concurrent;
using System.Diagnostics;
using System.Text.Json;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class VisualizerWorkerTests
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
    public async Task WorkerPropagatesLibraryErrors()
    {
        var op = CalculationCatalog.All.Single(o => o.Member?.Name == "GetConductor");
        var request = Request(op.Id) with { Equation = "y^2 = x^3" };
        var result = await new CalculationRunner(() => StartInfo()).RunAsync(request, _ => { });
        Assert.Equal("Failed", result.Status);
        Assert.Contains("singular", result.Message.ToLowerInvariant());
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

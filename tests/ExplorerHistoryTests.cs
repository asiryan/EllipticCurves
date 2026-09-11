using System.Diagnostics;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerHistoryTests
{
    [Fact]
    public void HistoryBranchesAndEvictsOnlyTheOldestUndoSteps()
    {
        var history = new MementoHistory<int>(2);
        history.Record(0);
        history.Record(1);
        history.Record(2);
        Assert.Equal(2, history.Undo(3));
        Assert.Equal(1, history.Undo(2));
        Assert.False(history.CanUndo);
        Assert.Equal(2, history.Redo(1));
        history.Record(2);
        Assert.False(history.CanRedo);
        Assert.Equal(2, history.Undo(4));
        history.Clear();
        Assert.False(history.CanUndo);
        Assert.False(history.CanRedo);
    }

    [Fact]
    public async Task EditorMementoRestoresExactInputSliderAnchorsAndComputedSamples()
    {
        using var model = new MainViewModel();
        model.Step.Text = "1/7";
        model.SimpleCoefficients[0].Text = "-2/7";
        model.SimpleCoefficients[0].SliderOffset = 3;
        model.FlushUpdate();
        var state = model.CaptureMemento();
        await model.PendingSamples; // The result enriches the existing memento.
        var points = model.Samples;
        var slider = model.SimpleCoefficients[0].CaptureMemento();
        model.ApplyPreset(CurvePreset.Classic);
        model.Equation.Text = "y^2 + xy + y = x^3 - x";
        model.FlushUpdate();
        model.RestoreMemento(state);
        Assert.Same(state.Snapshot, model.Snapshot);
        Assert.Same(points, model.Samples);
        Assert.Equal("1/7", model.Step.Text);
        Assert.Equal(slider, model.SimpleCoefficients[0].CaptureMemento());
        Assert.True(model.PendingSamples.IsCompletedSuccessfully);
        Assert.True(model.PendingUpdate.IsCompletedSuccessfully);
        await Task.Delay(350);
        Assert.Same(state.Snapshot, model.Snapshot); // No stale debounce can overwrite Undo.
        model.SimpleCoefficients[0].SliderOffset = 4;
        Assert.Equal(new BigRational(2, 7), model.SimpleCoefficients[0].ExactValue);
    }

    [Fact]
    public void IncompleteInputCanBeRestoredWithoutParsingOrLosingTheLastValidCurve()
    {
        using var model = new MainViewModel();
        model.Equation.Text = "y^2 =";
        model.Equation.CommitEdit();
        var state = model.CaptureMemento();
        model.ApplyPreset(CurvePreset.Classic);
        model.RestoreMemento(state);
        Assert.Equal("y^2 =", model.Equation.Text);
        Assert.True(model.HasInputError);
        Assert.True(model.HasIncompleteInput);
        Assert.Same(state.Snapshot, model.Snapshot);
    }

    private static CalculationSession Report(int index) => new(
        new("Q.TorsionStructure", "y^2 = x^3 - x", new() { ["input"] = index.ToString() }),
        new DateTime(2026, 1, 1).AddMinutes(index), CalculationStatus.Completed, "Done", TimeSpan.FromSeconds(index), 100, "Report " + index);

    [Fact]
    public void DeleteAndClearRestoreDetachedReportsOrderAndSelection()
    {
        using var workbench = new WorkbenchViewModel();
        workbench.RestoreHistory(new[] { Report(3), Report(2), Report(1) });
        workbench.Selected = workbench.Jobs[1];
        var state = workbench.CaptureMemento();
        var report = workbench.Selected.Report;
        workbench.Selected.Request.Arguments["input"] = "changed externally";
        workbench.Delete(workbench.Selected);
        workbench.ClearHistory();
        workbench.RestoreMemento(state);
        Assert.Equal(new[] { "Report 3", "Report 2", "Report 1" }, workbench.Jobs.Select(job => job.Result));
        Assert.Same(workbench.Jobs[1], workbench.Selected);
        Assert.Equal(report, workbench.Selected.Report);
        Assert.Equal("2", workbench.Selected.Request.Arguments["input"]);
        Assert.False(workbench.IsBusy);
        Assert.Null(workbench.Active);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public async Task RunningMementosReceiveFinalResultsAndRestoreWithoutStartingAWorker(bool cancel)
    {
        var launches = 0;
        var runner = new CalculationRunner(() =>
        {
            launches++;
            var start = new ProcessStartInfo("dotnet");
            start.ArgumentList.Add(Path.Combine(AppContext.BaseDirectory, "worker", "WorkerHost.dll"));
            if (cancel) start.ArgumentList.Add("--unresponsive");
            return start;
        });
        using var workbench = new WorkbenchViewModel(runner);
        workbench.RestoreHistory(Enumerable.Range(0, ExplorerSession.HistoryLimit).Select(Report).ToArray());
        var before = workbench.CaptureMemento();
        var running = workbench.RunAsync(new("Q.TorsionStructure", "y^2 = x^3 - x", new()));
        var during = workbench.CaptureMemento();
        Assert.Equal(ExplorerSession.HistoryLimit, workbench.Jobs.Count);
        Assert.Throws<InvalidOperationException>(() => workbench.RestoreMemento(before));
        // Deleting an older report while running must remain independently undoable.
        workbench.Delete(workbench.Jobs[1]);
        if (cancel) workbench.Cancel();
        await running.WaitAsync(TimeSpan.FromSeconds(15));
        var report = workbench.Jobs[0].Report;
        Assert.Equal(cancel ? CalculationStatus.Cancelled : CalculationStatus.Completed, workbench.Jobs[0].Status);
        workbench.RestoreMemento(before);
        Assert.Equal("Report 49", workbench.Jobs[^1].Result);
        workbench.RestoreMemento(during);
        Assert.Equal(report, workbench.Jobs[0].Report);
        Assert.Equal(ExplorerSession.HistoryLimit, workbench.Jobs.Count);
        Assert.False(workbench.IsBusy);
        Assert.Equal(1, launches);
    }

    [Fact]
    public async Task TorusMementoRetainsPreparedDataWithoutRecomputingPeriods()
    {
        var calls = 0;
        using var torus = new ComplexTorusViewModel((curve, token) => { calls++; return TorusLattice.Create(curve, token); });
        var curve = CurvePreset.Classic.CreateCurve();
        var samples = Array.Empty<EllipticCurvePoint>();
        torus.Update(curve, samples, true);
        var state = torus.CaptureMemento();
        await torus.PendingUpdate;
        Assert.NotNull(torus.Lattice);
        var lattice = torus.Lattice;
        torus.Update(new EllipticCurveQ(0, 0, 0, 0, 7), samples, false);
        torus.RestoreMemento(state, curve, samples, true);
        torus.Update(curve, samples, true);
        Assert.Same(lattice, torus.Lattice);
        Assert.False(torus.IsBusy);
        Assert.True(torus.PendingUpdate.IsCompletedSuccessfully);
        Assert.Equal(1, calls);
    }

    [Fact]
    public async Task RestoringAnUnchangedGraphKeepsItsOriginalBackgroundWork()
    {
        using var model = new MainViewModel();
        var samples = model.PendingSamples;
        model.RestoreMemento(model.CaptureMemento());
        Assert.Same(samples, model.PendingSamples);
        await samples;
        Assert.NotEmpty(model.Samples);

        var calls = 0;
        using var torus = new ComplexTorusViewModel((curve, token) => { calls++; return TorusLattice.Create(curve, token); });
        torus.Update(model.Snapshot.Curve, model.Samples, true);
        var pending = torus.PendingUpdate;
        torus.RestoreMemento(torus.CaptureMemento(), model.Snapshot.Curve, model.Samples, true);
        Assert.Same(pending, torus.PendingUpdate);
        await pending;
        Assert.NotNull(torus.Lattice);
        Assert.Equal(1, calls);
    }
}

#nullable enable
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class CurveSearchSetting(string label, string help, int value) : ObservableObject
{
    private string text = value.ToString(CultureInfo.InvariantCulture);
    public string Label { get; } = label;
    public string Help { get; } = help;
    public string Text { get => text; set { text = value; OnPropertyChanged(); } }
}

public sealed class CurveSearchViewModel : ObservableObject, IDisposable
{
    private readonly WorkbenchViewModel workbench;
    private readonly CalculationRunner runner;
    private CancellationTokenSource? running;
    private CurveSearchState? state;
    private CurveSearchCandidate? selected;
    private bool disposed, verified = true, dirty;
    private string status = "Ready. The family formulas and 17 initial points are built in.", error = "";
    public CurveSearchViewModel(WorkbenchViewModel workbench, CalculationRunner? runner = null)
    {
        this.workbench = workbench;
        this.runner = runner ?? new();
        foreach (var field in Settings) field.PropertyChanged += SettingChanged;
        workbench.PropertyChanged += WorkbenchChanged;
    }

    public CurveSearchSetting[] Settings { get; } =
    {
        new("Numerator from", "Smallest a in t = a/b.", -50),
        new("Numerator to", "Largest a in t = a/b.", 50),
        new("Denominator up to", "Positive b; duplicate fractions are skipped.", 20),
        new("CPU workers", "Concurrent workers, from 1 to 32. Start with 2.", 2),
        new("Score primes up to", "Local point counts rank candidates. A higher score is only a heuristic.", 127),
        new("Certificate primes up to", "More primes can improve the proved rank lower bound. Maximum 10000.", 1009),
        new("Keep top candidates", "Retain the best 1–32 candidates by heuristic score.", 8),
        new("Extra point search depth", "Try x = section x + k/d², |k| up to this value, d = 1…4. 0 disables this small search; 2 seconds per candidate.", 16)
    };
    public IEnumerable<CurveSearchSetting> BasicSettings => Settings.Take(4);
    public IEnumerable<CurveSearchSetting> AdvancedSettings => Settings.Skip(4);
    public ObservableCollection<CurveSearchCandidate> Candidates { get; } = new();
    public CurveSearchCandidate? Selected
    {
        get => selected;
        set { selected = value; OnPropertyChanged(); OnPropertyChanged(nameof(Details)); OnPropertyChanged(nameof(CanOpenCurve)); }
    }
    public string Details => Selected == null ? "Candidates will appear here after the first batch.\n\nEach curve starts with the 17 published family sections. The program checks their coordinates and certifies a rank lower bound. A small additional point search is also included."
        : (verified ? "" : "SAVED REPORT · Resume / recheck before treating these values as verified in this run.\n\n") + Selected.Details;
    public bool IsBusy => running != null;
    public bool IsIdle => !IsBusy;
    public bool CanEdit => !IsBusy && state == null;
    public bool CanStart => !disposed && !IsBusy && workbench.CanRun && InputError.Length == 0 && (state == null || !state.Complete || !verified);
    public bool CanOpenCurve => !IsBusy && workbench.CanRun && Selected != null;
    public bool IsDirty => dirty;
    public bool HasProgress => state?.NextSlot > 0;
    public string StartLabel => state == null ? "Start search" : state.Complete ? "Recheck" : "Resume";
    public string InputError { get { try { ReadOptions(); return ""; } catch (Exception e) when (e is ArgumentException or FormatException or OverflowException) { return e.Message; } } }
    public string Status => status;
    public string Error => error;
    public double Percent => state == null ? 0 : 100.0 * state.NextSlot / state.Options.Slots;
    public string ProgressText => state == null ? "t = a/b · duplicate fractions are skipped · no data files needed"
        : $"{state.Tested:N0} distinct parameters · {Percent:F1}% of the range · {Candidates.Count} retained by score";
    public string PersistenceText => dirty ? "Unsaved search · Save search keeps a checkpoint for the next Explorer session." : "Search files (.ecsearch) are separate from plot sessions (.ec).";
    public CurveSearchState Snapshot => state ?? new(ReadOptions());

    private CurveSearchOptions ReadOptions()
    {
        int[] values;
        try { values = Settings.Select(s => int.Parse(s.Text, NumberStyles.Integer, CultureInfo.InvariantCulture)).ToArray(); }
        catch (Exception e) when (e is FormatException or OverflowException) { throw new FormatException("Enter whole numbers in all search settings."); }
        var options = new CurveSearchOptions(values[0], values[1], values[2], values[4], values[5], values[6], values[7], values[3]);
        options.Validate();
        return options;
    }

    public async Task StartAsync()
    {
        if (!CanStart) return;
        state ??= new(ReadOptions());
        using var cancellation = new CancellationTokenSource();
        running = cancellation;
        workbench.BeginExternal(Pause);
        error = "";
        status = "Starting search worker…";
        Notify();
        var context = SynchronizationContext.Current;
        void Apply(CalculationUpdate update)
        {
            if (disposed || !ReferenceEquals(running, cancellation) || update.Kind != CalculationProtocol.Progress) return;
            try
            {
                if (update.Result != null) SetState(CurveSearchState.Parse(update.Result));
                if (!update.Message.StartsWith("Preparing", StringComparison.Ordinal)) verified = true;
                if (!cancellation.IsCancellationRequested) status = update.Message;
                Notify();
            }
            catch (Exception e) when (e is System.IO.IOException or JsonException or ArgumentException)
            { error = e.Message; Pause(); }
        }
        try
        {
            var request = new CalculationRequest(CurveSearchEngine.OperationId, "", new() { ["state"] = JsonSerializer.Serialize(state) }, TimeoutSeconds: 0);
            var outcome = await runner.RunAsync(request, update =>
            {
                if (context == null) Apply(update);
                else context.Post(_ => Apply(update), null);
            }, cancellation.Token);
            if (disposed) return;
            if (outcome.Status == CalculationStatus.Completed && outcome.Text != null)
            {
                SetState(CurveSearchState.Parse(outcome.Text));
                verified = true;
                status = "Finished. The shortlist is ordered by heuristic score; rank values are proved lower bounds.";
            }
            else if (outcome.Status == CalculationStatus.Cancelled)
                status = "Paused. Resume repeats any unfinished batch from the last checkpoint.";
            else { error = outcome.Message; status = "Search stopped. You can save the last checkpoint or resume."; }
        }
        catch (Exception e) when (e is System.IO.IOException or JsonException or ArgumentException or InvalidOperationException)
        { error = e.Message; status = "Search stopped; the last checkpoint is retained."; }
        finally
        {
            running = null;
            workbench.EndExternal();
            Notify();
        }
    }

    public void Pause()
    {
        if (running == null) return;
        status = "Pausing and stopping the worker…";
        running.Cancel();
        Notify();
    }

    public void NewSearch()
    {
        if (IsBusy) throw new InvalidOperationException("Pause the search first.");
        state = null; verified = true; dirty = false; error = "";
        Candidates.Clear(); Selected = null;
        status = "Ready. Adjust the range or start with the current settings.";
        Notify();
    }

    public void Load(CurveSearchState saved)
    {
        if (IsBusy) throw new InvalidOperationException("Pause the search first.");
        saved.Validate();
        var o = saved.Options;
        int[] values = { o.NumeratorMin, o.NumeratorMax, o.DenominatorMax, o.Workers, o.ScorePrimeBound, o.CertificatePrimeBound, o.KeepBest, o.ExtraSearchDepth };
        for (int i = 0; i < values.Length; i++) Settings[i].Text = values[i].ToString(CultureInfo.InvariantCulture);
        SetState(saved);
        verified = false; dirty = false; error = "";
        status = "Saved report loaded. Resume / recheck reconstructs the curves and verifies the retained points.";
        Notify();
    }

    public void MarkSaved() { dirty = false; Notify(); }
    public void ShowError(string message) { error = message; Notify(); }
    private void SetState(CurveSearchState next)
    {
        dirty |= state == null || state.NextSlot != next.NextSlot || !state.Results.SequenceEqual(next.Results);
        state = next;
        if (Candidates.SequenceEqual(next.Results)) return;
        var previous = Selected;
        Candidates.Clear();
        foreach (var candidate in next.Results) Candidates.Add(candidate);
        Selected = Candidates.FirstOrDefault(c => c.Numerator == previous?.Numerator && c.Denominator == previous?.Denominator) ?? Candidates.FirstOrDefault();
    }

    private void SettingChanged(object? sender, PropertyChangedEventArgs e) => Notify();
    private void WorkbenchChanged(object? sender, PropertyChangedEventArgs e) { if (e.PropertyName == nameof(WorkbenchViewModel.CanRun)) Notify(); }
    private void Notify()
    {
        foreach (var name in new[] { nameof(IsBusy), nameof(IsIdle), nameof(CanEdit), nameof(CanStart), nameof(CanOpenCurve), nameof(StartLabel), nameof(InputError), nameof(Status), nameof(Error), nameof(Percent), nameof(ProgressText), nameof(Details), nameof(IsDirty), nameof(PersistenceText) }) OnPropertyChanged(name);
    }
    public void Dispose()
    {
        disposed = true; Pause();
        workbench.PropertyChanged -= WorkbenchChanged;
        foreach (var field in Settings) field.PropertyChanged -= SettingChanged;
    }
}

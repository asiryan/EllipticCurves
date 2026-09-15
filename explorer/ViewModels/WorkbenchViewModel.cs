#nullable enable
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class WorkbenchViewModel(CalculationRunner? runner = null) : ObservableObject, IDisposable
{
    private readonly CalculationRunner runner = runner ?? new();
    private CancellationTokenSource? running;
    private bool disposed;
    private int historyChangeDepth;
    private CalculationJobViewModel? selected, active;
    private readonly ConditionalWeakTable<CalculationJobViewModel, ResultMemento> savedResults = new();
    internal event Action? HistoryChanging;
    internal event Action? HistoryChanged;
    internal event Action? HistoryReplaced;
    public ObservableCollection<CalculationJobViewModel> Jobs { get; } = new();
    public CalculationJobViewModel? Selected
    {
        get => selected;
        set
        {
            if (ReferenceEquals(selected, value)) return;
            ChangeResults(() =>
            {
                selected = value;
                OnPropertyChanged(nameof(Selected));
                OnPropertyChanged(nameof(HasSelection));
            });
        }
    }

    public CalculationJobViewModel? Active
    {
        get => active;
        private set
        {
            active = value;
            OnPropertyChanged();
        }
    }

    public bool IsBusy => running != null;
    public bool CanRun => !IsBusy && !disposed;
    public bool HasSelection => Selected != null;
    public bool HasResults => Jobs.Count > 0;
    public bool CanClearHistory => CanRun && HasResults;
    public string HistoryHeading => $"SESSION HISTORY · LAST {ExplorerSession.HistoryLimit}";
    public bool CanDelete(CalculationJobViewModel? job) => job != null && job != Active && Jobs.Contains(job);
    public string Summary => IsBusy ? "Calculation in progress" : Jobs.Count == 0 ? "Choose a calculation in Tools" : Jobs.Count + " calculations this session";
    public RelayCommand CancelCommand => new(_ => Cancel());

    public void Delete(CalculationJobViewModel? job)
    {
        if (!CanDelete(job)) return;
        ChangeResults(() =>
        {
            var index = Jobs.IndexOf(job!);
            var wasSelected = job == Selected;
            Jobs.RemoveAt(index);
            if (wasSelected) Selected = Jobs.Count == 0 ? null : Jobs[Math.Min(index, Jobs.Count - 1)];
            NotifyState();
        });
    }

    public void ClearHistory()
    {
        if (!CanClearHistory) return;
        ChangeResults(() =>
        {
            Jobs.Clear();
            Selected = null;
            NotifyState();
        });
    }

    public void RestoreHistory(IReadOnlyList<CalculationSession> history)
    {
        if (!CanRun) throw new InvalidOperationException(SessionMessages.StopCalculationBeforeOpen);
        var restored = history.Select(CalculationJobViewModel.FromSession).ToArray();
        ChangeResults(() =>
        {
            Jobs.Clear();
            foreach (var job in restored) Jobs.Add(job);
            // History is stored newest first; browsing another report is temporary.
            Selected = Jobs.FirstOrDefault();
            NotifyState();
        }, recordHistory: false);
        HistoryReplaced?.Invoke();
    }

    public async Task RunAsync(CalculationRequest request)
    {
        if (!CanRun) throw new InvalidOperationException("A calculation is already running.");
        request = request with { Arguments = new Dictionary<string, string>(request.Arguments) };
        using var cancellation = new CancellationTokenSource();
        using var clock = new CancellationTokenSource();
        var job = new CalculationJobViewModel(request, CalculationCatalog.Get(request.OperationId).Title);
        ChangeResults(() =>
        {
            running = cancellation;
            Jobs.Insert(0, job);
            // Keep the current session bounded; each result can contain up to 2 MB of text.
            if (Jobs.Count > ExplorerSession.HistoryLimit) Jobs.RemoveAt(Jobs.Count - 1);
            Selected = Active = job;
            NotifyState();
        });
        var watch = Stopwatch.StartNew();
        var timer = UpdateClockAsync(job, watch, clock.Token);
        var progress = new Progress<CalculationUpdate>(update =>
        {
            if (disposed || job.Status != CalculationStatus.Running) return;
            job.Stage = update.Message;
            job.Percent = update.Percent;
        });
        try
        {
            var outcome = await runner.RunAsync(request, update => ((IProgress<CalculationUpdate>)progress).Report(update), cancellation.Token);
            job.Status = outcome.Status;
            job.Stage = outcome.Message;
            job.Elapsed = watch.Elapsed;
            job.Percent = outcome.Status == CalculationStatus.Completed ? 100 : 0;
            job.Result = outcome.Text ?? outcome.Message;
        }
        finally
        {
            clock.Cancel();
            await timer;
            // Seal the shared result before undo becomes available. Every
            // checkpoint containing this job now sees its final report.
            if (savedResults.TryGetValue(job, out var saved)) saved.Complete(job);
            running = null;
            Active = null;
            NotifyState();
        }
    }

    private static async Task UpdateClockAsync(CalculationJobViewModel job, Stopwatch watch, CancellationToken token)
    {
        try
        {
            while (true)
            {
                await Task.Delay(500, token);
                job.Elapsed = watch.Elapsed;
            }
        }
        catch (OperationCanceledException) { }
    }

    public void Cancel()
    {
        if (Active != null) Active.Stage = "Stopping calculation…";
        running?.Cancel();
    }

    private void ChangeResults(Action change, bool recordHistory = true)
    {
        // WPF may change SelectedItem while Jobs is being updated. Keep these
        // automatic selections inside the deletion, clear, run or restore action.
        var record = recordHistory && historyChangeDepth == 0;
        historyChangeDepth++;
        try
        {
            if (record) HistoryChanging?.Invoke();
            change();
        }
        finally { historyChangeDepth--; }
        if (record) HistoryChanged?.Invoke();
    }

    private void NotifyState()
    {
        OnPropertyChanged(nameof(IsBusy));
        OnPropertyChanged(nameof(CanRun));
        OnPropertyChanged(nameof(HasResults));
        OnPropertyChanged(nameof(Summary));
        OnPropertyChanged(nameof(CanClearHistory));
    }

    public void Dispose()
    {
        disposed = true;
        Cancel();
    }

    internal sealed class ResultMemento
    {
        private CalculationSession saved;
        public ResultMemento(CalculationJobViewModel job) => saved = job.CaptureSession();
        internal void Complete(CalculationJobViewModel job) => saved = job.CaptureSession();
        internal CalculationJobViewModel Restore() => CalculationJobViewModel.FromSession(saved);
    }

    internal sealed record Memento(ResultMemento[] Results, int Selection)
    {
        public bool SameEdit(Memento other) => Selection == other.Selection && Results.SequenceEqual(other.Results);
    }

    internal Memento CaptureMemento() => new(Jobs.Select(job => savedResults.GetValue(job, value => new(value))).ToArray(),
        Selected == null ? -1 : Jobs.IndexOf(Selected));

    internal void RestoreMemento(Memento state)
    {
        if (!CanRun) throw new InvalidOperationException(SessionMessages.StopCalculationBeforeOpen);
        ChangeResults(() =>
        {
            Jobs.Clear();
            foreach (var result in state.Results)
            {
                var job = result.Restore();
                savedResults.Add(job, result);
                Jobs.Add(job);
            }
            Selected = state.Selection >= 0 && state.Selection < Jobs.Count ? Jobs[state.Selection] : null;
            NotifyState();
        }, recordHistory: false);
    }
}

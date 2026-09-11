#nullable enable
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Text;
using EllipticCurves.Explorer.Computations;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class CalculationJobViewModel(CalculationRequest request, string title) : ObservableObject
{
    public CalculationRequest Request { get; } = request;
    public string Title { get; } = title;
    public DateTime StartedAt { get; } = DateTime.Now;
    public string Equation => Request.Equation;
    private string status = "Running", stage = "Starting calculation", result = "Waiting for the calculation to finish…";
    private TimeSpan elapsed;
    private double? percent;
    public string Status { get => status; internal set { status = value; OnPropertyChanged(); OnPropertyChanged(nameof(HistoryLabel)); } }
    public string Stage { get => stage; internal set { stage = value; OnPropertyChanged(); } }
    public string Result { get => result; internal set { result = value; OnPropertyChanged(); OnPropertyChanged(nameof(Report)); } }
    public string HistoryLabel => StartedAt.ToString("HH:mm:ss") + " · " + Title + " · " + Status;
    public TimeSpan Elapsed { get => elapsed; internal set { elapsed = value; OnPropertyChanged(nameof(Timing)); } }
    public string Timing => elapsed.ToString(@"hh\:mm\:ss") + (Request.TimeoutSeconds > 0 ? " · limit " + Request.TimeoutSeconds + " s" : " · no time limit");
    public double? Percent { get => percent; internal set { percent = value; OnPropertyChanged(nameof(IsIndeterminate)); OnPropertyChanged(nameof(ProgressValue)); } }
    public bool IsIndeterminate => !Percent.HasValue;
    public double ProgressValue => Percent ?? 0;
    public string Report
    {
        get
        {
            var operation = CalculationCatalog.Get(Request.OperationId);
            var text = new StringBuilder().AppendLine(Title).AppendLine(StartedAt.ToString("yyyy-MM-dd HH:mm:ss zzz"))
                .AppendLine("Status: " + Status).AppendLine("Elapsed: " + Timing)
                .AppendLine(operation.UsesPlot ? "Captured plot: " + Equation : "Independent of the plotted curve.")
                .AppendLine(operation.Description).AppendLine();
            foreach (var parameter in operation.Parameters)
                text.Append(parameter.Label).Append(": ").AppendLine(Request.Arguments.GetValueOrDefault(parameter.Key, parameter.Default));
            return text.AppendLine().AppendLine(Result).ToString();
        }
    }
}

public sealed class WorkbenchViewModel(CalculationRunner? runner = null) : ObservableObject, IDisposable
{
    private readonly CalculationRunner runner = runner ?? new();
    private CancellationTokenSource? running;
    private bool disposed;
    private CalculationJobViewModel? selected, active;
    public ObservableCollection<CalculationJobViewModel> Jobs { get; } = new();
    public CalculationJobViewModel? Selected { get => selected; set { selected = value; OnPropertyChanged(); OnPropertyChanged(nameof(HasSelection)); } }
    public CalculationJobViewModel? Active { get => active; private set { active = value; OnPropertyChanged(); } }
    public bool IsBusy => running != null;
    public bool CanRun => !IsBusy && !disposed;
    public bool HasSelection => Selected != null;
    public bool HasResults => Jobs.Count > 0;
    public bool CanClearHistory => CanRun && HasResults;
    public bool CanDelete(CalculationJobViewModel? job) => job != null && job != Active && Jobs.Contains(job);
    public string Summary => IsBusy ? "Calculation in progress" : Jobs.Count == 0 ? "Choose a calculation in Explorer" : Jobs.Count + " calculations this session";
    public RelayCommand CancelCommand => new(_ => Cancel());

    public void Delete(CalculationJobViewModel? job)
    {
        if (!CanDelete(job)) return;
        var index = Jobs.IndexOf(job!);
        var wasSelected = job == Selected;
        Jobs.RemoveAt(index);
        if (wasSelected) Selected = Jobs.Count == 0 ? null : Jobs[Math.Min(index, Jobs.Count - 1)];
        NotifyState();
    }

    public void ClearHistory()
    {
        if (!CanClearHistory) return;
        Jobs.Clear();
        Selected = null;
        NotifyState();
    }

    public async Task RunAsync(CalculationRequest request)
    {
        if (!CanRun) throw new InvalidOperationException("A calculation is already running.");
        request = request with { Arguments = new Dictionary<string, string>(request.Arguments) };
        using var cancellation = new CancellationTokenSource();
        using var clock = new CancellationTokenSource();
        running = cancellation;
        var job = new CalculationJobViewModel(request, CalculationCatalog.Get(request.OperationId).Title);
        Jobs.Insert(0, job);
        // Keep the current session bounded; each result can contain up to 2 MB of text.
        if (Jobs.Count > 50) Jobs.RemoveAt(Jobs.Count - 1);
        Selected = Active = job;
        NotifyState();
        var watch = Stopwatch.StartNew();
        var timer = UpdateClockAsync(job, watch, clock.Token);
        var progress = new Progress<CalculationUpdate>(update =>
        {
            if (disposed || job.Status != "Running") return;
            job.Stage = update.Message;
            job.Percent = update.Percent;
        });
        try
        {
            var outcome = await runner.RunAsync(request, update => ((IProgress<CalculationUpdate>)progress).Report(update), cancellation.Token);
            job.Status = outcome.Status;
            job.Stage = outcome.Message;
            job.Elapsed = watch.Elapsed;
            job.Percent = outcome.Status == "Completed" ? 100 : 0;
            job.Result = outcome.Text ?? outcome.Message;
        }
        finally
        {
            clock.Cancel();
            await timer;
            running = null;
            Active = null;
            NotifyState();
        }
    }

    private static async Task UpdateClockAsync(CalculationJobViewModel job, Stopwatch watch, CancellationToken token)
    {
        try { while (true) { await Task.Delay(500, token); job.Elapsed = watch.Elapsed; } }
        catch (OperationCanceledException) { }
    }
    public void Cancel() { if (Active != null) Active.Stage = "Stopping calculation…"; running?.Cancel(); }
    private void NotifyState()
    {
        OnPropertyChanged(nameof(IsBusy)); OnPropertyChanged(nameof(CanRun)); OnPropertyChanged(nameof(HasResults)); OnPropertyChanged(nameof(Summary));
        OnPropertyChanged(nameof(CanClearHistory));
    }
    public void Dispose() { disposed = true; Cancel(); }
}

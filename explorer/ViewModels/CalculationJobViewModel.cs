#nullable enable
using System.Text;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class CalculationJobViewModel(CalculationRequest request, string title, DateTime? startedAt = null) : ObservableObject
{
    public CalculationRequest Request { get; } = request;
    public string Title { get; } = title;
    public DateTime StartedAt { get; } = startedAt ?? DateTime.Now;
    public string Equation => Request.Equation;

    private string status = CalculationStatus.Running, stage = "Starting calculation", result = "Waiting for the calculation to finish…";
    private TimeSpan elapsed;
    private double? percent;

    public string Status
    {
        get => status;
        internal set
        {
            status = value;
            OnPropertyChanged();
            OnPropertyChanged(nameof(HistoryLabel));
        }
    }

    public string Stage
    {
        get => stage;
        internal set
        {
            stage = value;
            OnPropertyChanged();
        }
    }

    public string Result
    {
        get => result;
        internal set
        {
            result = value;
            OnPropertyChanged();
            OnPropertyChanged(nameof(Report));
        }
    }

    public string HistoryLabel => StartedAt.ToString("HH:mm:ss") + " · " + Title + " · " + Status;

    public CalculationSession CaptureSession() => new(
        Request with { Arguments = new Dictionary<string, string>(Request.Arguments) },
        StartedAt, Status, Stage, Elapsed, Percent, Result);

    public static CalculationJobViewModel FromSession(CalculationSession saved) =>
        new(saved.Request with { Arguments = new Dictionary<string, string>(saved.Request.Arguments) },
            CalculationCatalog.Get(saved.Request.OperationId).Title, saved.StartedAt)
        {
            Status = saved.Status == CalculationStatus.Running ? CalculationStatus.Interrupted : saved.Status,
            Stage = saved.Status == CalculationStatus.Running ? "Saved during a calculation. Use Repeat to run it again." : saved.Stage,
            Elapsed = saved.Elapsed,
            Percent = saved.Status == CalculationStatus.Running ? null : saved.Percent,
            Result = saved.Status == CalculationStatus.Running ? "This calculation was still running when the session was saved. Use Repeat to run it again." : saved.Result
        };

    public TimeSpan Elapsed
    {
        get => elapsed;
        internal set
        {
            elapsed = value;
            OnPropertyChanged(nameof(Timing));
        }
    }

    public string Timing => elapsed.ToString(@"hh\:mm\:ss")
        + (Request.TimeoutSeconds > 0 ? " · limit " + Request.TimeoutSeconds + " s" : " · no time limit");

    public double? Percent
    {
        get => percent;
        internal set
        {
            percent = value;
            OnPropertyChanged(nameof(IsIndeterminate));
            OnPropertyChanged(nameof(ProgressValue));
        }
    }

    public bool IsIndeterminate => !Percent.HasValue;
    public double ProgressValue => Percent ?? 0;

    public string Report
    {
        get
        {
            var operation = CalculationCatalog.Get(Request.OperationId);
            var text = new StringBuilder()
                .AppendLine(Title)
                .AppendLine(StartedAt.ToString("yyyy-MM-dd HH:mm:ss zzz"))
                .AppendLine("Status: " + Status)
                .AppendLine("Elapsed: " + Timing)
                .AppendLine(operation.UsesPlot ? "Captured plot: " + Equation : "Independent of the plotted curve.")
                .AppendLine(operation.Description)
                .AppendLine();
            foreach (var parameter in operation.Parameters)
                text.Append(CalculationFormatter.FormatInput(parameter,
                    Request.Arguments.GetValueOrDefault(parameter.Key, parameter.Default)));
            return text.AppendLine().AppendLine(Result).ToString();
        }
    }
}

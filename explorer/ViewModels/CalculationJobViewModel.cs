#nullable enable
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
                text.Append(parameter.Label).Append(": ")
                    .AppendLine(Request.Arguments.GetValueOrDefault(parameter.Key, parameter.Default));
            return text.AppendLine().AppendLine(Result).ToString();
        }
    }
}

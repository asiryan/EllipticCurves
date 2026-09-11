#nullable enable
using System.ComponentModel;
using System.Globalization;
using EllipticCurves.Visualizer.Computations;

namespace EllipticCurves.Visualizer.ViewModels;

public sealed class CalculationFieldViewModel(CalculationParameter parameter) : ObservableObject
{
    private string text = parameter.Default;
    private bool relevant = true;
    public CalculationParameter Parameter { get; } = parameter;
    public string Label => Parameter.Label;
    public string Help => Parameter.Help;
    public bool IsBoolean => Parameter.Kind == ParameterKind.Boolean;
    public bool IsMultiline => Parameter.Kind == ParameterKind.Multiline;
    public bool IsText => !IsBoolean;
    public bool IsRelevant { get => relevant; set { relevant = value; OnPropertyChanged(); OnPropertyChanged(nameof(Error)); } }
    public string Text
    {
        get => text;
        set { if (text == value) return; text = value; OnPropertyChanged(); OnPropertyChanged(nameof(Checked)); OnPropertyChanged(nameof(Error)); }
    }
    public bool Checked { get => bool.TryParse(Text, out var value) && value; set => Text = value.ToString(); }
    public string Error => IsRelevant ? CalculationInput.Validate(Parameter, Text) : "";
}

public sealed class CalculationFormViewModel : ObservableObject, IDisposable
{
    private readonly WorkbenchViewModel workbench;
    public CalculationOperation Operation { get; }
    public string Equation { get; }
    public IReadOnlyList<CalculationFieldViewModel> Fields { get; }
    public IEnumerable<CalculationFieldViewModel> BasicFields => Fields.Where(f => !f.Parameter.Advanced);
    public IEnumerable<CalculationFieldViewModel> AdvancedFields => Fields.Where(f => f.Parameter.Advanced);
    public bool HasAdvancedFields => AdvancedFields.Any();
    private string timeout = "120", maxItems = "1000";
    public string Timeout { get => timeout; set { timeout = value; Changed(); } }
    public string MaxItems { get => maxItems; set { maxItems = value; Changed(); } }
    public string LimitError => !int.TryParse(Timeout, out var seconds) || seconds is < 0 or > 86400
        ? "Time limit: enter 0–86,400 seconds. 0 means unlimited."
        : !int.TryParse(MaxItems, out var count) || count is < 1 or > 100_000 ? "Result limit: enter 1–100,000 items." : "";
    public bool CanRun => workbench.CanRun && LimitError.Length == 0 && Fields.All(f => f.Error.Length == 0);
    public string RunHint => !workbench.CanRun ? "A calculation is running. Stop it or wait for it to finish."
        : CanRun ? "The result will appear in the Results panel. You can keep editing the plot."
        : "Correct the highlighted input before running.";

    public CalculationFormViewModel(CalculationOperation operation, string equation, WorkbenchViewModel workbench, CalculationRequest? previous = null)
    {
        Operation = operation; Equation = equation; this.workbench = workbench;
        Fields = operation.Parameters.Select(p => new CalculationFieldViewModel(p)).ToArray();
        if (previous != null)
        {
            timeout = previous.TimeoutSeconds.ToString(CultureInfo.InvariantCulture);
            maxItems = previous.MaxItems.ToString(CultureInfo.InvariantCulture);
            foreach (var field in Fields) if (previous.Arguments.TryGetValue(field.Parameter.Key, out var value)) field.Text = value;
        }
        foreach (var field in Fields) field.PropertyChanged += FieldChanged;
        workbench.PropertyChanged += WorkbenchChanged;
        UpdateRelevance();
    }

    private void FieldChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName != nameof(CalculationFieldViewModel.Text)) return;
        UpdateRelevance(); Changed();
    }
    private void WorkbenchChanged(object? sender, PropertyChangedEventArgs e)
    { if (e.PropertyName == nameof(WorkbenchViewModel.CanRun)) Changed(); }
    private void UpdateRelevance()
    {
        var values = Fields.ToDictionary(f => f.Parameter.Key, f => f.Text);
        foreach (var field in Fields) field.IsRelevant = !CalculationEngine.IsIgnoredCoordinate(field.Parameter.Key, values);
    }
    private void Changed()
    {
        OnPropertyChanged(nameof(Timeout)); OnPropertyChanged(nameof(MaxItems));
        OnPropertyChanged(nameof(LimitError)); OnPropertyChanged(nameof(CanRun)); OnPropertyChanged(nameof(RunHint));
    }
    public CalculationRequest CreateRequest() => new(Operation.Id, Equation, Fields.ToDictionary(f => f.Parameter.Key, f => f.Text),
        int.Parse(Timeout, CultureInfo.InvariantCulture), int.Parse(MaxItems, CultureInfo.InvariantCulture));
    public void Dispose()
    {
        foreach (var field in Fields) field.PropertyChanged -= FieldChanged;
        workbench.PropertyChanged -= WorkbenchChanged;
    }
}

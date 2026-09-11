#nullable enable
using EllipticCurves.Explorer.Computations;

namespace EllipticCurves.Explorer.ViewModels;

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

    public bool IsRelevant
    {
        get => relevant;
        set
        {
            relevant = value;
            OnPropertyChanged();
            OnPropertyChanged(nameof(Error));
        }
    }

    public string Text
    {
        get => text;
        set
        {
            if (text == value) return;
            text = value;
            OnPropertyChanged();
            OnPropertyChanged(nameof(Checked));
            OnPropertyChanged(nameof(Error));
        }
    }

    public bool Checked
    {
        get => bool.TryParse(Text, out var value) && value;
        set => Text = value.ToString();
    }

    public string Error => IsRelevant ? CalculationInput.Validate(Parameter, Text) : "";
}

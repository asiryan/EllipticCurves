#nullable enable
using System.ComponentModel;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class EquationViewModel(Action changed) : ObservableObject, IDataErrorInfo
{
    private string text = "y^2 = x^3 - x", parseError = "";
    private bool editing;
    public EllipticCurveQ Curve { get; private set; } = new(0, 0, 0, -1, 0);
    public bool IsValid => parseError.Length == 0;
    public string Error => editing ? "" : parseError;
    public string this[string columnName] => columnName == nameof(Text) ? Error : "";

    public string Text
    {
        get => text;
        set
        {
            text = value;
            editing = true;
            if (CurveEquationText.TryParse(value, out var curve, out parseError)) Curve = curve!;
            Notify();
        }
    }

    public void CommitEdit()
    {
        editing = false;
        Notify();
    }

    public void SetCurve(EllipticCurveQ curve)
    {
        Curve = curve;
        text = CurveEquationText.Format(curve);
        parseError = "";
        editing = false;
        Notify();
    }

    private void Notify()
    {
        OnPropertyChanged(nameof(Text));
        OnPropertyChanged(nameof(IsValid));
        OnPropertyChanged(nameof(Error));
        changed();
    }
}

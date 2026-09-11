#nullable enable
using System.ComponentModel;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class EquationViewModel(Action changed) : ObservableObject, IDataErrorInfo
{
    private string text = CurvePreset.ClassicEquation, parseError = "";
    private bool editing;
    public EllipticCurveQ Curve { get; private set; } = CurvePreset.Classic.CreateCurve();
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

    internal sealed record Memento(string Text, EllipticCurveQ Curve, string Error, bool Editing);
    internal Memento CaptureMemento() => new(text, Curve, parseError, editing);
    internal void RestoreMemento(Memento state)
    {
        (text, Curve, parseError, editing) = (state.Text, state.Curve, state.Error, state.Editing);
        Notify();
    }
}

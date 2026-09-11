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
        text = CurveEquationText.Format(curve);
        // Slider arithmetic can cross the coefficient limit too. Never mark an
        // equation as valid if calculations or session loading would reject it.
        if (CurveEquationText.TryParse(text, out var parsed, out parseError)) Curve = parsed!;
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

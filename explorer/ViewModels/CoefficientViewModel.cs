#nullable enable
using System.ComponentModel;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class CoefficientViewModel : ObservableObject, IDataErrorInfo
{
    private readonly Action changed;
    private readonly Func<BigRational?> getStep;
    private readonly bool requirePositive;
    private BigRational exactValue, sliderAnchor;
    private int sliderOffset;
    private string text = "0", parseError = "";
    private bool editing;
    public string Name { get; }
    public string Term { get; }
    public BigRational ExactValue => exactValue;
    public bool IsValid => parseError.Length == 0;
    public string Error => editing ? "" : parseError;
    public string this[string columnName] => columnName == nameof(Text) ? Error : "";
    public RelayCommand ResetCommand { get; }

    public CoefficientViewModel(string name, string term, Action changed, Func<BigRational?>? getStep = null, bool requirePositive = false)
    {
        Name = name;
        Term = term;
        this.changed = changed;
        this.getStep = getStep ?? (() => new BigRational(1, 10));
        this.requirePositive = requirePositive;
        ResetCommand = new RelayCommand(_ => SetExact(BigRational.Zero));
    }

    public string Text
    {
        get => text;
        set
        {
            text = value;
            editing = true;
            if (RationalText.TryParse(value, out var parsed) && (!requirePositive || parsed.Sign > 0))
            {
                exactValue = parsed;
                parseError = "";
                RecenterSlider();
            }
            else parseError = requirePositive ? "Enter a positive step, such as 0.01 or 1/7."
                : "Enter a number or fraction, such as 8.325, -2/7 or 1e-5 (up to 4096 characters; exponent ±4096).";
            Notify();
        }
    }

    public void CommitEdit()
    {
        editing = false;
        Notify();
    }

    public void SetExact(BigRational value)
    {
        exactValue = value;
        text = RationalText.Format(value);
        parseError = "";
        editing = false;
        RecenterSlider();
        Notify();
    }

    // The slider changes an exact anchor by an integer number of exact steps.
    public double SliderOffset
    {
        get => sliderOffset;
        set
        {
            if (!double.IsFinite(value) || !IsValid || getStep() is not BigRational step) return;
            var next = (int)Math.Round(Math.Clamp(value, -50, 50), MidpointRounding.AwayFromZero);
            if (next == sliderOffset) return;
            sliderOffset = next;
            exactValue = sliderAnchor + next * step;
            text = RationalText.Format(exactValue);
            editing = false;
            Notify();
        }
    }

    public string SliderMinimum => RationalText.Format(sliderAnchor - 50 * (getStep() ?? BigRational.Zero));
    public string SliderMaximum => RationalText.Format(sliderAnchor + 50 * (getStep() ?? BigRational.Zero));

    public void RecenterSlider()
    {
        sliderAnchor = exactValue;
        sliderOffset = 0;
        OnPropertyChanged(nameof(SliderOffset));
        OnPropertyChanged(nameof(SliderMinimum));
        OnPropertyChanged(nameof(SliderMaximum));
    }

    private void Notify()
    {
        OnPropertyChanged(nameof(Text));
        OnPropertyChanged(nameof(ExactValue));
        OnPropertyChanged(nameof(SliderOffset));
        OnPropertyChanged(nameof(IsValid));
        OnPropertyChanged(nameof(Error));
        changed();
    }
}

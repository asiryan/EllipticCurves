using System.ComponentModel;
using System.Globalization;

namespace EllipticCurves.Visualizer.ViewModels;

public sealed class CoefficientViewModel(string name, string term, Action changed) : ObservableObject, IDataErrorInfo
{
    private int tenths;
    private string text = "0";
    public string Name { get; } = name;
    public string Term { get; } = term;
    public BigRational ExactValue => new(tenths, 10);
    public bool IsValid => Error.Length == 0;
    public string Error { get; private set; } = "";
    public string this[string columnName] => columnName == nameof(Text) ? Error : "";

    public double Value
    {
        get => tenths / 10.0;
        set
        {
            if (!double.IsFinite(value)) return;
            var next = (int)Math.Round(Math.Clamp(value, -100, 100) * 10, MidpointRounding.AwayFromZero);
            if (next == tenths && IsValid) return;
            tenths = next;
            text = (tenths / 10m).ToString("0.#", CultureInfo.InvariantCulture);
            Error = "";
            Notify();
        }
    }

    public string Text
    {
        get => text;
        set
        {
            text = value;
            if (decimal.TryParse(value.Replace(',', '.'), NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint,
                    CultureInfo.InvariantCulture, out var parsed) && parsed is >= -100 and <= 100 && parsed * 10 == decimal.Truncate(parsed * 10))
            {
                tenths = (int)(parsed * 10);
                Error = "";
            }
            else Error = "Enter −100 to 100, with at most one decimal place.";
            Notify();
        }
    }

    private void Notify()
    {
        OnPropertyChanged(nameof(Value));
        OnPropertyChanged(nameof(Text));
        OnPropertyChanged(nameof(IsValid));
        OnPropertyChanged(nameof(Error));
        changed();
    }
}

#nullable enable
using System.Globalization;

namespace EllipticCurves.Explorer.Models;

public sealed record LmfdbConductorRange(int Minimum, int Maximum)
{
    public const int Limit = 500_000;

    public static bool TryParse(string text, out LmfdbConductorRange? range)
    {
        range = null;
        var parts = text.Trim().Replace(" ", "").Replace("\u00a0", "").Replace('–', '-').Split('-');
        if (parts.Length is < 1 or > 2
            || !int.TryParse(parts[0], NumberStyles.None, CultureInfo.InvariantCulture, out var minimum)) return false;
        var maximum = minimum;
        if (parts.Length == 2 && !int.TryParse(parts[1], NumberStyles.None, CultureInfo.InvariantCulture, out maximum)) return false;
        if (minimum < 1 || maximum < minimum || maximum > Limit) return false;
        range = new(minimum, maximum);
        return true;
    }

    public override string ToString() => Minimum == Maximum ? Minimum.ToString(CultureInfo.InvariantCulture)
        : FormattableString.Invariant($"{Minimum}–{Maximum}");
}

// The import keeps only the equation. Labels identify choices in the online picker.
public sealed record LmfdbCurveFormula(string Label, int Conductor, string Equation);

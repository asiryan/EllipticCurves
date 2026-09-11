using System.Globalization;
using System.Numerics;
using System.Text.RegularExpressions;

namespace EllipticCurves.Visualizer.Models;

/// <summary>Exact text conversion for the editor, without a floating-point intermediate.</summary>
public static class RationalText
{
    // Bound pasted input and exponent expansion, not the slider's coefficient range.
    private const int MaxDigits = 4096;
    private static readonly Regex Number = new(@"^([+-]?)([0-9]*)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?$",
        RegexOptions.CultureInvariant | RegexOptions.NonBacktracking);

    public static bool TryParse(string text, out BigRational value)
    {
        value = BigRational.Zero;
        if (string.IsNullOrWhiteSpace(text) || text.Length > MaxDigits) return false;
        var parts = text.Trim().Replace('−', '-').Replace(',', '.').Split('/');
        if (parts.Length > 2 || !TryNumber(parts[0].Trim(), out var numerator)) return false;
        if (parts.Length == 1) { value = numerator; return true; }
        if (!TryNumber(parts[1].Trim(), out var denominator) || denominator.IsZero) return false;
        value = numerator / denominator;
        return true;
    }

    private static bool TryNumber(string text, out BigRational value)
    {
        value = BigRational.Zero;
        var match = Number.Match(text);
        if (!match.Success) return false;
        var digits = match.Groups[2].Value + match.Groups[3].Value;
        if (digits.Length == 0) return false;
        var exponent = 0;
        if (match.Groups[4].Success && (!int.TryParse(match.Groups[4].Value, NumberStyles.AllowLeadingSign,
                CultureInfo.InvariantCulture, out exponent) || Math.Abs((long)exponent) > MaxDigits)) return false;
        var numerator = BigInteger.Parse(digits, CultureInfo.InvariantCulture);
        if (match.Groups[1].Value == "-") numerator = -numerator;
        var scale = match.Groups[3].Length - exponent;
        value = scale >= 0 ? new BigRational(numerator, BigInteger.Pow(10, scale))
            : new BigRational(numerator * BigInteger.Pow(10, -scale));
        return true;
    }

    public static string Format(BigRational value)
    {
        var denominator = value.Den;
        var twos = 0;
        var fives = 0;
        while (denominator % 2 == 0) { denominator /= 2; twos++; }
        while (denominator % 5 == 0) { denominator /= 5; fives++; }
        if (!denominator.IsOne) return value.ToString();
        var places = Math.Max(twos, fives);
        if (places == 0) return value.ToString();
        var scaled = BigInteger.Abs(value.Num) * BigInteger.Pow(2, places - twos) * BigInteger.Pow(5, places - fives);
        var digits = scaled.ToString(CultureInfo.InvariantCulture).PadLeft(places + 1, '0');
        return (value.Sign < 0 ? "-" : "") + digits.Insert(digits.Length - places, ".");
    }
}

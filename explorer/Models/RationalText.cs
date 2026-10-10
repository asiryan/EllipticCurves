using System.Globalization;
using System.Numerics;

namespace EllipticCurves.Explorer.Models;

/// <summary>Exact text conversion for the editor, without a floating-point intermediate.</summary>
public static class RationalText
{
    // Leave room for the exact numerator and denominator of an editor coefficient.
    public const int MaxTextLength = 20_000;
    internal const int MaxValueBits = 32768;
    private const int MaxExponent = 4096;
    private static readonly RationalParsingOptions ParsingLimits = new()
    {
        MaxTextLength = MaxTextLength, MaxExponent = MaxExponent,
        MaxValueBits = MaxValueBits, MaxIntermediateBits = 131072
    };

    public static bool TryParse(string text, out BigRational value)
        => BigRational.TryParse(text, ParsingLimits, out value);

    public static string Format(BigRational value)
    {
        var denominator = value.Den;
        var twos = 0;
        var fives = 0;
        while (denominator % 2 == 0) { denominator /= 2; twos++; }
        while (denominator % 5 == 0) { denominator /= 5; fives++; }
        if (!denominator.IsOne) return value.ToString();
        var places = Math.Max(twos, fives);
        // A power-of-two denominator can need far more decimal places than its
        // rational representation. Keep generated text inside the parser budget.
        if (places > MaxExponent) return value.ToString();
        if (places == 0) return value.ToString();
        var scaled = BigInteger.Abs(value.Num) * BigInteger.Pow(2, places - twos) * BigInteger.Pow(5, places - fives);
        var digits = scaled.ToString(CultureInfo.InvariantCulture).PadLeft(places + 1, '0');
        return (value.Sign < 0 ? "-" : "") + digits.Insert(digits.Length - places, ".");
    }

    internal static bool IsWithinLimit(BigRational value) =>
        BigInteger.Abs(value.Num).GetBitLength() <= MaxValueBits && value.Den.GetBitLength() <= MaxValueBits;
}

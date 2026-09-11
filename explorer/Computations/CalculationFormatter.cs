#nullable enable
using System.Collections;
using System.Globalization;
using System.Numerics;
using System.Reflection;
using System.Text;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Computations;

public static class CalculationFormatter
{
    public static string Format(object? value, int maxItems, Action<CalculationUpdate>? report = null, CancellationToken token = default)
    {
        var writer = new ResultWriter(maxItems, report, token);
        writer.Write("Result", value, 0);
        return writer.Finish();
    }

    private sealed class ResultWriter(int maxItems, Action<CalculationUpdate>? report, CancellationToken token)
    {
        private const int MaxCharacters = 2_000_000;
        private readonly StringBuilder text = new();
        private readonly HashSet<object> ancestors = new(ReferenceEqualityComparer.Instance);
        private bool shortened;
        private void Line(int depth, string line)
        {
            if (text.Length + line.Length > MaxCharacters) { shortened = true; return; }
            text.Append(' ', depth * 2).AppendLine(line);
        }
        public string Finish() => text + (shortened ? "\nOUTPUT TRUNCATED by the item or text limit. This is not a complete list. Increase the result limit and run again.\n" : "");
        public void Write(string label, object? value, int depth)
        {
            token.ThrowIfCancellationRequested();
            if (text.Length >= MaxCharacters || depth > 14) { shortened = true; return; }
            if (value == null) { Line(depth, label + ": unavailable / not applicable"); return; }
            if (value is EllipticCurveQ curve) { Line(depth, label + ": " + CurveEquationText.Format(curve)); return; }
            if (value is EllipticCurvePoint or EllipticCurvePointFp or EllipticCurvePointFq or FiniteFieldElement)
            { Line(depth, label + ": " + value); return; }
            if (value is string or BigInteger or BigRational || value.GetType().IsPrimitive || value.GetType().IsEnum)
            {
                var scalar = value is double real ? real.ToString("G17", CultureInfo.InvariantCulture) : Convert.ToString(value, CultureInfo.InvariantCulture);
                Line(depth, label + ": " + scalar); return;
            }
            if (!value.GetType().IsValueType && !ancestors.Add(value)) { Line(depth, label + ": (already shown)"); return; }
            try
            {
                Line(depth, label + (value is RealEnclosure ? " · certified enclosure" : "") + ":");
                if (value is LmfdbRealValue stored)
                {
                    Line(depth + 1, "Stored decimal approximation, not a certified error interval.");
                    Write("Exact rational representation of stored decimal", stored.AsRational(), depth + 1);
                }
                if (value is IDictionary dictionary)
                {
                    foreach (DictionaryEntry item in dictionary) Write(item.Key.ToString()!, item.Value, depth + 1);
                }
                else if (value is Array { Rank: 2 } matrix)
                {
                    Line(depth + 1, "Dimensions: " + matrix.GetLength(0) + " × " + matrix.GetLength(1) + " (indices start at 0)");
                    var items = 0;
                    for (var row = 0; row < matrix.GetLength(0); row++)
                        for (var column = 0; column < matrix.GetLength(1); column++)
                        {
                            if (items++ >= maxItems) { shortened = true; return; }
                            Write("[" + row + ", " + column + "]", matrix.GetValue(row, column), depth + 1);
                        }
                }
                else if (value is IEnumerable sequence)
                {
                    var index = 0;
                    var count = value is ICollection collection ? collection.Count : (int?)null;
                    foreach (var item in sequence)
                    {
                        if (index >= maxItems || text.Length >= MaxCharacters) { shortened = true; break; }
                        Write("[" + index + "]", item, depth + 1);
                        index++;
                        if (index % 100 == 0) report?.Invoke(new("progress",
                            count.HasValue ? "Formatting items · " + index + " / " + count : "Collecting items · " + index,
                            count is > 0 ? 100.0 * index / count.Value : null));
                    }
                    Line(depth + 1, "Items displayed: " + index + (count.HasValue ? " / " + count : ""));
                }
                else
                {
                    // Result DTOs and finite-field objects only. Never inspect the expensive
                    // computed properties of a rational curve (handled above).
                    foreach (var property in value.GetType().GetProperties(BindingFlags.Public | BindingFlags.Instance).Where(p => p.CanRead && p.GetIndexParameters().Length == 0))
                        Write(CalculationOperation.Humanize(property.Name), property.GetValue(value), depth + 1);
                }
            }
            finally { ancestors.Remove(value); }
        }
    }
}

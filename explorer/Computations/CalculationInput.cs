#nullable enable
using System.Net.Http;
using System.Globalization;
using System.Numerics;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Computations;

public static class CalculationInput
{
    public const string ElementHelp = "A scalar or coefficients in ascending powers of t separated by ;. Example: 0; 1 means t.";
    private static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;
    public static int DefaultRankWorkers => Math.Min(4, Math.Max(1, Environment.ProcessorCount - 1));
    public const string RankWorkersHelp = "Maximum parallel workers for general rank descent. 1 runs sequentially. All workers share the work limits; the 2-isogeny method stays sequential.";
    public const string FactorizationWorkersKey = "options.MaxDegreeOfParallelism";
    // A default, not a ceiling: more workers can lose time to memory traffic.
    public static int DefaultFactorizationWorkers => Math.Min(4, Environment.ProcessorCount);
    public static string FactorizationWorkersHelp => $"Maximum parallel workers for integer factorization. 1 runs sequentially. Up to {Environment.ProcessorCount} CPU workers are available; smaller inputs may use fewer workers.";
    public static bool IsOptions(Type type) => type == typeof(RankComputationOptions) || type == typeof(AnalyticRankOptions)
        || type == typeof(RealComputationOptions) || type == typeof(PointDivisionOptions) || type == typeof(SaturationOptions);
    public static bool IsPoint(Type type) => type == typeof(EllipticCurvePoint) || type == typeof(EllipticCurvePointFp) || type == typeof(EllipticCurvePointFq);
    internal static bool IsPoints(Type type) => type == typeof(IReadOnlyList<EllipticCurvePoint>) || type == typeof(IEnumerable<EllipticCurvePoint>);
    private static bool IsIntegers(Type type) => type == typeof(IReadOnlyList<int>) || type == typeof(IReadOnlyList<BigInteger>) || type == typeof(BigInteger[]);

    public static IEnumerable<CalculationParameter> Describe(Type type, string key, object? value = null, bool advanced = false)
    {
        if (IsOptions(type))
        {
            var defaults = value ?? (type == typeof(RankComputationOptions)
                ? new RankComputationOptions { MaxDegreeOfParallelism = DefaultRankWorkers } : Activator.CreateInstance(type)!);
            foreach (var property in type.GetProperties().Where(p => p.CanWrite))
                foreach (var field in Describe(property.PropertyType, key + "." + property.Name, property.GetValue(defaults), true))
                    yield return field;
            yield break;
        }
        if (key.EndsWith("." + nameof(RankComputationOptions.MaxDegreeOfParallelism), StringComparison.Ordinal))
            advanced = false;
        var label = CalculationOperation.Humanize(key.Replace(".", " · "));
        if (IsPoint(type))
        {
            yield return new(key + ".infinity", label + " · point at infinity", "False", "O is the neutral element.", typeof(bool), advanced, ParameterKind.Boolean);
            var coordinate = type == typeof(EllipticCurvePoint) ? typeof(BigRational) : type == typeof(EllipticCurvePointFp) ? typeof(BigInteger) : typeof(FiniteFieldElement);
            foreach (var axis in new[] { "x", "y" })
                yield return new(key + "." + axis, label + " · " + axis, "0", coordinate == typeof(FiniteFieldElement) ? ElementHelp : "Exact coordinate; ignored when point at infinity is checked.", coordinate, advanced);
            yield break;
        }
        var initial = value != null ? Convert.ToString(value, Invariant)! : Default(type, key);
        var help = key.EndsWith("." + nameof(RankComputationOptions.MaxDegreeOfParallelism), StringComparison.Ordinal) ? RankWorkersHelp
            : type == typeof(BigRational) ? "Exact decimal, fraction or scientific notation. Decimal commas are accepted."
            : type == typeof(EllipticCurveQ) ? "Full Weierstrass equation. Use ^ for powers."
            : type == typeof(FiniteFieldElement) ? ElementHelp
            : IsPoints(type) ? "One point per line: x, y or (x, y). Fractions are accepted. Use x; y for decimal commas (1,5; 2,5). Use O for infinity. An empty list is allowed."
            : IsIntegers(type) ? "Integers separated by ; or spaces."
            : type == typeof(bool) ? "Enable this option."
            : type == typeof(double) ? "Finite decimal or scientific notation."
            : "Integer value. Work and precision limits are enforced by the library.";
        yield return new(key, label, initial, help, type, advanced,
            type == typeof(bool) ? ParameterKind.Boolean : IsPoints(type) ? ParameterKind.Multiline : ParameterKind.Text);
    }

    private static string Default(Type type, string key) => type == typeof(EllipticCurveQ) ? CurvePreset.ClassicEquation
        : IsPoints(type) ? "0; 0" : IsIntegers(type) ? "2; 3; 5"
        : type == typeof(bool) ? "False" : type == typeof(FiniteFieldElement) ? (key is "x" or "y" ? "0" : "1")
        : key is "u" or "d" or "den" or "a" or "b" ? "1" : key is "n" or "k" or "exponent" ? "2"
        : key is "prime" ? "5" : key is "numMax" or "xmax" ? "12" : key == "denMax" ? "4"
        : key is "count" or "index" ? "10" : "0";

    public static string Validate(CalculationParameter parameter, string text)
    {
        try
        {
            var value = ParseScalar(parameter.ValueType, text);
            if (parameter.Key.EndsWith("." + nameof(RankComputationOptions.MaxDegreeOfParallelism), StringComparison.Ordinal)
                && (int)value < 1) return parameter.Label + ": enter at least 1 worker.";
            return "";
        }
        catch (Exception error) when (error is FormatException or OverflowException or ArgumentException)
        { return parameter.Label + ": " + error.Message; }
    }

    public static object ParseScalar(Type type, string text)
    {
        if (text.Length > 100_000) throw new FormatException("Input is limited to 100,000 characters.");
        text = text.Trim();
        if (type == typeof(string)) return text.Length > 0 ? text : throw new FormatException("Enter a value.");
        if (type == typeof(bool)) return bool.Parse(text);
        if (type == typeof(int)) return int.Parse(text, NumberStyles.Integer, Invariant);
        if (type == typeof(long)) return long.Parse(text, NumberStyles.Integer, Invariant);
        if (type == typeof(BigInteger)) return BigInteger.Parse(text, NumberStyles.Integer, Invariant);
        if (type == typeof(double))
        {
            var number = double.Parse(text.Replace(',', '.'), NumberStyles.Float, Invariant);
            return double.IsFinite(number) ? number : throw new FormatException("Enter a finite number.");
        }
        if (type == typeof(BigRational))
            return RationalText.TryParse(text, out var rational) ? rational : throw new FormatException("Enter a decimal, fraction or scientific notation.");
        if (type == typeof(EllipticCurveQ))
            return CurveEquationText.TryParse(text, out var curve, out var error) ? curve! : throw new FormatException(error);
        if (type == typeof(FiniteFieldElement))
            return text.Split(';').Select(part => (BigRational)ParseScalar(typeof(BigRational), part)).ToArray();
        if (IsPoints(type))
            return text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                .Select(ParseRationalPoint).ToArray();
        if (IsIntegers(type))
        {
            var values = text.Split(new[] { ';', ' ', ',', '\r', '\n', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            if (type == typeof(IReadOnlyList<int>)) return values.Select(v => (int)ParseScalar(typeof(int), v)).ToArray();
            return values.Select(v => (BigInteger)ParseScalar(typeof(BigInteger), v)).ToArray();
        }
        throw new NotSupportedException("Unsupported parameter type: " + type.Name);
    }

    private static EllipticCurvePoint ParseRationalPoint(string text)
    {
        if (text.Equals("O", StringComparison.OrdinalIgnoreCase)) return EllipticCurvePoint.Infinity;
        if (text.StartsWith('(') && text.EndsWith(')')) text = text[1..^1].Trim();
        // A semicolon keeps decimal commas inside each coordinate unambiguous.
        var coordinates = text.Split(text.Contains(';') ? ';' : ',');
        if (coordinates.Length != 2)
            throw new FormatException("Use x, y, (x, y), or x; y on each line, or O. Use ; between coordinates with decimal commas.");
        return new EllipticCurvePoint((BigRational)ParseScalar(typeof(BigRational), coordinates[0]),
            (BigRational)ParseScalar(typeof(BigRational), coordinates[1]));
    }

    public static FiniteFieldElement Element(FiniteField field, string text)
    {
        var coefficients = (BigRational[])ParseScalar(typeof(FiniteFieldElement), text);
        var result = field.Zero;
        for (var i = coefficients.Length - 1; i >= 0; i--)
            result = result * field.Generator + field.CreateElement(coefficients[i].Num) / field.CreateElement(coefficients[i].Den);
        return result;
    }

    public static object? Read(Type type, string key, IReadOnlyDictionary<string, string> values, object target, CancellationToken token)
    {
        if (type == typeof(CancellationToken)) return token;
        if (type == typeof(HttpClient)) return null;
        if (IsOptions(type))
        {
            var options = Activator.CreateInstance(type)!;
            foreach (var property in type.GetProperties().Where(p => p.CanWrite))
                property.SetValue(options, Read(property.PropertyType, key + "." + property.Name, values, target, token));
            return options;
        }
        if (IsPoint(type))
        {
            var infinity = bool.Parse(values[key + ".infinity"]);
            if (type == typeof(EllipticCurvePoint))
                return infinity ? EllipticCurvePoint.Infinity : new EllipticCurvePoint(
                    (BigRational)ParseScalar(typeof(BigRational), values[key + ".x"]), (BigRational)ParseScalar(typeof(BigRational), values[key + ".y"]));
            if (target is EllipticCurveFp fp)
                return infinity ? EllipticCurvePointFp.Infinity : fp.CreatePoint(
                    (BigInteger)ParseScalar(typeof(BigInteger), values[key + ".x"]), (BigInteger)ParseScalar(typeof(BigInteger), values[key + ".y"]));
            var fq = (EllipticCurveFq)target;
            return infinity ? EllipticCurvePointFq.Infinity : fq.CreatePoint(Element(fq.Field, values[key + ".x"]), Element(fq.Field, values[key + ".y"]));
        }
        if (type == typeof(FiniteFieldElement))
            return Element(target is FiniteField field ? field : ((EllipticCurveFq)target).Field, values[key]);
        return ParseScalar(type, values[key]);
    }
}

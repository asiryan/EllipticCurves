#nullable enable
namespace EllipticCurves.Explorer.Computations;

public sealed record CalculationRequest(string OperationId, string Equation, Dictionary<string, string> Arguments,
    int TimeoutSeconds = 120, int MaxItems = 1000);

public sealed record CalculationUpdate(string Kind, string Message, double? Percent = null, string? Result = null);

public enum CalculationContext { RationalCurve, PrimeCurve, ExtensionCurve, FiniteField, Database, RationalNumber }
public enum ParameterKind { Text, Boolean, Multiline }

public sealed record CalculationParameter(string Key, string Label, string Default, string Help,
    Type ValueType, bool Advanced = false, ParameterKind Kind = ParameterKind.Text);

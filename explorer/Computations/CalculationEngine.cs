#nullable enable
using System.Numerics;
using System.Reflection;
using System.Runtime.ExceptionServices;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Computations;

/// <summary>Invoked in the calculation process only; never on the window's dispatcher.</summary>
public static class CalculationEngine
{
    public static async Task<string> ExecuteAsync(CalculationRequest request, Action<CalculationUpdate>? report = null, CancellationToken token = default)
    {
        if (request.MaxItems is < 1 or > 100_000) throw new ArgumentOutOfRangeException(nameof(request.MaxItems), "Result limit must be between 1 and 100,000.");
        var operation = CalculationCatalog.Get(request.OperationId);
        var arguments = operation.Parameters.ToDictionary(p => p.Key, p => p.Default);
        foreach (var pair in request.Arguments) arguments[pair.Key] = pair.Value;
        foreach (var parameter in operation.Parameters)
        {
            if (IsIgnoredCoordinate(parameter.Key, arguments)) continue;
            var error = CalculationInput.Validate(parameter, arguments[parameter.Key]);
            if (error.Length > 0) throw new FormatException(error);
        }
        if (!CurveEquationText.TryParse(request.Equation, out var curve, out var equationError)) throw new FormatException(equationError);
        token.ThrowIfCancellationRequested();
        report?.Invoke(new(CalculationProtocol.Progress, "Preparing curve and field"));
        var target = CreateTarget(operation.Context, curve!, arguments, token);
        report?.Invoke(new(CalculationProtocol.Progress, operation.UsesNetwork ? "Requesting LMFDB data" : "Computing · " + operation.Title));
        object? result;
        if (operation.Context == CalculationContext.Database)
        {
            var database = operation.Id == "database.import" ? LmfdbEllipticCurve.FromStoredDataJson(arguments["json"])
                : await LmfdbEllipticCurve.FetchAsync(curve!, cancellationToken: token).ConfigureAwait(false);
            result = operation.Id == "database.generators"
                ? new Dictionary<string, object?> { ["LMFDB label"] = database.Label, ["Model"] = curve, ["Generators"] = database.GetGeneratorsOnModel(curve!) }
                : database;
        }
        else if (operation.Member is MethodInfo method)
        {
            var parameters = method.GetParameters();
            // Finite-field point constructors validate membership. For this predicate, a
            // rejected affine point means false, while malformed coordinates still throw.
            object?[] inputs;
            try
            {
                inputs = parameters.Select(p => p.IsOut ? null : CalculationInput.Read(p.ParameterType, p.Name!, arguments, target, token)).ToArray();
            }
            catch (ArgumentException error) when (method.Name == "IsOnCurve" && target is EllipticCurveFp or EllipticCurveFq
                && error.ParamName == "point")
            {
                return CalculationFormatter.Format(false, request.MaxItems, report, token);
            }
            if (target is EllipticCurveQ inputCurve && method.Name is "Negate" or "Add" or "Double" or "Subtract" or "Multiply" or "IsTorsionPoint" or "TorsionOrder")
            {
                if (inputCurve.IsSingular) throw new InvalidOperationException("This operation requires a nonsingular elliptic curve.");
                foreach (var point in inputs.OfType<EllipticCurvePoint>())
                    if (!inputCurve.IsOnCurve(point)) throw new ArgumentException("Point " + point + " is not on the captured curve.");
            }
            try
            {
                if (target is EllipticCurveQ conductorCurve && method.Name == nameof(EllipticCurveQ.GetConductor))
                {
                    var conductor = conductorCurve.GetConductor(new FactorizationOptions
                    {
                        MaxDegreeOfParallelism = (int)CalculationInput.ParseScalar(typeof(int), arguments[CalculationInput.FactorizationWorkersKey])
                    }, out var factorization, token);
                    result = new ConductorCalculationResult(conductor, factorization);
                }
                else if (target is EllipticCurveQ rankCurve && method.Name == nameof(EllipticCurveQ.GetRankBounds)
                    && parameters[0].ParameterType == typeof(int))
                    result = rankCurve.GetRankBounds(new RankComputationOptions
                    {
                        SearchBound = (int)inputs[0]!, MaxSquareClasses = (int)inputs[1]!,
                        MaxDegreeOfParallelism = (int)CalculationInput.ParseScalar(typeof(int), arguments["execution.MaxDegreeOfParallelism"])
                    }, token);
                else result = method.Invoke(method.IsStatic ? null : target, inputs);
            }
            catch (TargetInvocationException error) when (error.InnerException != null)
            { ExceptionDispatchInfo.Capture(error.InnerException).Throw(); throw; }
            if (operation.Id.StartsWith("map.", StringComparison.Ordinal))
                result = MapResult(operation.Id, result!, arguments, target, token);
            else if (parameters.Any(p => p.IsOut))
            {
                var output = new Dictionary<string, object?> { ["Succeeded"] = result };
                for (var index = 0; index < parameters.Length; index++)
                    if (parameters[index].IsOut) output[CalculationOperation.Humanize(parameters[index].Name!)] = result is false ? null : inputs[index];
                result = output;
            }
        }
        else if (operation.Member is PropertyInfo property)
        {
            try { result = property.GetValue(target); }
            catch (TargetInvocationException error) when (error.InnerException != null)
            { ExceptionDispatchInfo.Capture(error.InnerException).Throw(); throw; }
        }
        else result = target is EllipticCurveQ rational ? CurveOverview(rational) : target is BigRational number
            ? new Dictionary<string, object?> { ["Reduced rational"] = number, ["Numerator"] = number.Num, ["Denominator"] = number.Den, ["Sign"] = number.Sign, ["Is zero"] = number.IsZero }
            : target;
        token.ThrowIfCancellationRequested();
        report?.Invoke(new(CalculationProtocol.Progress, "Collecting and formatting results"));
        return CalculationFormatter.Format(result, request.MaxItems, report, token);
    }

    public static bool IsIgnoredCoordinate(string key, IReadOnlyDictionary<string, string> values) =>
        (key.EndsWith(".x") || key.EndsWith(".y")) && values.TryGetValue(key[..^2] + ".infinity", out var flag)
            && bool.TryParse(flag, out var infinity) && infinity;

    private static object CreateTarget(CalculationContext context, EllipticCurveQ curve, IReadOnlyDictionary<string, string> values, CancellationToken token)
    {
        if (context is CalculationContext.RationalCurve or CalculationContext.Database) return curve;
        if (context == CalculationContext.RationalNumber) return CalculationInput.ParseScalar(typeof(BigRational), values.GetValueOrDefault("number", "0"));
        var prime = (BigInteger)CalculationInput.ParseScalar(typeof(BigInteger), values["field.prime"]);
        if (context == CalculationContext.PrimeCurve)
        {
            // A degree-one field supplies exact modular inversion of rational coefficients.
            var scalars = new FiniteField(prime, new BigInteger[] { 0, 1 }, token);
            BigInteger Reduce(BigRational value)
            {
                var element = scalars.CreateElement(value.Num) / scalars.CreateElement(value.Den);
                return element.IsZero ? BigInteger.Zero : element.Coefficients[0];
            }
            return new EllipticCurveFp(prime, Reduce(curve.A1), Reduce(curve.A2), Reduce(curve.A3), Reduce(curve.A4), Reduce(curve.A6), token);
        }
        var field = new FiniteField(prime, (BigInteger[])CalculationInput.ParseScalar(typeof(BigInteger[]), values["field.modulus"]), token);
        return context == CalculationContext.FiniteField ? field : new EllipticCurveFq(field,
            CalculationInput.Element(field, values["field.a1"]), CalculationInput.Element(field, values["field.a2"]),
            CalculationInput.Element(field, values["field.a3"]), CalculationInput.Element(field, values["field.a4"]),
            CalculationInput.Element(field, values["field.a6"]));
    }

    private static object MapResult(string id, object map, IReadOnlyDictionary<string, string> values, object target, CancellationToken token)
    {
        if (id == "map.change.inverse") return ((WeierstrassIsomorphism)map).Inverse();
        var point = (EllipticCurvePoint)CalculationInput.Read(typeof(EllipticCurvePoint), "point", values, target, token)!;
        var mapped = map switch
        {
            WeierstrassIsomorphism isomorphism => id == "map.change.backward" ? isomorphism.MapBack(point) : isomorphism.Map(point),
            RationalIsogeny isogeny => isogeny.Map(point),
            TwoIsogenyPair pair => id == "map.dual" ? pair.Dual.Map(point) : pair.Forward.Map(point),
            _ => throw new InvalidOperationException("Unknown map.")
        };
        return new Dictionary<string, object?> { ["Map"] = map, ["Input point"] = point, ["Mapped point"] = mapped };
    }

    private static object CurveOverview(EllipticCurveQ curve) => new Dictionary<string, object?>
    {
        ["Equation"] = CurveEquationText.Format(curve),
        ["a1"] = curve.A1, ["a2"] = curve.A2, ["a3"] = curve.A3, ["a4"] = curve.A4, ["a6"] = curve.A6,
        ["b2"] = curve.B2, ["b4"] = curve.B4, ["b6"] = curve.B6, ["b8"] = curve.B8,
        ["c4"] = curve.C4, ["c6"] = curve.C6, ["Discriminant"] = curve.Discriminant,
        ["j invariant"] = curve.IsSingular ? null : curve.JInvariant,
        ["Singular"] = curve.IsSingular, ["Real components"] = curve.IsSingular ? null : curve.NumberOfRealComponents,
        ["Integer coefficients"] = curve.HasIntegerCoefficients,
        ["Complex multiplication"] = curve.IsSingular ? null : curve.HasComplexMultiplication,
        ["CM discriminant (0 means non-CM)"] = curve.IsSingular ? null : curve.CmDiscriminant
    };
}

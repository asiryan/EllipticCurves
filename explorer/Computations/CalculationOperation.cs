#nullable enable
using System.Net.Http;
using System.Reflection;
using System.Text.RegularExpressions;

namespace EllipticCurves.Explorer.Computations;

public sealed class CalculationOperation
{
    public string Id { get; }
    public string Title { get; }
    public string Group { get; }
    public string Description { get; }
    public CalculationContext Context { get; }
    public MemberInfo? Member { get; }
    public IReadOnlyList<CalculationParameter> Parameters { get; }
    public bool UsesNetwork => Context == CalculationContext.Database && Id != "database.import";
    public bool UsesPlot => Context == CalculationContext.PrimeCurve
        || (Context == CalculationContext.RationalCurve && Member?.Name != "FromJInvariant")
        || (Context == CalculationContext.Database && Id != "database.import");

    public CalculationOperation(string id, string title, string group, string description, CalculationContext context,
        MemberInfo? member = null, IEnumerable<CalculationParameter>? extra = null)
    {
        Id = id; Title = title; Group = group; Description = description; Context = context; Member = member;
        var parameters = new List<CalculationParameter>();
        if (context == CalculationContext.RationalNumber && member is not MethodInfo { IsStatic: true })
            parameters.Add(new("number", "Rational number", "8.325", "Exact decimal, fraction or scientific notation.", typeof(BigRational)));
        if (context is CalculationContext.PrimeCurve or CalculationContext.ExtensionCurve or CalculationContext.FiniteField)
            parameters.Add(new("field.prime", "Prime characteristic p", "5", "A prime integer, including 2 or 3.", typeof(System.Numerics.BigInteger)));
        if (context is CalculationContext.ExtensionCurve or CalculationContext.FiniteField)
            parameters.Add(new("field.modulus", "Defining polynomial", "2; 0; 1", "Coefficients in ascending powers of t, separated by ;. Example: 2; 0; 1 means t^2 + 2. Must be irreducible modulo p.", typeof(System.Numerics.BigInteger[])));
        if (context == CalculationContext.ExtensionCurve)
            foreach (var name in new[] { "a1", "a2", "a3", "a4", "a6" })
                parameters.Add(new("field." + name, name, name == "a4" ? "-1" : "0", CalculationInput.ElementHelp, typeof(FiniteFieldElement)));
        if (member is MethodInfo method)
        {
            foreach (var parameter in method.GetParameters())
                if (!parameter.IsOut && parameter.ParameterType != typeof(CancellationToken) && parameter.ParameterType != typeof(HttpClient))
                    parameters.AddRange(CalculationInput.Describe(parameter.ParameterType, parameter.Name!,
                        parameter.HasDefaultValue ? parameter.DefaultValue : null));
            // The convenience overload has no options object. Expose the same Explorer
            // execution setting and route it through the options overload in the worker.
            if (method.DeclaringType == typeof(EllipticCurveQ) && method.Name == nameof(EllipticCurveQ.GetRankBounds)
                && method.GetParameters()[0].ParameterType == typeof(int))
                parameters.AddRange(CalculationInput.Describe(typeof(int), "execution.MaxDegreeOfParallelism",
                    CalculationInput.DefaultRankWorkers, true));
        }
        if (extra != null) parameters.AddRange(extra);
        Parameters = parameters;
    }

    public override string ToString() => Title;
    public static string Humanize(string name) => Regex.Replace(name, "(?<=[a-z0-9])(?=[A-Z])", " ");
}

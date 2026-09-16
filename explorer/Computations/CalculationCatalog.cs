#nullable enable
using System.Reflection;

namespace EllipticCurves.Explorer.Computations;

/// <summary>Public mathematical operations; object identity and duplicate property aliases are not menu actions.</summary>
public static class CalculationCatalog
{
    public static IReadOnlyList<CalculationOperation> All { get; } = Build();
    public static CalculationOperation Get(string id) => All.SingleOrDefault(item => item.Id == id)
        ?? throw new ArgumentException("Unknown calculation.");

    private static IReadOnlyList<CalculationOperation> Build()
    {
        var result = new List<CalculationOperation>
        {
            new("curve.overview", "Coefficients and invariants", "Curve and models", "Exact a, b and c coefficients, discriminant, j-invariant, real components and complex multiplication.", CalculationContext.RationalCurve),
            new("prime.overview", "Curve invariants over Fp", "Prime fields", "Reduce the entered coefficients modulo p and inspect the resulting curve. This does not change to the minimal model.", CalculationContext.PrimeCurve),
            new("extension.overview", "Curve invariants over Fq", "Extension fields", "Construct a separate curve over the specified finite extension. Its coefficients below can involve t.", CalculationContext.ExtensionCurve),
            new("field.overview", "Field and defining polynomial", "Finite-field arithmetic", "Construct the field and prove primality and irreducibility.", CalculationContext.FiniteField),
            new("number.overview", "Exact rational representation", "Rational arithmetic", "Convert a decimal or fraction to its reduced numerator and denominator, without floating-point rounding.", CalculationContext.RationalNumber),
            new("database.fetch", "Fetch LMFDB data", "LMFDB · internet", "Contact LMFDB for the captured curve: ranks, conductor, generators, torsion, heights, periods, isogenies and stored metadata. Internet access is required.", CalculationContext.Database),
            new("database.import", "Read stored LMFDB JSON", "LMFDB · internet", "Read an LMFDB data response pasted below, without contacting the network.", CalculationContext.Database,
                extra: new[] { new CalculationParameter("json", "Stored JSON", "", "Paste the complete stored LMFDB data response.", typeof(string), Kind: ParameterKind.Multiline) }),
            new("database.generators", "LMFDB generators on this model", "LMFDB · internet", "Fetch the database record and map its rational generators to the captured curve. Internet access is required.", CalculationContext.Database)
        };
        AddMethods(result, typeof(EllipticCurveQ), CalculationContext.RationalCurve);
        AddMethods(result, typeof(EllipticCurveFp), CalculationContext.PrimeCurve);
        AddMethods(result, typeof(EllipticCurveFq), CalculationContext.ExtensionCurve);
        AddMethods(result, typeof(FiniteField), CalculationContext.FiniteField);
        AddMethods(result, typeof(BigRational), CalculationContext.RationalNumber);
        foreach (var (name, title) in new[] { ("op_Addition", "Add rational numbers"), ("op_Subtraction", "Subtract rational numbers"), ("op_Multiply", "Multiply rational numbers"), ("op_Division", "Divide rational numbers") })
            result.Add(new("number." + name, title, "Rational arithmetic", "Exact rational arithmetic, independent of the plotted curve.", CalculationContext.RationalNumber,
                typeof(BigRational).GetMethod(name, new[] { typeof(BigRational), typeof(BigRational) })));
        foreach (var name in new[] { "ShortWeierstrass", "TorsionPoints", "TorsionStructure", "TamagawaProduct" })
        {
            var property = typeof(EllipticCurveQ).GetProperty(name)!;
            result.Add(new("Q." + name, Title(name), Group(name), Description(name), CalculationContext.RationalCurve, property));
        }
        // Returned maps are computations in their own right, with inputs for the point to map.
        var point = CalculationInput.Describe(typeof(EllipticCurvePoint), "point");
        var change = typeof(EllipticCurveQ).GetMethod(nameof(EllipticCurveQ.ChangeModel))!;
        result.Add(new("map.change.forward", "Map point to changed model", "Isomorphisms and isogenies", "Apply the chosen coordinate change to a rational point.", CalculationContext.RationalCurve, change, point));
        result.Add(new("map.change.backward", "Map point back from changed model", "Isomorphisms and isogenies", "The point is supplied in the changed model's coordinates; map it back to the captured curve.", CalculationContext.RationalCurve, change, point));
        result.Add(new("map.change.inverse", "Inverse coordinate change", "Isomorphisms and isogenies", "Construct the inverse isomorphism, with its source, target and exact coordinate change.", CalculationContext.RationalCurve, change));
        result.Add(new("map.minimal", "Map point to minimal model", "Isomorphisms and isogenies", "Find the minimal-model isomorphism and map the point.", CalculationContext.RationalCurve,
            typeof(EllipticCurveQ).GetMethod(nameof(EllipticCurveQ.GetMinimalModelIsomorphism)), point));
        result.Add(new("map.isogeny", "Map point through a Velu isogeny", "Isomorphisms and isogenies", "Construct the isogeny from rational kernel generators and map the supplied point.", CalculationContext.RationalCurve,
            typeof(EllipticCurveQ).GetMethod(nameof(EllipticCurveQ.CreateIsogeny)), point));
        result.Add(new("map.dual", "Map point through the dual 2-isogeny", "Isomorphisms and isogenies", "Construct a 2-isogeny. The point must be on its target; map it back using the dual.", CalculationContext.RationalCurve,
            typeof(EllipticCurveQ).GetMethod(nameof(EllipticCurveQ.CreateTwoIsogeny)), point));
        result.Add(new("map.two", "Map point through a 2-isogeny", "Isomorphisms and isogenies", "Construct the 2-isogeny and map a point from the captured curve to its target.", CalculationContext.RationalCurve,
            typeof(EllipticCurveQ).GetMethod(nameof(EllipticCurveQ.CreateTwoIsogeny)), point));
        return result.OrderBy(item => GroupOrder(item.Group)).ThenBy(item => item.Title).ToArray();
    }

    public static bool IsMathematicalMethod(MethodInfo method) => method.IsPublic && !method.IsSpecialName
        && method.DeclaringType != typeof(object) && method.Name is not ("Equals" or "GetHashCode" or "ToString");

    private static void AddMethods(List<CalculationOperation> target, Type type, CalculationContext context)
    {
        var methods = type.GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly)
            .Where(IsMathematicalMethod).ToArray();
        foreach (var method in methods)
        {
            // One Conductor action keeps its saved operation ID and exposes the
        // worker setting and factorization output through CalculationEngine.
            if (method.DeclaringType == typeof(EllipticCurveQ) && method.Name == nameof(EllipticCurveQ.GetConductor)
                && method.GetParameters()[0].ParameterType == typeof(FactorizationOptions)) continue;
            var signature = string.Join(",", method.GetParameters().Select(p => p.ParameterType.Name));
            var title = Title(method.Name);
            if (methods.Count(other => other.Name == method.Name) > 1)
            {
                var inputs = method.GetParameters().Where(p => p.ParameterType != typeof(CancellationToken)).ToArray();
                if (inputs.Length > 0)
                    title += " (" + string.Join(", ", inputs.Select(p => p.IsOut ? "with details" : p.ParameterType.Name.EndsWith("Options") ? "work limits" : CalculationOperation.Humanize(p.Name!))) + ")";
            }
            var group = context switch
            {
                CalculationContext.PrimeCurve => "Prime fields", CalculationContext.ExtensionCurve => "Extension fields",
                CalculationContext.FiniteField => "Finite-field arithmetic", _ => Group(method.Name)
            };
            if (context == CalculationContext.RationalNumber) group = "Rational arithmetic";
            var description = Description(method.Name);
            if (context == CalculationContext.RationalNumber) description = "Exact rational arithmetic, independent of the plotted curve. Square testing also returns an exact root when one exists; comparison returns -1, 0 or 1.";
            if (context == CalculationContext.PrimeCurve) description += " Uses the entered coefficients modulo p, without a minimal-model change.";
            if (context == CalculationContext.ExtensionCurve) description += " Uses the separate extension-field curve entered below.";
            target.Add(new(type.Name + "." + method.Name + "(" + signature + ")", title, group, description, context, method));
        }
    }

    private static int GroupOrder(string group) => group switch
    {
        "Curve and models" => 0, "Rational points and torsion" => 1, "Ranks and arithmetic" => 2,
        "Heights and periods" => 3, "Isomorphisms and isogenies" => 4, "Fourier coefficients and reduction" => 5,
        "Prime fields" => 6, "Extension fields" => 7, "Finite-field arithmetic" => 8, "Rational arithmetic" => 9, _ => 10
    };

    private static string Group(string name)
    {
        if (name.Contains("Height") || name is "Regulator" or "GetPeriods" or "RealEllipticLogarithm") return "Heights and periods";
        if (name.Contains("Isogen") || name.Contains("Isomorph") || name is "ChangeModel" or "GetDivisionPoints" or "TryDividePoint" or "Saturate") return "Isomorphisms and isogenies";
        if (name.Contains("Rank") || name is "GetConductor" or "GetRootNumber" or "GetLocalData" or "TamagawaProduct") return "Ranks and arithmetic";
        if (name.Contains("Fourier") || name.Contains("Frobenius") || name.Contains("Reduce") || name == "CountPoints") return "Fourier coefficients and reduction";
        if (name.Contains("Torsion") || name.Contains("Points") || name is "IsOnCurve" or "Negate" or "Add" or "Double" or "Subtract" or "Multiply") return "Rational points and torsion";
        return "Curve and models";
    }

    private static string Title(string name) => name switch
    {
        "EstimateAnalyticRank" => "Analytic rank and certification", "GetRankBounds" => "Proved rank bounds",
        "GetRankLowerBound" => "Verify rank from supplied points",
        "FaltingsHeight" => "Faltings height", "StableFaltingsHeight" => "Stable Faltings height",
        "TorsionStructure" => "Torsion group", "TorsionPoints" => "All rational torsion points",
        "RationalPoints" => "Search rational points", "IntegralPoints" => "Search integral points",
        "ShortWeierstrass" => "Short Weierstrass model", "GetGlobalMinimalModel" => "Global minimal model",
        "CreateIsogeny" => "Construct Velu isogeny", "CreateTwoIsogeny" => "Construct 2-isogeny and dual",
        "GetTwoIsogenies" => "All rational 2-isogenies", "Saturate" => "Saturate a subgroup",
        "Abs" => "Absolute value", "Pow" => "Power", "IsSquare" => "Square test and exact root", "CompareTo" => "Compare rational numbers",
        _ => CalculationOperation.Humanize(name.StartsWith("Get") ? name[3..] : name)
    };

    private static string Description(string name) => name switch
    {
        "RationalPoints" => "Search x = m/n within the chosen numerator and denominator bounds. This is a bounded search, not all rational points.",
        "IntegralPoints" => "Search integral points with |x| within the chosen bound. This does not prove that no larger integral points exist.",
        "EstimateAnalyticRank" => "Return the estimate, error diagnostics, certification status and reason. A numerical estimate is not a proof; low ranks may be certified.",
        "GetRankBounds" => "Unconditional lower and upper bounds. Matching bounds prove the rank; an unfinished descent may leave the upper bound unknown.",
        "GetRankLowerBound" => "Prove a rank lower bound from supplied rational points. Checks independence with exact arithmetic and no full descent. An inconclusive result does not prove dependence.",
        "Saturate" => "Certify saturation at the supplied primes only. The result records unresolved primes, work limits and independence certification.",
        "RealEllipticLogarithm" => "Numerical elliptic logarithm in minimal-model coordinates. This double-valued result is not a certified enclosure.",
        "Regulator" => "Certified regulator of the supplied points; this is not a claim that they form a full Mordell-Weil basis.",
        "CreateIsogeny" => "Construct a Velu isogeny from rational kernel generators. This does not search for nonrational Galois-stable kernels.",
        "FromJInvariant" => "Construct a rational curve with the supplied j-invariant. This operation does not depend on the plotted curve.",
        _ when name.Contains("Height") || name == "GetPeriods" => "Compute certified enclosures, including exact lower and upper bounds. Precision and work limits are editable below.",
        _ => "Compute " + Title(name).ToLowerInvariant() + ". The result keeps the input curve and parameters, even if you edit the plot during the calculation."
    };
}

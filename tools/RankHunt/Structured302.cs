using System.Diagnostics;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;

// Recovered sections of the published ICARM #302 family. Their identities were
// checked independently by RankStructureAudit; every specialization is checked
// again over Q. A generic rank does not guarantee that rank at every parameter.
static class Structured302
{
    static readonly JsonDocument Data = JsonDocument.Parse(typeof(Structured302).Assembly
        .GetManifestResourceStream("icarm302-sections.json")!);
    static readonly (BigRational[] X, BigRational[] Y)[] Sections = Data.RootElement.GetProperty("sections")
        .EnumerateArray().Select(p => (Coefficients(p, "x_coeffs"), Coefficients(p, "y_coeffs"))).ToArray();
    static readonly Dictionary<string, BigRational[]> Conic = Data.RootElement.GetProperty("base_change")
        .EnumerateObject().ToDictionary(p => p.Name, p => p.Value.EnumerateArray().Select(x => Q(x.GetString()!)).ToArray());

    static BigRational[] Coefficients(JsonElement p, string name) => p.GetProperty(name)
        .EnumerateArray().Select(x => Q(x.GetString()!)).ToArray();
    static BigRational Q(string s)
    {
        var parts = s.Split('/');
        return new(BigInteger.Parse(parts[0]), parts.Length == 1 ? BigInteger.One : BigInteger.Parse(parts[1]));
    }
    static BigRational Evaluate(BigRational[] coefficients, BigRational t)
    {
        BigRational value = 0;
        for (int i = coefficients.Length - 1; i >= 0; i--) value = value * t + coefficients[i];
        return value;
    }

    public static (EllipticCurveQ, EllipticCurvePoint[]) Create(BigInteger n, BigInteger d)
    {
        var t = new BigRational(n, d);
        var curve = Family302.Curve(t.Num, t.Den);
        if (curve.IsSingular) throw new InvalidOperationException("Singular family specialization.");
        var sx = new BigRational(BigInteger.Pow(t.Den, 4));
        var sy = new BigRational(BigInteger.Pow(t.Den, 6));
        var points = Sections.Select(p => new EllipticCurvePoint(Evaluate(p.X, t) * sx, Evaluate(p.Y, t) * sy)).ToArray();
        if (points.Any(p => !curve.IsOnCurve(p))) throw new InvalidOperationException("Section specialization is off the curve.");
        return (curve, points);
    }

    public static (EllipticCurveQ, EllipticCurvePoint[]) BaseChange(BigInteger n, BigInteger d)
    {
        var w = new BigRational(n, d); n = w.Num; d = w.Den;
        var denominator = 1876775*n*n + 83273802*n*d - 27990104273L*d*d;
        if (denominator.IsZero) throw new ArgumentException("This parameter lies outside the finite T chart.");
        var t = new BigRational(-423901*n*n + 39079980*n*d + 6293092321L*d*d, denominator);
        var u = new BigRational(-54322088703051L*n*n + 54314850684198L*n*d - 808951219546714491L*d*d, denominator);
        var h = new BigRational(54264759471843376L)*t*t + new BigRational(22048582273831768L)*t + new BigRational(3049459337164321L);
        if (u*u != h) throw new InvalidOperationException("Conic parametrization failed.");
        BigRational At(string name) => Evaluate(Conic[name], t);
        var l = At("l");
        if (l.IsZero) throw new ArgumentException("The extra section has a pole in this affine chart.");
        var (curve, seeds) = Create(t.Num, t.Den);
        var x = (-At("quadratic_u") + At("sqrt_factor")*u)/2;
        var parameters = Family302.Parameters(t.Num, t.Den);
        var lFamily = new BigRational(parameters[0], BigInteger.Pow(t.Den, 2));
        var bFamily = new BigRational(parameters[5], BigInteger.Pow(t.Den, 6));
        var y = (At("a")*x - At("c"))/l + (lFamily*x+bFamily)/2;
        var extra = new EllipticCurvePoint(x * new BigRational(BigInteger.Pow(t.Den, 4)),
            y * new BigRational(BigInteger.Pow(t.Den, 6)));
        if (!curve.IsOnCurve(extra)) throw new InvalidOperationException("Extra section specialization is off the curve.");
        return (curve, seeds.Append(extra).ToArray());
    }

    // Preserve existing points and scores; enrich saved grid candidates in a
    // separate directory so previous experiments remain reproducible.
    public static void Enrich(string input, string output, CancellationToken token)
    {
        if (Path.GetFullPath(input).TrimEnd(Path.DirectorySeparatorChar)
            .Equals(Path.GetFullPath(output).TrimEnd(Path.DirectorySeparatorChar), StringComparison.OrdinalIgnoreCase))
            throw new ArgumentException("Use a separate output directory for enriched candidates.");
        var candidates = JsonSerializer.Deserialize<Candidate[]>(File.ReadAllText(Path.Combine(input, "candidates.json")))!;
        Directory.CreateDirectory(output);
        var timer = Stopwatch.StartNew();
        var rows = new List<object>();
        var counts = new Dictionary<int, int>();
        foreach (var candidate in candidates.Where(c => !c.Control))
        {
            token.ThrowIfCancellationRequested();
            var watch = Stopwatch.StartNew();
            var (curve, sections) = Create(candidate.U, candidate.V);
            string name = $"curve_{candidate.U}_{candidate.V}.json";
            using var original = JsonDocument.Parse(File.ReadAllText(Path.Combine(input, name)));
            var a = original.RootElement.GetProperty("ainvs").EnumerateArray().Select(x => Q(x.GetString()!));
            if (!a.SequenceEqual(new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 }))
                throw new ArgumentException("Candidate equation differs from the section model.");
            var old = original.RootElement.GetProperty("points").EnumerateArray()
                .Select(p => new EllipticCurvePoint(Q(p[0].GetString()!), Q(p[1].GetString()!))).ToArray();
            var points = old.Concat(sections).Distinct().ToArray();
            var before = curve.GetRankLowerBound(old, 1009, token);
            var after = curve.GetRankLowerBound(points, 1009, token);
            counts[after.LowerBound] = counts.GetValueOrDefault(after.LowerBound) + 1;
            Hunt.Save(Path.Combine(output, name), Hunt.CurveData(curve, points,
                $"icarm302 recovered sections at {candidate.U}/{candidate.V}; enriched from {input}", after.LowerBound));
            rows.Add(new { candidate.U, candidate.V, before = before.LowerBound, after = after.LowerBound,
                point_count = points.Length, seconds = watch.Elapsed.TotalSeconds });
        }
        Hunt.Save(Path.Combine(output, "candidates.json"), candidates);
        Hunt.Save(Path.Combine(output, "enrichment.json"), new { input, seconds = timer.Elapsed.TotalSeconds, counts, rows });
        Console.WriteLine($"Enriched {rows.Count} candidates in {timer.Elapsed.TotalSeconds:F3}s; bounds: {JsonSerializer.Serialize(counts)}");
    }
}

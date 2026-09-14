using System.Globalization;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;
using EllipticCurves.Explorer.Models;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
using var cancellation = new CancellationTokenSource();
Console.CancelKeyPress += (_, e) => { e.Cancel = true; cancellation.Cancel(); };
string mode = args.FirstOrDefault() ?? "help";
var options = new Dictionary<string, string>();
for (int i = 1; i < args.Length; i += 2)
{
    if (i + 1 >= args.Length || !args[i].StartsWith("--")) throw new ArgumentException("Use --name value options.");
    options.Add(args[i][2..], args[i + 1]);
}
string Get(string name, string fallback) => options.GetValueOrDefault(name, fallback);
int Number(string name, int fallback, int min, int max)
{
    int n = int.Parse(Get(name, fallback.ToString()));
    if (n < min || n > max) throw new ArgumentOutOfRangeException(name);
    return n;
}
try
{
    if (mode == "grid")
    {
        var cfg = new HuntOptions(Number("height", 10000, 1, 1000000), Number("keep", 8192, 1, 100000),
            Number("refine-keep", 256, 1, 10000), Number("final-keep", 48, 1, 1000),
            Number("prime-bound", 65521, 16382, 262139), Number("workers", Math.Min(4, Environment.ProcessorCount), 1, 32));
        if (cfg.FinalKeep > cfg.RefineKeep || cfg.RefineKeep > cfg.Keep) throw new ArgumentException("Require keep >= refine-keep >= final-keep.");
        Hunt.Run(cfg, Get("output", "artifacts/rank-hunt/h10000"), cancellation.Token);
    }
    else if (mode == "export")
    {
        var u = BigInteger.Parse(Get("u", "0"));
        var v = BigInteger.Parse(Get("v", "1"));
        if (v <= 0) throw new ArgumentOutOfRangeException("v");
        string family = Get("family", "icarm302");
        var pair = family == "elkies17" ? ElkiesSearchFamily.Create(checked((int)u), checked((int)v), cancellation.Token) :
            family == "icarm302" ? (Family302.Curve(u, v), Family302.Points(u, v)) :
            family == "icarm302-17" ? Structured302.Create(u, v) :
            family == "icarm302-18" ? Structured302.BaseChange(u, v) :
            throw new ArgumentException("Unknown family");
        var cert = pair.Item1.GetRankLowerBound(pair.Item2, 1009, cancellation.Token);
        string output = Get("output", "artifacts/rank-hunt/curve.json");
        Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(output))!);
        Hunt.Save(output, Hunt.CurveData(pair.Item1, pair.Item2, $"{family} sections at {u}/{v}", cert.LowerBound));
        Console.WriteLine($"{family} {u}/{v}: certified lower bound {cert.LowerBound}");
    }
    else if (mode == "bisections")
        BisectionSearch.Certify(Get("input-dir", "artifacts/bisection-hunt"),
            Get("output", "artifacts/bisection-hunt/certified"), Number("limit",1000,1,10000), cancellation.Token);
    else if (mode == "enrich")
        Structured302.Enrich(Get("input-dir", "artifacts/rank-hunt/h10000"),
            Get("output", "artifacts/rank-structure-audit/enriched-h10000"), cancellation.Token);
    else if (mode is "verify" or "relations" or "basis")
    {
        using var data = JsonDocument.Parse(File.ReadAllText(Get("input", "")));
        var a = data.RootElement.GetProperty("ainvs").EnumerateArray().Select(x => Q(x.GetString()!)).ToArray();
        var curve = new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]);
        var points = data.RootElement.GetProperty("points").EnumerateArray()
            .Select(x => new EllipticCurvePoint(Q(x[0].GetString()!), Q(x[1].GetString()!))).ToArray();
        var cert = curve.GetRankLowerBound(points, Number("prime-bound", 1009, 3, 10000), cancellation.Token);
        if (mode == "basis")
        {
            var relations = KummerRelations.Find(curve, points, Number("prime-bound", 1009, 3, 10000));
            var nonpivots = relations.Select(r => r.Max()).ToHashSet();
            var indices = Enumerable.Range(0, points.Length).Where(i => !nonpivots.Contains(i)).ToArray();
            var selected = indices.Select(i => points[i]).ToArray();
            var selectedCert = curve.GetRankLowerBound(selected, Number("prime-bound", 1009, 3, 10000), cancellation.Token);
            // Character independence alone may include rational 2-torsion.
            // Only expose a basis for the height search after proving every
            // selected point independent modulo torsion by the lower bound.
            bool independent = selectedCert.LowerBound == selected.Length;
            Console.WriteLine(JsonSerializer.Serialize(new { point_count = points.Length,
                cert.LowerBound, cert.ImageDimension, cert.NoTwoTorsionPrime,
                selected_indices = indices, selected_lower_bound = selectedCert.LowerBound,
                all_selected_independent = independent,
                points = independent ? selected.Select(p => new[] { p.X.ToString(), p.Y.ToString() }).ToArray() : [],
                hypotheses = Array.Empty<string>() }, Hunt.JsonOptions));
            return;
        }
        if (mode == "relations")
        {
            var relations = KummerRelations.Find(curve, points, Number("prime-bound", 1009, 3, 10000));
            Hunt.Save(Get("output", "artifacts/rank-hunt/relations.json"), new
            {
                ainvs = a.Select(x => x.ToString()).ToArray(),
                points = points.Select(p => new[] { p.X.ToString(), p.Y.ToString() }).ToArray(), relations,
                source = Get("input", ""), rank_lower_bound = cert.LowerBound,
                relation_warning = "These are local character relations, not claims of rational dependence or divisibility."
            });
            Console.WriteLine($"{relations.Length} candidate relations, certified lower bound {cert.LowerBound}");
            return;
        }
        Console.WriteLine(JsonSerializer.Serialize(new { point_count = points.Length, cert.LowerBound, cert.ImageDimension,
            cert.NoTwoTorsionPrime, hypotheses = Array.Empty<string>() }, Hunt.JsonOptions));
    }
    else Console.WriteLine("RankHunt grid [--height 10000 --keep 8192 --refine-keep 256 --final-keep 48 --prime-bound 65521 --workers 4 --output directory]\nRankHunt export --family icarm302|icarm302-17|icarm302-18|elkies17 --u 0 --v 1 --output curve.json\nRankHunt enrich --input-dir candidates-directory --output enriched-directory\nRankHunt bisections --input-dir bisection-data-directory --output certified-directory --limit 1000\nRankHunt verify --input points.json");
}
catch (OperationCanceledException) { Console.WriteLine("Stopped. Completed grid batches are checkpointed; rerun with the same options."); }
catch (Exception error)
{
    Console.Error.WriteLine($"RankHunt failed: {error.Message}");
    Environment.ExitCode = 1;
}

static BigRational Q(string s)
{
    var parts = s.Split('/');
    return new(BigInteger.Parse(parts[0]), parts.Length == 1 ? BigInteger.One : BigInteger.Parse(parts[1]));
}

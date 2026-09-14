using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;

// Experimental bounded search. No score is interpreted as a rank bound.
sealed record HuntOptions(int Height, int Keep, int RefineKeep, int FinalKeep, int PrimeBound, int Workers);
sealed class Candidate
{
    public int U { get; set; }
    public int V { get; set; }
    public double Score { get; set; }
    public double ScreeningScore { get; set; }
    public double ValidationScore { get; set; }
    public bool Control { get; set; }
}
sealed class GridCheckpoint
{
    public HuntOptions Options { get; set; } = null!;
    public int NextDenominator { get; set; } = 1;
    public long PrimitiveCount { get; set; }
    public Candidate[] Retained { get; set; } = [];
}
static class Hunt
{
    const int Scale = 1 << 20;
    public static readonly JsonSerializerOptions JsonOptions = new() { WriteIndented = true };
    public static void Save(string path, object data)
    {
        File.WriteAllText(path + ".tmp", JsonSerializer.Serialize(data, JsonOptions));
        File.Move(path + ".tmp", path, true);
    }
    static int Gcd(int a, int b) { a = Math.Abs(a); while (b != 0) (a, b) = (b, a % b); return a; }

    public static void Run(HuntOptions options, string directory, CancellationToken token)
    {
        Directory.CreateDirectory(directory);
        string checkpointPath = Path.Combine(directory, "grid-checkpoint.json");
        var checkpoint = File.Exists(checkpointPath)
            ? JsonSerializer.Deserialize<GridCheckpoint>(File.ReadAllText(checkpointPath))!
            : new GridCheckpoint { Options = options };
        if (checkpoint.Options != options) throw new ArgumentException("Checkpoint options differ. Use a different output directory.");
        bool resumed = checkpoint.NextDenominator > 1;
        var queue = new PriorityQueue<Candidate, double>();
        foreach (var c in checkpoint.Retained) queue.Enqueue(c, c.Score);
        var timer = Stopwatch.StartNew();
        var tables = Local302.Primes(1021).Select(p => new Local302(p)).ToArray();
        Parallel.ForEach(tables, new ParallelOptions { MaxDegreeOfParallelism = options.Workers, CancellationToken = token },
            t => { for (int j = 0; j <= t.P; j++) t.Score(j); });
        int h = options.Height;
        var scores = new int[2 * h + 1];
        foreach (int v in Enumerable.Range(checkpoint.NextDenominator, h + 1 - checkpoint.NextDenominator))
        {
            if (token.IsCancellationRequested) break;
            Array.Clear(scores);
            foreach (var table in tables)
            {
                int p = table.P;
                if (v % p == 0)
                {
                    int add = (int)Math.Round(Scale * table.Score(p), MidpointRounding.AwayFromZero);
                    for (int j = 0; j < scores.Length; j++) scores[j] += add;
                }
                else
                {
                    int residue = table.Slot(-h, v), inverse = table.Slot(1, v);
                    var pattern = new int[p];
                    for (int j = 0; j < p; j++)
                    {
                        pattern[j] = (int)Math.Round(Scale * table.Score(residue), MidpointRounding.AwayFromZero);
                        residue += inverse; if (residue >= p) residue -= p;
                    }
                    for (int offset = 0; offset < scores.Length; offset += p)
                        for (int j = 0, n = Math.Min(p, scores.Length - offset); j < n; j++) scores[offset + j] += pattern[j];
                }
            }
            for (int j = 0; j < scores.Length; j++)
            {
                int u = j - h;
                if (Gcd(u, v) != 1) continue;
                checkpoint.PrimitiveCount++;
                double score = (double)scores[j] / Scale;
                if (queue.Count < options.Keep) queue.Enqueue(new() { U = u, V = v, Score = score }, score);
                else if (score > queue.Peek().Score)
                {
                    queue.Dequeue(); queue.Enqueue(new() { U = u, V = v, Score = score }, score);
                }
            }
            checkpoint.NextDenominator = v + 1;
            if (v % 500 == 0 || v == h)
            {
                checkpoint.Retained = queue.UnorderedItems.Select(e => e.Element).ToArray();
                Save(checkpointPath, checkpoint);
                Console.WriteLine($"grid {v}/{h}: {checkpoint.PrimitiveCount:N0} primitive parameters; this session {timer.Elapsed.TotalSeconds:F1}s");
            }
        }
        checkpoint.Retained = queue.UnorderedItems.Select(e => e.Element).ToArray();
        Save(checkpointPath, checkpoint);
        if (token.IsCancellationRequested) return;
        bool recordRetained = checkpoint.Retained.Any(c => c.U == 164518 && c.V == 924945);
        var candidates = checkpoint.Retained.Where(c => c.U != 164518 || c.V != 924945)
            .Select(c => new Candidate { U = c.U, V = c.V }).ToList();
        candidates.Add(new() { U = 164518, V = 924945, Control = true });
        // The record is evaluated separately and consumes no finalist slot.
        ScoreStage(candidates, 0, 4093, options.Workers, token);
        candidates = Select(candidates, options.RefineKeep);
        Save(Path.Combine(directory, "stage-4093.json"), candidates);
        ScoreStage(candidates, 4093, 16381, options.Workers, token);
        candidates = Select(candidates, options.FinalKeep);
        foreach (var c in candidates) c.ScreeningScore = c.Score;
        Save(Path.Combine(directory, "stage-16381.json"), candidates);
        // This set is frozen before any prime above 16381 is used.
        ScoreStage(candidates, 16381, options.PrimeBound, options.Workers, token);
        foreach (var c in candidates) c.ValidationScore = c.Score - c.ScreeningScore;
        candidates = candidates.OrderByDescending(c => c.Score).ToList();
        Save(Path.Combine(directory, "candidates.json"), candidates);
        Save(Path.Combine(directory, "run.json"), new { options, checkpoint.PrimitiveCount,
            resumed_grid = resumed, seconds_this_session = timer.Elapsed.TotalSeconds,
            record_retained_by_grid = recordRetained,
            validation_primes_start_exclusive = 16381, scores_are_rank_bounds = false,
            candidate_sections = "icarm302 recovered 17 sections; specialized bounds recomputed" });
        using var writer = new StreamWriter(Path.Combine(directory, "candidates.tsv"));
        writer.WriteLine("u\tv\tscore_16381\tfull_score\tvalidation_score\trole");
        foreach (var c in candidates)
        {
            writer.WriteLine(FormattableString.Invariant($"{c.U}\t{c.V}\t{c.ScreeningScore:G17}\t{c.Score:G17}\t{c.ValidationScore:G17}\t{(c.Control ? "published_record_control" : "candidate")}"));
            if (!c.Control)
            {
                var curve = Family302.Curve(c.U, c.V);
                if (curve.IsSingular) continue;
                var points = Structured302.Create(c.U, c.V).Item2;
                var cert = curve.GetRankLowerBound(points, 1009, token);
                Save(Path.Combine(directory, $"curve_{c.U}_{c.V}.json"), CurveData(curve, points,
                    $"icarm302 recovered 17 sections at {c.U}/{c.V}", cert.LowerBound));
            }
        }
        foreach (var c in candidates.Take(8))
            Console.WriteLine($"{c.U}/{c.V}: score={c.Score:F9}, new primes={c.ValidationScore:F9}, control={c.Control}");
    }
    public static object CurveData(EllipticCurveQ curve, IEnumerable<EllipticCurvePoint> points, string source, int lower)
        => new { ainvs = new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 }.Select(a => a.ToString()).ToArray(),
            points = points.Select(p => new[] { p.X.ToString(), p.Y.ToString() }).ToArray(), source, rank_lower_bound = lower };
    static List<Candidate> Select(List<Candidate> candidates, int count)
        => candidates.Where(c => !c.Control).OrderByDescending(c => c.Score).ThenBy(c => c.V).ThenBy(c => c.U)
            .Take(count).Concat(candidates.Where(c => c.Control)).ToList();
    public static void ScoreStage(List<Candidate> candidates, int start, int end, int workers, CancellationToken token)
    {
        var primes = Local302.Primes(end).Where(p => p > start).ToArray();
        var increments = new double[primes.Length][];
        int done = 0;
        Console.WriteLine($"scoring {candidates.Count} curves on {primes.Length} primes in ({start},{end}]");
        Parallel.For(0, primes.Length, new ParallelOptions { MaxDegreeOfParallelism = workers, CancellationToken = token }, i =>
        {
            var table = new Local302(primes[i]);
            increments[i] = candidates.Select(c => table.Score(table.Slot(c.U, c.V))).ToArray();
            int completed = Interlocked.Increment(ref done);
            if (completed % 1000 == 0) Console.WriteLine($"  {completed}/{primes.Length} primes finished");
        });
        // Deterministic summation order, irrespective of thread scheduling.
        for (int i = 0; i < primes.Length; i++)
            for (int j = 0; j < candidates.Count; j++) candidates[j].Score += increments[i][j];
    }
}

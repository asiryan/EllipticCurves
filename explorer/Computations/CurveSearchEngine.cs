#nullable enable
using System.Diagnostics;
using System.Numerics;
using System.Text.Json;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Computations;

/// <summary>Bounded family search. Local point-count scores are heuristics; only exact point certificates are rank claims.</summary>
public static class CurveSearchEngine
{
    public const string OperationId = "explorer.search-elkies";
    private sealed record PrimeTable(int Prime, double[] Scores);

    public static CurveSearchState Run(CurveSearchState initial, Action<CalculationUpdate>? report = null, CancellationToken token = default)
    {
        initial.Validate();
        var options = initial.Options;
        var state = initial;
        var stopwatch = Stopwatch.StartNew();
        var lastReport = TimeSpan.Zero;
        void Report(string message)
        {
            report?.Invoke(new(CalculationProtocol.Progress, message, 100.0 * state.NextSlot / options.Slots, JsonSerializer.Serialize(state)));
            lastReport = stopwatch.Elapsed;
        }
        Report("Preparing local point-count tables…");
        var primes = Primes(options.ScorePrimeBound);
        var tables = new PrimeTable[primes.Length];
        var parallel = new ParallelOptions { MaxDegreeOfParallelism = options.Workers, CancellationToken = token };
        Parallel.For(0, primes.Length, parallel, i => tables[i] = MakeTable(primes[i], token));
        // Saved reports are not proof inputs. Reconstruct and recheck retained candidates.
        if (state.Results.Length > 0)
        {
            var restored = state.Results.Select(c => (c.Numerator, c.Denominator)).Distinct().ToArray();
            var checkedResults = new CurveSearchCandidate[restored.Length];
            Parallel.For(0, restored.Length, parallel, i => checkedResults[i] = Inspect(restored[i].Numerator, restored[i].Denominator,
                Score(restored[i].Numerator, restored[i].Denominator, tables), options, token));
            state = state with { Candidates = checkedResults };
        }
        Report("Searching the built-in Elkies family…");
        long startTested = state.Tested;
        while (!state.Complete)
        {
            token.ThrowIfCancellationRequested();
            int count = (int)Math.Min(256, options.Slots - state.NextSlot);
            var scores = new (int A, int B, double Score, bool Valid)[count];
            long start = state.NextSlot;
            long width = (long)options.NumeratorMax - options.NumeratorMin + 1;
            Parallel.For(0, count, parallel, i =>
            {
                long slot = start + i;
                int a = options.NumeratorMin + (int)(slot % width), b = 1 + (int)(slot / width);
                if (BigInteger.GreatestCommonDivisor(a, b) != 1) return;
                scores[i] = (a, b, Score(a, b, tables), true);
            });
            var retained = state.Results.ToDictionary(c => (c.Numerator, c.Denominator));
            var best = scores.Where(s => s.Valid).Select(s => (s.A, s.B, s.Score))
                .Concat(state.Results.Select(c => (A: c.Numerator, B: c.Denominator, c.Score)))
                .OrderByDescending(s => s.Score).ThenBy(s => s.B).ThenBy(s => s.A).Take(options.KeepBest).ToArray();
            var newCandidates = best.Where(s => !retained.ContainsKey((s.A, s.B))).ToArray();
            var inspected = new CurveSearchCandidate[newCandidates.Length];
            Parallel.For(0, inspected.Length, parallel, i =>
            {
                var entry = newCandidates[i];
                inspected[i] = Inspect(entry.A, entry.B, entry.Score, options, token);
            });
            foreach (var entry in inspected) retained[(entry.Numerator, entry.Denominator)] = entry;
            state = state with
            {
                NextSlot = start + count, Tested = state.Tested + scores.Count(s => s.Valid),
                Candidates = best.Select(s => retained[(s.A, s.B)]).ToArray()
            };
            if (start == initial.NextSlot || state.Complete || stopwatch.Elapsed - lastReport >= TimeSpan.FromMilliseconds(200))
            {
                double speed = (state.Tested - startTested) / Math.Max(0.001, stopwatch.Elapsed.TotalSeconds);
                Report($"{state.Tested:N0} distinct parameters checked · {speed:N0}/s · top {state.Results.Length} candidates");
            }
        }
        return state;
    }

    private static PrimeTable MakeTable(int p, CancellationToken token)
    {
        var squares = new int[p];
        Array.Fill(squares, -1); squares[0] = 0;
        for (long x = 1; x < p; x++) squares[x * x % p] = 1;
        var s = ElkiesSearchFamily.S.Select(c => Mod(c, p)).ToArray();
        var t = ElkiesSearchFamily.T.Select(c => Mod(c, p)).ToArray();
        var scores = new double[p + 1];
        for (int parameter = 0; parameter <= p; parameter++)
        {
            token.ThrowIfCancellationRequested();
            long a = Mod(-432L * (parameter == p ? s[0] : Evaluate(s, parameter, p)), p);
            long b = Mod(432L * (parameter == p ? t[0] : Evaluate(t, parameter, p)), p);
            if (Mod(4 * (a * a % p) * a + 27 * b * b, p) == 0) continue;
            long points = p + 1;
            for (long x = 0; x < p; x++) points += squares[((x * x % p) * x + a * x + b) % p];
            scores[parameter] = Math.Log((double)points / p);
        }
        return new(p, scores);
    }

    private static double Score(int a, int b, PrimeTable[] tables)
    {
        double score = 0;
        foreach (var table in tables)
        {
            int p = table.Prime;
            int parameter = b % p == 0 ? p : (int)(Mod(a, p) * Pow(b % p, p - 2, p) % p);
            score += table.Scores[parameter];
        }
        return score;
    }

    private static CurveSearchCandidate Inspect(int a, int b, double score, CurveSearchOptions options, CancellationToken token)
    {
        var (curve, sections) = ElkiesSearchFamily.Create(a, b, token);
        if (curve.IsSingular) return new(a, b, score, CurveEquationText.Format(curve), "", 0, 0, false);
        var points = sections.Distinct().ToList();
        bool limited = false;
        var timer = Stopwatch.StartNew();
        // This first search is explicitly bounded. Small rational offsets can add
        // points, but do not replace a search on reduced 2-coverings for record work.
        var seen = new HashSet<BigRational>(points.Select(p => p.X));
        foreach (var section in sections)
        {
            if (limited) break;
            for (int denominator = 1; denominator <= 4 && !limited; denominator++)
            for (int offset = -options.ExtraSearchDepth; offset <= options.ExtraSearchDepth; offset++)
            {
                token.ThrowIfCancellationRequested();
                if (timer.Elapsed > TimeSpan.FromSeconds(2) || points.Count >= 1024) { limited = true; break; }
                if (offset == 0) continue;
                var x = section.X + new BigRational(offset, denominator * denominator);
                if (!seen.Add(x)) continue;
                if (BigRational.IsSquare(x * x * x + curve.A4 * x + curve.A6, out var y)) points.Add(new(x, y));
            }
        }
        var certificate = curve.GetRankLowerBound(points, options.CertificatePrimeBound, token);
        return new(a, b, score, CurveEquationText.Format(curve), string.Join("\n", points.Select(p => p.X + "; " + p.Y)),
            certificate.LowerBound, points.Count, limited);
    }

    private static int[] Primes(int bound)
    {
        var composite = new bool[bound + 1]; var primes = new List<int>();
        for (int p = 2; p <= bound; p++)
        {
            if (composite[p]) continue;
            if (p > 3) primes.Add(p);
            for (int n = p * 2; n <= bound; n += p) composite[n] = true;
        }
        return primes.ToArray();
    }
    private static long Evaluate(long[] coefficients, int x, int p)
    { long value = 0; foreach (var c in coefficients) value = (value * x + c) % p; return value; }
    private static long Pow(long x, int power, int p)
    { long result = 1; while (power > 0) { if ((power & 1) != 0) result = result * x % p; x = x * x % p; power >>= 1; } return result; }
    private static long Mod(long value, int p) => (value % p + p) % p;
    private static long Mod(BigInteger value, int p) => (long)((value % p + p) % p);
}

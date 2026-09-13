#nullable enable
using System.Globalization;
using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace EllipticCurves.Explorer.Models;

public sealed record CurveSearchOptions(int NumeratorMin = -50, int NumeratorMax = 50,
    int DenominatorMax = 20, int ScorePrimeBound = 127, int CertificatePrimeBound = 1009,
    int KeepBest = 8, int ExtraSearchDepth = 16, int Workers = 2)
{
    [JsonIgnore] public long Slots => ((long)NumeratorMax - NumeratorMin + 1) * DenominatorMax;
    public void Validate()
    {
        if (NumeratorMin < -1_000_000 || NumeratorMax > 1_000_000 || NumeratorMin > NumeratorMax
            || DenominatorMax is < 1 or > 100_000 || ScorePrimeBound is < 7 or > 1009
            || CertificatePrimeBound is < 3 or > 10000 || KeepBest is < 1 or > 32
            || ExtraSearchDepth is < 0 or > 1000 || Workers is < 1 or > 32)
            throw new ArgumentException("Check the search limits: numerators ±1,000,000; denominators 1–100,000; score primes 7–1009; certificate primes 3–10,000; keep 1–32; extra depth 0–1000; workers 1–32.");
    }
}

public sealed record CurveSearchCandidate(int Numerator, int Denominator, double Score,
    string Equation, string PointsText, int LowerBound, int PointCount, bool ExtraSearchLimited)
{
    [JsonIgnore] public string Parameter => Denominator == 1 ? Numerator.ToString(CultureInfo.InvariantCulture) : $"{Numerator}/{Denominator}";
    [JsonIgnore] public string RankText => "≥ " + LowerBound;
    [JsonIgnore] public string ScoreText => Score.ToString("F4", CultureInfo.InvariantCulture);
    [JsonIgnore] public string Details => $"Elkies family · t = {Parameter}\nHeuristic score: {ScoreText}\nProved rank ≥ {LowerBound}\nRational points: {PointCount}\n\n{Equation}\n\nPoints (x; y):\n{PointsText}\n\nOnly a rank lower bound is proved. Extra point search is bounded to nearby rational coordinates; failure to find more points does not give an upper bound."
        + (ExtraSearchLimited ? "\nThe additional point search reached its time or point limit." : "");
}

public sealed record CurveSearchState(CurveSearchOptions Options, long NextSlot = 0, long Tested = 0,
    CurveSearchCandidate[]? Candidates = null, string Family = ElkiesSearchFamily.Id, int Version = 1)
{
    [JsonIgnore] public CurveSearchCandidate[] Results => Candidates ?? Array.Empty<CurveSearchCandidate>();
    [JsonIgnore] public bool Complete => NextSlot == Options.Slots;
    public void Validate()
    {
        if (Options == null || Family != ElkiesSearchFamily.Id || Version != 1) throw new InvalidDataException("Unsupported search file.");
        Options.Validate();
        if (NextSlot < 0 || NextSlot > Options.Slots || Tested < 0 || Tested > NextSlot || Results.Length > Options.KeepBest)
            throw new InvalidDataException("Invalid search progress.");
        foreach (var result in Results)
            if (result == null || result.Numerator < Options.NumeratorMin || result.Numerator > Options.NumeratorMax
                || result.Denominator < 1 || result.Denominator > Options.DenominatorMax || !double.IsFinite(result.Score)
                || result.Equation == null || result.Equation.Length > 10000 || result.PointsText == null || result.PointsText.Length > 100000
                || result.PointCount < 0 || result.PointCount > 1024 || result.LowerBound < 0 || result.LowerBound > result.PointCount)
                throw new InvalidDataException("Invalid saved candidate.");
        if (Results.Select(c => (c.Numerator, c.Denominator)).Distinct().Count() != Results.Length
            || Results.Any(c => System.Numerics.BigInteger.GreatestCommonDivisor(c.Numerator, c.Denominator) != 1
                || (c.Denominator - 1L) * ((long)Options.NumeratorMax - Options.NumeratorMin + 1) + c.Numerator - Options.NumeratorMin >= NextSlot))
            throw new InvalidDataException("Saved candidates do not match the search progress.");
    }

    public static CurveSearchState Parse(string text)
    {
        if (text.Length > 1_800_000) throw new InvalidDataException("The search snapshot is too large.");
        var state = JsonSerializer.Deserialize<CurveSearchState>(text) ?? throw new InvalidDataException("Empty search snapshot.");
        state.Validate();
        return state;
    }

    public static CurveSearchState Load(string path)
    {
        if (new FileInfo(path).Length > 1_800_000) throw new InvalidDataException("The search file is too large.");
        return Parse(File.ReadAllText(path));
    }

    public void Save(string path)
    {
        Validate();
        var destination = Path.GetFullPath(path);
        var temporary = destination + "." + Guid.NewGuid().ToString("N") + ".tmp";
        try
        {
            var text = JsonSerializer.Serialize(this);
            if (text.Length > 1_800_000) throw new InvalidDataException("The search snapshot is too large.");
            File.WriteAllText(temporary, text);
            File.Move(temporary, destination, true);
        }
        finally { if (File.Exists(temporary)) File.Delete(temporary); }
    }
}

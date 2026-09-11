using System.Numerics;

namespace EllipticCurves.ConsoleApp;

internal static class CurveReport
{
    public static void Write(EllipticCurveQ curve, bool useLmfdb, TextWriter output)
    {
        WriteInvariants(curve, output);
        var native = WriteNativeArithmetic(curve, output);
        if (!useLmfdb) return;

        output.WriteLine();
        output.WriteLine("Fetching LMFDB data...");
        var stored = new LmfdbEllipticCurve(curve);
        WriteLmfdbComparison(curve, native, stored, output);
    }

    private static void WriteInvariants(EllipticCurveQ curve, TextWriter output)
    {
        output.WriteLine($"E: {curve}");
        output.WriteLine($"Short Weierstrass: {curve.ShortWeierstrass}");
        output.WriteLine($"b2 = {curve.B2}");
        output.WriteLine($"b4 = {curve.B4}");
        output.WriteLine($"b6 = {curve.B6}");
        output.WriteLine($"b8 = {curve.B8}");
        output.WriteLine($"D  = {curve.Discriminant}");
        output.WriteLine($"c4 = {curve.C4}");
        output.WriteLine($"c6 = {curve.C6}");
        output.WriteLine($"j  = {curve.JInvariant}");

        output.WriteLine();
        output.WriteLine("Computing rational torsion...");
        output.WriteLine($"Torsion: {curve.TorsionStructure}");
        output.WriteLine("Torsion points:");
        foreach (var point in curve.TorsionPoints) output.WriteLine(point);
    }

    private static NativeResults WriteNativeArithmetic(EllipticCurveQ curve, TextWriter output)
    {
        output.WriteLine();
        output.WriteLine("Computing the native minimal model and conductor...");
        var minimal = curve.GlobalMinimalModel;
        var conductor = curve.Conductor;
        output.WriteLine($"Native minimal Weierstrass model: {minimal}");
        output.WriteLine($"Native Cond(E) = {conductor}");

        output.WriteLine("Computing native rank bounds...");
        var rank = curve.GetRankBounds();
        output.WriteLine($"Native rank bounds(E) = {rank}");
        output.WriteLine($"Exact native rank proved: {rank.IsExact}");
        output.WriteLine($"Rank bounds reason: {rank.Reason}");

        output.WriteLine("Estimating the analytic rank and attempting certification...");
        var analytic = curve.EstimateAnalyticRank();
        output.WriteLine($"Native analytic rank(E) = {analytic}");
        output.WriteLine($"Root number(E) = {analytic.RootNumber}");
        output.WriteLine($"Proved rank from L-series = {analytic.ProvenRank?.ToString() ?? "Unknown"}");
        output.WriteLine($"Analytic rank reason: {analytic.Reason}");
        return new(minimal, conductor, rank, analytic);
    }

    private static void WriteLmfdbComparison(EllipticCurveQ curve, NativeResults native,
        LmfdbEllipticCurve stored, TextWriter output)
    {
        output.WriteLine($"LMFDB: {stored.Label}");
        output.WriteLine($"Url: {stored.Url}");
        output.WriteLine($"Minimal Weierstrass model: {stored.GlobalMinimalModel}");
        output.WriteLine($"Torsion: {stored.TorsionStructure}");
        output.WriteLine($"Rank(E) = {stored.Rank}");
        output.WriteLine($"Analytic rank(E) = {stored.AnalyticRank?.ToString() ?? "Unknown"}");
        output.WriteLine($"Cond(E) = {stored.Conductor}");
        output.WriteLine($"Isomorphic to E: {curve.IsIsomorphic(stored.GlobalMinimalModel)}");

        var minimalMatches = native.MinimalModel.Equals(stored.GlobalMinimalModel);
        var conductorMatches = native.Conductor == stored.Conductor;
        var rankMatches = native.Rank.LowerBound <= stored.Rank &&
            (!native.Rank.UpperBound.HasValue || stored.Rank <= native.Rank.UpperBound.Value);
        output.WriteLine($"Native minimal model matches LMFDB: {minimalMatches}");
        output.WriteLine($"Native conductor matches LMFDB: {conductorMatches}");
        output.WriteLine($"LMFDB rank is within native bounds: {rankMatches}");
        if (!minimalMatches || !conductorMatches || !rankMatches)
            throw new InvalidOperationException("Native arithmetic results do not match LMFDB.");

        // Missing estimates are unknown, not a match or a contradiction.
        var analyticMatches = native.Analytic.EstimatedRank.HasValue && stored.AnalyticRank.HasValue
            ? (native.Analytic.EstimatedRank.Value == stored.AnalyticRank.Value).ToString()
            : "Unknown";
        output.WriteLine($"Native analytic rank matches LMFDB: {analyticMatches}");
    }

    private sealed record NativeResults(EllipticCurveQ MinimalModel, BigInteger Conductor,
        RankBounds Rank, AnalyticRankResult Analytic);
}

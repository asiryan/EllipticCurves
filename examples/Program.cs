using EllipticCurves;

public static class Program
{
    public static void Main()
    {
        // Irreducibility of the Cuboid Polynomial Pa,u(t) via a Rank-Zero Elliptic Curve
        // https://arxiv.org/abs/2510.11768
        // E: Y^2 = X^3 - 17 X^2 + 72 X
        var E = new EllipticCurveQ(0, -17, 0, 72, 0);

        Console.WriteLine("E: " + E);
        Console.WriteLine($"Short Weierstrass: {E.ShortWeierstrass}");
        Console.WriteLine($"Torsion: {E.TorsionStructure}");

        Console.WriteLine($"b2 = {E.B2}");
        Console.WriteLine($"b4 = {E.B4}");
        Console.WriteLine($"b6 = {E.B6}");
        Console.WriteLine($"b8 = {E.B8}");
        Console.WriteLine($"D  = {E.Discriminant}");
        Console.WriteLine($"c4 = {E.C4}");
        Console.WriteLine($"c6 = {E.C6}");
        Console.WriteLine($"j  = {E.JInvariant}");

        Console.WriteLine("Torsion points:");
        foreach (var P in E.TorsionPoints) Console.WriteLine(P);

        // Compute via LMFDB
        var E_LMFDB = new LmfdbEllipticCurve(E);
        Console.WriteLine($"LMFDB: {E_LMFDB.Label}");
        Console.WriteLine($"Url: {E_LMFDB.Url}");
        Console.WriteLine($"Minimal Weirstrass model: {E_LMFDB.GlobalMinimalModel}");
        Console.WriteLine($"Torsion: {E_LMFDB.TorsionStructure}");
        Console.WriteLine($"Rank(E) = {E_LMFDB.Rank}");
        Console.WriteLine($"Analytic rank(E) = {E_LMFDB.AnalyticRank}");
        Console.WriteLine($"Cond(E) = {E_LMFDB.Conductor}");
        Console.WriteLine($"Isomorphic to E: {E.IsIsomorphic(E_LMFDB.GlobalMinimalModel)}");

        // Compute the same invariants locally and check them against LMFDB
        var nativeMinimal = E.GlobalMinimalModel;
        var nativeRank = E.GetRankBounds();
        var nativeConductor = E.Conductor;
        Console.WriteLine($"Native minimal Weierstrass model: {nativeMinimal}");
        Console.WriteLine($"Native rank bounds(E) = {nativeRank}");
        Console.WriteLine($"Exact native rank proved: {nativeRank.IsExact}");
        Console.WriteLine($"Native Cond(E) = {nativeConductor}");

        var minimalMatches = nativeMinimal.Equals(E_LMFDB.GlobalMinimalModel);
        var conductorMatches = nativeConductor == E_LMFDB.Conductor;
        var rankMatches = nativeRank.LowerBound <= E_LMFDB.Rank && (!nativeRank.UpperBound.HasValue || E_LMFDB.Rank <= nativeRank.UpperBound.Value);
        Console.WriteLine($"Native minimal model matches LMFDB: {minimalMatches}");
        Console.WriteLine($"Native conductor matches LMFDB: {conductorMatches}");
        Console.WriteLine($"LMFDB rank is within native bounds: {rankMatches}");

        if (!minimalMatches || !conductorMatches || !rankMatches) Console.WriteLine("Native arithmetic results do not match LMFDB.");

        // Estimate the analytic rank locally and attempt a rigorous rank 0/1 certificate
        var analytic = E.EstimateAnalyticRank();
        Console.WriteLine($"Native analytic rank(E) = {analytic}");
        Console.WriteLine($"Root number(E) = {analytic.RootNumber}");
        Console.WriteLine($"Proved rank from L-series = {analytic.ProvenRank?.ToString() ?? "Unknown"}");
        var analyticMatches = analytic.EstimatedRank == E_LMFDB.AnalyticRank;
        Console.WriteLine($"Native analytic rank matches LMFDB: {analyticMatches}");
    }
}

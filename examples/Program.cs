using EllipticCurves;

public static class Program
{
    public static void Main()
    {
        // Y^2 = X^3 - 17 X^2 + 72 X
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

        var E_LMFDB = new EllipticCurveLMFDB(E);
        Console.WriteLine($"LMFDB: {E_LMFDB.Label}");
        Console.WriteLine($"Url: {E_LMFDB.Url}");
        Console.WriteLine($"Minimal Weirstrass model: {E_LMFDB.GlobalMinimalModel}");
        Console.WriteLine($"Torsion: {E_LMFDB.TorsionStructure}");
        Console.WriteLine($"Rank(E) = {E_LMFDB.Rank}");
        Console.WriteLine($"Analytic rank(E) = {E_LMFDB.AnalyticRank}");
        Console.WriteLine($"Cond(E) = {E_LMFDB.Conductor}");
        Console.WriteLine($"Isomorphic to E: {E.IsIsomorphic(E_LMFDB.GlobalMinimalModel)}");
    }
}
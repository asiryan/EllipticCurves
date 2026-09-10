namespace EllipticCurves.Visualizer.Models;

public sealed record CurvePreset(string Name, string Description, int A1, int A2, int A3, int A4, int A6)
{
    public override string ToString() => Name;

    public static IReadOnlyList<CurvePreset> All { get; } = new[]
    {
        new CurvePreset("The classic", "Two real components", 0, 0, 0, -1, 0),
        new CurvePreset("37.a1", "A general Weierstrass model", 0, 0, 1, -1, 0),
        new CurvePreset("48.a3", "A curve with eight torsion points", 0, -17, 0, 72, 0),
        new CurvePreset("The cusp", "A singular cubic: Δ = 0", 0, 0, 0, 0, 0)
    };
}

namespace EllipticCurves.Explorer.Models;

public sealed record CurvePreset(string Name, string Description, int A1, int A2, int A3, int A4, int A6)
{
    public const string ClassicEquation = "y^2 = x^3 - x";
    public static CurvePreset Classic { get; } = new("The classic", "Two real components", 0, 0, 0, -1, 0);
    public EllipticCurveQ CreateCurve() => new(A1, A2, A3, A4, A6);
    public override string ToString() => Name;

    public static IReadOnlyList<CurvePreset> All { get; } = new[]
    {
        Classic,
        new CurvePreset("The cusp", "A singular cubic: Δ = 0", 0, 0, 0, 0, 0),
        new CurvePreset("37.a1", "A general Weierstrass model", 0, 0, 1, -1, 0),
        new CurvePreset("48.a3", "A curve with eight torsion points", 0, -17, 0, 72, 0)
    };
}

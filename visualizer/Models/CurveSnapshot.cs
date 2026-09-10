namespace EllipticCurves.Visualizer.Models;

/// <summary>Exact, inexpensive library results for one immutable set of coefficients.</summary>
public sealed class CurveSnapshot
{
    public EllipticCurveQ Curve { get; }
    public CurvePlotData Plot { get; }
    public string Equation { get; }
    public string Discriminant { get; }
    public string JInvariant { get; }
    public string Components { get; }
    public string Status { get; }
    public string ComponentNote { get; }
    public string ShortEquation { get; }
    public string C4 { get; }
    public string C6 { get; }
    public bool IsSingular { get; }

    public CurveSnapshot(EllipticCurveQ curve)
    {
        Curve = curve;
        IsSingular = curve.IsSingular;
        Equation = Pretty(curve.ToString());
        Discriminant = curve.Discriminant.ToString();
        JInvariant = IsSingular ? "Undefined" : curve.JInvariant.ToString();
        Components = IsSingular ? "—" : curve.NumberOfRealComponents.ToString();
        Status = IsSingular ? "Singular cubic" : "Smooth curve";
        ComponentNote = IsSingular ? "Not an elliptic curve" : "Connected components of E(ℝ)";
        ShortEquation = Pretty(curve.ShortWeierstrass.ToString());
        C4 = curve.C4.ToString();
        C6 = curve.C6.ToString();
        Plot = new CurvePlotData(curve);
    }

    public string Summary => $"{Equation}\nΔ = {Discriminant}\nj = {JInvariant}\nc₄ = {C4}\nc₆ = {C6}\nReal components = {Components}\nShort model: {ShortEquation}";

    private static string Pretty(string value) => value.Replace("^2", "²").Replace("^3", "³").Replace("*", " ").Replace("-", "−");
}

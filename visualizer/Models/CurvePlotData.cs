namespace EllipticCurves.Visualizer.Models;

/// <summary>Floating-point coordinates for the real locus; never used for arithmetic certificates.</summary>
public sealed class CurvePlotData
{
    private readonly double a1, a3, b, c, d;
    public IReadOnlyList<double> Roots { get; }

    public CurvePlotData(EllipticCurveQ curve)
    {
        a1 = ToDouble(curve.A1);
        a3 = ToDouble(curve.A3);
        // Completing the square gives (y + (a1*x+a3)/2)^2 = x^3+b*x^2+c*x+d.
        b = ToDouble(curve.B2) / 4;
        c = ToDouble(curve.B4) / 2;
        d = ToDouble(curve.B6) / 4;
        Roots = FindRoots();
    }

    public bool TryEvaluate(double x, out double upper, out double lower)
    {
        var square = Polynomial(x);
        var error = 1e-13 * (1 + Math.Abs(x * x * x) + Math.Abs(b * x * x) + Math.Abs(c * x) + Math.Abs(d));
        if (!double.IsFinite(square) || square < -error)
        {
            upper = lower = double.NaN;
            return false;
        }
        var radius = Math.Sqrt(Math.Max(0, square));
        var center = -(a1 * x + a3) / 2;
        upper = center + radius;
        lower = center - radius;
        return double.IsFinite(upper) && double.IsFinite(lower);
    }

    public double CenterY(double x) => -(a1 * x + a3) / 2;

    public static double ToDouble(BigRational value) => (double)value.Num / (double)value.Den;

    private double Polynomial(double x) => ((x + b) * x + c) * x + d;

    private double[] FindRoots()
    {
        // Isolate roots between derivative zeros, then bisect monotone intervals.
        // Including stationary zeros also retains isolated real points of singular cubics.
        var bound = 1 + Math.Max(Math.Abs(b), Math.Max(Math.Abs(c), Math.Abs(d)));
        var cuts = new List<double> { -bound };
        var derivative = b * b - 3 * c;
        if (derivative >= 0)
        {
            var radius = Math.Sqrt(derivative);
            cuts.Add((-b - radius) / 3);
            if (radius > 0) cuts.Add((-b + radius) / 3);
        }
        cuts.Add(bound);
        cuts.Sort();
        var roots = new List<double>();
        foreach (var cut in cuts)
        {
            var tolerance = 1e-13 * (1 + Math.Abs(cut * cut * cut) + Math.Abs(b * cut * cut) + Math.Abs(c * cut) + Math.Abs(d));
            if (Math.Abs(Polynomial(cut)) <= tolerance) roots.Add(cut);
        }
        for (var i = 1; i < cuts.Count; i++)
        {
            var left = cuts[i - 1];
            var right = cuts[i];
            var sign = Math.Sign(Polynomial(left));
            if (sign == 0 || Math.Sign(Polynomial(right)) == 0 || sign == Math.Sign(Polynomial(right))) continue;
            for (var iteration = 0; iteration < 80; iteration++)
            {
                var middle = left + (right - left) / 2;
                if (Math.Sign(Polynomial(middle)) == sign) left = middle;
                else right = middle;
            }
            roots.Add(left + (right - left) / 2);
        }
        roots.Sort();
        var distinct = new List<double>();
        foreach (var root in roots)
            if (distinct.Count == 0 || Math.Abs(root - distinct[^1]) > 1e-9 * (1 + Math.Abs(root))) distinct.Add(root);
        return distinct.ToArray();
    }
}

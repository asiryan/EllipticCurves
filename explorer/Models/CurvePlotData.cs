namespace EllipticCurves.Explorer.Models;

/// <summary>Floating-point coordinates for the real locus; never used for arithmetic certificates.</summary>
public sealed class CurvePlotData
{
    internal const double MaximumPlotMagnitude = 1e300;
    private readonly double a1, a3, b, c, d;
    public IReadOnlyList<double> Roots { get; }
    public bool IsDrawable { get; }
    public double CharacteristicScale { get; }

    public CurvePlotData(EllipticCurveQ curve)
    {
        a1 = ToDouble(curve.A1);
        a3 = ToDouble(curve.A3);
        // Completing the square gives (y + (a1*x+a3)/2)^2 = x^3+b*x^2+c*x+d.
        b = ToDouble(curve.B2 / 4);
        c = ToDouble(curve.B4 / 2);
        d = ToDouble(curve.B6 / 4);
        CharacteristicScale = Math.Max(Math.Abs(b), Math.Max(Math.Sqrt(Math.Abs(c)), Math.Cbrt(Math.Abs(d))));
        var exact = new[] { curve.A1, curve.A3, curve.B2 / 4, curve.B4 / 2, curve.B6 / 4 };
        IsDrawable = exact.All(value => double.IsFinite(ToDouble(value)) && (value.IsZero || ToDouble(value) != 0));
        Roots = IsDrawable ? FindRoots() : Array.Empty<double>();
        IsDrawable = IsDrawable && Roots.Count > 0 && Roots.All(root =>
            double.IsFinite(root) && Math.Abs(root) <= MaximumPlotMagnitude
            && double.IsFinite(CenterY(root)) && Math.Abs(CenterY(root)) <= MaximumPlotMagnitude);
    }

    public bool TryEvaluate(double x, out double upper, out double lower)
    {
        if (!IsDrawable || !double.IsFinite(x)) { upper = lower = double.NaN; return false; }
        // Evaluate after rescaling x. An absolute error floor would accept entire
        // negative intervals of small curves; unscaled cubes can also underflow.
        var scale = Math.Max(Math.Abs(x), CharacteristicScale);
        var z = scale == 0 ? 0 : x / scale;
        var nb = scale == 0 ? 0 : b / scale;
        var nc = scale == 0 ? 0 : c / scale / scale;
        var nd = scale == 0 ? 0 : d / scale / scale / scale;
        var square = ((z + nb) * z + nc) * z + nd;
        var error = 1e-13 * (Math.Abs(z * z * z) + Math.Abs(nb * z * z) + Math.Abs(nc * z) + Math.Abs(nd));
        if (!double.IsFinite(square) || !double.IsFinite(error) || square < -error)
        {
            upper = lower = double.NaN;
            return false;
        }
        var radius = Math.Sqrt(Math.Max(0, square)) * Math.Sqrt(scale) * scale;
        var center = CenterY(x);
        upper = center + radius;
        lower = center - radius;
        return double.IsFinite(upper) && double.IsFinite(lower);
    }

    public double CenterY(double x) => -(a1 / 2 * x + a3 / 2);

    public static double ToDouble(BigRational value)
    {
        // Large numerator and denominator can have a small ratio; avoid Infinity / Infinity.
        var numerator = System.Numerics.BigInteger.Abs(value.Num);
        var numeratorShift = (int)Math.Max(0, numerator.GetBitLength() - 54);
        var denominatorShift = (int)Math.Max(0, value.Den.GetBitLength() - 54);
        return value.Sign * Math.ScaleB((double)(numerator >> numeratorShift) / (double)(value.Den >> denominatorShift),
            numeratorShift - denominatorShift);
    }

    private double[] FindRoots()
    {
        // Isolate roots between derivative zeros, then bisect monotone intervals.
        // Including stationary zeros also retains isolated real points of singular cubics.
        // Normalize the x scale before isolating roots, so large and tiny coefficients
        // do not overflow the derivative or collapse roots under an absolute tolerance.
        var scale = CharacteristicScale;
        if (scale == 0) return new[] { 0.0 };
        var nb = b / scale;
        var nc = c / scale / scale;
        var nd = d / scale / scale / scale;
        double P(double x) => ((x + nb) * x + nc) * x + nd;
        var bound = 1 + Math.Max(Math.Abs(nb), Math.Max(Math.Abs(nc), Math.Abs(nd)));
        var cuts = new List<double> { -bound };
        var derivative = nb * nb - 3 * nc;
        if (derivative >= 0)
        {
            var radius = Math.Sqrt(derivative);
            cuts.Add((-nb - radius) / 3);
            if (radius > 0) cuts.Add((-nb + radius) / 3);
        }
        cuts.Add(bound);
        cuts.Sort();
        var roots = new List<double>();
        foreach (var cut in cuts)
        {
            var tolerance = 1e-13 * (1 + Math.Abs(cut * cut * cut) + Math.Abs(nb * cut * cut) + Math.Abs(nc * cut) + Math.Abs(nd));
            if (Math.Abs(P(cut)) <= tolerance) roots.Add(cut);
        }
        for (var i = 1; i < cuts.Count; i++)
        {
            var left = cuts[i - 1];
            var right = cuts[i];
            var sign = Math.Sign(P(left));
            if (sign == 0 || Math.Sign(P(right)) == 0 || sign == Math.Sign(P(right))) continue;
            for (var iteration = 0; iteration < 80; iteration++)
            {
                var middle = left + (right - left) / 2;
                if (Math.Sign(P(middle)) == sign) left = middle;
                else right = middle;
            }
            roots.Add(left + (right - left) / 2);
        }
        roots.Sort();
        var distinct = new List<double>();
        foreach (var root in roots)
            if (distinct.Count == 0 || Math.Abs(root - distinct[^1]) > 1e-9 * (1 + Math.Abs(root))) distinct.Add(root);
        return distinct.Select(root => root * scale).ToArray();
    }
}

#nullable enable
using System.Globalization;

namespace EllipticCurves.Explorer.Models;

/// <summary>Numerical coordinates in the minimal model's period basis. No drawing geometry is used for arithmetic.</summary>
public sealed class TorusLattice
{
    public const int MaximumSamples = 32;
    public EllipticCurveQ Curve { get; }
    public PeriodResult Periods { get; }
    public double Omega1 { get; }
    public double Omega2Real { get; }
    public double Omega2Imaginary { get; }
    public double TauReal => Omega2Real / Omega1;
    public double TauImaginary => Omega2Imaginary / Omega1;
    public string BasisText => $"ω₁ ≈ {Number(Omega1)}    ω₂ ≈ {Number(Omega2Real)} + {Number(Omega2Imaginary)}i";
    public string TauText => $"τ ≈ {Number(TauReal)} + {Number(TauImaginary)}i";
    public string ModelNote => "Periods refer to the global minimal model: " + Periods.MinimalModel;

    private TorusLattice(EllipticCurveQ curve, PeriodResult periods)
    {
        Curve = curve;
        Periods = periods;
        Omega1 = periods.PrimitiveRealPeriod.Approximation;
        Omega2Real = periods.SecondPeriodRealPart.Approximation;
        Omega2Imaginary = periods.SecondPeriodImaginaryPart.Approximation;
        if (!double.IsFinite(Omega1) || Omega1 <= 0 || !double.IsFinite(Omega2Real)
            || !double.IsFinite(Omega2Imaginary) || Omega2Imaginary <= 0
            || !double.IsFinite(TauImaginary) || TauImaginary <= 0)
            throw new ArithmeticException("The period basis exceeds the view's numerical range.");
        foreach (var period in new[] { periods.PrimitiveRealPeriod, periods.SecondPeriodImaginaryPart })
            if (period.LowerBound <= 0 || period.Width > period.LowerBound / 10000000000L)
                throw new ArithmeticException("The period basis needs more precision than this interactive view provides.");
    }

    private static RealComputationOptions Options() => new()
    {
        DecimalDigits = 12,
        PrecisionBits = 256,
        MaxIterations = 256,
        MaxRootWork = 30000
    };

    public static TorusLattice Create(EllipticCurveQ curve, CancellationToken token)
    {
        token.ThrowIfCancellationRequested();
        if (curve.IsSingular) throw new ArgumentException("A singular cubic has no smooth complex torus.", nameof(curve));
        return new TorusLattice(curve, curve.GetPeriods(Options(), token));
    }

    public TorusPointSet MapPoints(IReadOnlyList<EllipticCurvePoint> samples, CancellationToken token)
    {
        var points = new List<TorusPoint> { TorusPoint.Origin };
        var skipped = 0;
        var candidates = samples.Where(point => !point.IsInfinity).Distinct().ToArray();
        foreach (var point in candidates.Take(MaximumSamples))
        {
            token.ThrowIfCancellationRequested();
            try
            {
                var log = Curve.RealEllipticLogarithm(point, Options(), token);
                if (Math.Abs(log.PrimitiveRealPeriod / Omega1 - 1) > 1e-9)
                    throw new ArithmeticException("The point and the view use different period bases.");
                var coordinates = TorusCoordinates.FromLogarithm(log.RealPart, log.ImaginaryPart,
                    Omega1, Omega2Real, Omega2Imaginary);
                points.Add(new TorusPoint($"P{points.Count}", point, coordinates));
            }
            catch (Exception error) when (error is ArithmeticException or ArgumentException or InvalidOperationException)
            {
                // A numerical logarithm can fail even when the exact rational point is valid.
                skipped++;
            }
        }
        return new TorusPointSet(points.ToArray(), candidates.Length, skipped);
    }

    internal static string Number(double value) => value.ToString("G7", CultureInfo.InvariantCulture);
}

public readonly record struct TorusCoordinates(double U, double V)
{
    public static TorusCoordinates FromLogarithm(double real, double imaginary,
        double omega1, double omega2Real, double omega2Imaginary)
    {
        if (!double.IsFinite(real) || !double.IsFinite(imaginary) || !double.IsFinite(omega1) || omega1 <= 0
            || !double.IsFinite(omega2Real) || !double.IsFinite(omega2Imaginary) || omega2Imaginary <= 0)
            throw new ArithmeticException("Invalid period coordinates.");
        var v = imaginary / omega2Imaginary;
        // Subtract the shear before reducing modulo 1: omega2 need not be purely imaginary.
        var u = real / omega1 - v * (omega2Real / omega1);
        if (!double.IsFinite(u) || !double.IsFinite(v)) throw new ArithmeticException("Period coordinates overflowed.");
        return new TorusCoordinates(Unit(u), Unit(v));
    }

    private static double Unit(double value) => value - Math.Floor(value);
}

public sealed record TorusPoint(string Name, EllipticCurvePoint Point, TorusCoordinates Coordinates)
{
    public static TorusPoint Origin { get; } = new("O", EllipticCurvePoint.Infinity, new(0, 0));
    public string DisplayName => Point.IsInfinity ? "O · point at infinity" : $"{Name} · {Point}";
    public string CoordinateText => $"z ≡ uω₁ + vω₂    u ≈ {TorusLattice.Number(Coordinates.U)}    v ≈ {TorusLattice.Number(Coordinates.V)}";
}

public sealed record TorusPointSet(IReadOnlyList<TorusPoint> Points, int SampleCount, int SkippedCount)
{
    public string Summary => $"{Points.Count - 1} / {SampleCount} rational samples mapped · O included"
        + (SampleCount > TorusLattice.MaximumSamples ? $" · first {TorusLattice.MaximumSamples} considered" : "")
        + (SkippedCount > 0 ? $" · {SkippedCount} numerical mappings unavailable" : "");
}

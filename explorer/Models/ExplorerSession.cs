#nullable enable
using System.IO;
using System.Text.Json;
using System.Text.Json.Serialization;
using EllipticCurves.Explorer.Computations;

namespace EllipticCurves.Explorer.Models;

public sealed record PlotViewState(double CenterX, double CenterY, double VerticalSpan)
{
    public static PlotViewState Default { get; } = new(0.3, 0, 3.4);
}
public sealed record TorusCameraState(double Azimuth, double Elevation, double Span)
{
    public static TorusCameraState Default { get; } = new(35, 32, 5.8);
}
public sealed record SidebarSession(bool Visible, double Width);
public sealed record CalculationSession(CalculationRequest Request, DateTime StartedAt, string Status,
    string Stage, TimeSpan Elapsed, double? Percent, string Result);

public sealed record ExplorerSession
{
    [JsonRequired] public string Format { get; init; } = "EllipticCurves.Explorer.Session";
    [JsonRequired] public int Version { get; init; } = 1;
    public required string Equation { get; init; }
    public required string SliderStep { get; init; }
    public string? Preset { get; init; }
    public int[] SliderOffsets { get; init; } = Array.Empty<int>();
    public bool ShowGrid { get; init; } = true;
    public bool ShowPoints { get; init; } = true;
    public bool ComplexView { get; init; }
    public bool FitRealViewWhenShown { get; init; }
    public bool CoefficientsExpanded { get; init; }
    public double EquationScrollOffset { get; init; }
    public double TorusScrollOffset { get; init; }
    public string? SelectedTorusPoint { get; init; }
    public required PlotViewState Plot { get; init; }
    public required TorusCameraState TorusCamera { get; init; }
    public required SidebarSession EquationPanel { get; init; }
    public required SidebarSession ResultsPanel { get; init; }
    public required List<CalculationSession> History { get; init; }

    public static ExplorerSession New() => new()
    {
        Equation = "y^2 = x^3 - x", SliderStep = "0.01", Preset = "The classic", SliderOffsets = new[] { 0, 0 },
        Plot = PlotViewState.Default, FitRealViewWhenShown = true, TorusCamera = TorusCameraState.Default,
        SelectedTorusPoint = EllipticCurvePoint.Infinity.ToString(),
        EquationPanel = new(true, 238), ResultsPanel = new(true, 300), History = new()
    };
}

public static class SessionFile
{
    private const long MaxFileBytes = 256L * 1024 * 1024;
    private static readonly JsonSerializerOptions Options = new() { WriteIndented = true, MaxDepth = 32 };

    public static ExplorerSession Load(string path)
    {
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        if (stream.Length > MaxFileBytes) throw new InvalidDataException("The session file exceeds 256 MB.");
        try
        {
            var session = JsonSerializer.Deserialize<ExplorerSession>(stream, Options)
                ?? throw new InvalidDataException("The session file is empty.");
            Validate(session);
            return session;
        }
        catch (JsonException error) { throw new InvalidDataException("This is not a valid Explorer session file.", error); }
    }

    public static void Save(string path, ExplorerSession session)
    {
        Validate(session);
        var destination = Path.GetFullPath(path);
        var temporary = Path.Combine(Path.GetDirectoryName(destination)!, ".ec-session-" + Guid.NewGuid().ToString("N") + ".tmp");
        try
        {
            using (var stream = new FileStream(temporary, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            {
                JsonSerializer.Serialize(stream, session, Options);
                if (stream.Length > MaxFileBytes) throw new InvalidDataException("The session file exceeds 256 MB.");
                stream.Flush(flushToDisk: true);
            }
            // Replace only after the complete snapshot has been written successfully.
            File.Move(temporary, destination, overwrite: true);
        }
        finally { if (File.Exists(temporary)) File.Delete(temporary); }
    }

    public static void Validate(ExplorerSession session)
    {
        if (session.Format != "EllipticCurves.Explorer.Session" || session.Version != 1)
            throw new InvalidDataException("This session format or version is not supported.");
        if (session.Equation == null || !CurveEquationText.TryParse(session.Equation, out var curve, out _))
            throw new InvalidDataException("The session contains an invalid curve equation.");
        if (session.SliderStep == null || !RationalText.TryParse(session.SliderStep, out var step) || step.Sign <= 0)
            throw new InvalidDataException("The session contains an invalid slider step.");
        var count = curve!.A1.IsZero && curve.A2.IsZero && curve.A3.IsZero ? 2 : 5;
        if (session.SliderOffsets == null || (session.SliderOffsets.Length != 0 && session.SliderOffsets.Length != count)
            || session.SliderOffsets.Any(value => value is < -50 or > 50))
            throw new InvalidDataException("The session contains invalid slider positions.");
        if (session.Plot is not { } plot || !Finite(plot.CenterX, plot.CenterY, plot.VerticalSpan)
            || plot.VerticalSpan is < 1e-12 or > 1e300 || Math.Abs(plot.CenterX) > 1e300 || Math.Abs(plot.CenterY) > 1e300
            || session.TorusCamera is not { } camera || !Finite(camera.Azimuth, camera.Elevation, camera.Span)
            || camera.Elevation is < -80 or > 80 || camera.Span is < 3.5 or > 18)
            throw new InvalidDataException("The session contains an invalid graph view.");
        if (!PanelValid(session.EquationPanel, 238) || !PanelValid(session.ResultsPanel, 300)
            || !Finite(session.EquationScrollOffset, session.TorusScrollOffset)
            || session.EquationScrollOffset < 0 || session.TorusScrollOffset < 0)
            throw new InvalidDataException("The session contains an invalid panel layout.");
        if (session.SelectedTorusPoint?.Length > 100_000)
            throw new InvalidDataException("The selected point is too large.");
        if (session.History == null || session.History.Count > 50)
            throw new InvalidDataException("The session contains an invalid calculation history.");
        foreach (var job in session.History)
        {
            if (job?.Request is not { } request || request.OperationId == null
                || !CalculationCatalog.All.Any(operation => operation.Id == request.OperationId)
                || request.Equation == null || request.Equation.Length > 100_000
                || request.Arguments == null || request.Arguments.Count > 512
                || request.Arguments.Any(pair => pair.Key.Length > 1024 || pair.Value == null || pair.Value.Length > 100_000)
                || request.TimeoutSeconds is < 0 or > 86_400 || request.MaxItems is < 1 or > 100_000
                || job.Status is not ("Running" or "Completed" or "Cancelled" or "Timed out" or "Failed" or "Interrupted")
                || job.Stage == null || job.Stage.Length > 100_000 || job.Result == null || job.Result.Length > 2_100_000
                || job.Elapsed < TimeSpan.Zero || job.Elapsed.TotalDays > 3650
                || (job.Percent is { } percent && (!double.IsFinite(percent) || percent is < 0 or > 100)))
                throw new InvalidDataException("The session contains an invalid or unsupported calculation.");
        }
    }

    private static bool Finite(params double[] values) => values.All(double.IsFinite);
    private static bool PanelValid(SidebarSession? panel, double minimum) =>
        panel != null && double.IsFinite(panel.Width) && panel.Width >= minimum && panel.Width <= 100_000;
}

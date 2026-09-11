#nullable enable
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
    public const int HistoryLimit = 50;
    public const string DefaultSliderStep = "0.01";
    [JsonRequired] public string Format { get; init; } = SessionFile.FormatName;
    [JsonRequired] public int Version { get; init; } = SessionFile.CurrentVersion;
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

    // Files restore data into a fresh workspace. Presentation belongs to the
    // current window and its undo history, not to the document.
    internal ExplorerSession DataOnly() => New() with
    {
        Format = Format, Version = Version, Equation = Equation, History = History,
        Preset = null, SliderOffsets = Array.Empty<int>()
    };

    public static ExplorerSession New() => new()
    {
        Equation = CurvePreset.ClassicEquation, SliderStep = DefaultSliderStep, Preset = CurvePreset.Classic.Name, SliderOffsets = new[] { 0, 0 },
        Plot = PlotViewState.Default, FitRealViewWhenShown = true, TorusCamera = TorusCameraState.Default,
        SelectedTorusPoint = EllipticCurvePoint.Infinity.ToString(),
        EquationPanel = new(true, 238), ResultsPanel = new(true, 300), History = new()
    };
}

#nullable enable
using System.IO;
using System.Globalization;
using System.Text.Json;
using EllipticCurves.Explorer.Computations;

namespace EllipticCurves.Explorer.Models;

public static class SessionFile
{
    public const string Extension = ".ec";
    public const string DefaultFileName = "untitled" + Extension;
    public const string DialogFilter = "Explorer session (*" + Extension + ")|*" + Extension;
    public const string FormatName = "EllipticCurves.Explorer.Session";
    public const int CurrentVersion = 1;
    private const long MaxFileBytes = 256L * 1024 * 1024;
    private const string FileTooLargeMessage = "The session file exceeds 256 MB.";
    private static readonly JsonSerializerOptions Options = new() { WriteIndented = true, MaxDepth = 32 };

    public static string GetFileName(string? path) => path == null ? DefaultFileName : Path.GetFileName(path);

    public static bool IsValidFileName(string name)
    {
        name = name.Trim();
        return name.Length > 0 && !name.EndsWith('.') && name.IndexOfAny(Path.GetInvalidFileNameChars()) < 0;
    }

    public static string? SuggestSavePath(string? currentPath, string? fileName)
    {
        if (string.IsNullOrWhiteSpace(fileName)) return currentPath;
        fileName = fileName.Trim();
        if (!fileName.EndsWith(Extension, StringComparison.OrdinalIgnoreCase)) fileName += Extension;
        return currentPath == null ? fileName : Path.Combine(Path.GetDirectoryName(currentPath)!, fileName);
    }

    public static ExplorerSession Load(string path)
    {
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        if (stream.Length > MaxFileBytes) throw new InvalidDataException(FileTooLargeMessage);
        try
        {
            var session = JsonSerializer.Deserialize<ExplorerSession>(stream, Options)
                ?? throw new InvalidDataException("The session file is empty.");
            Validate(session);
            return session;
        }
        catch (JsonException error) { throw new InvalidDataException("This is not a valid Explorer session file.", error); }
    }

    public static string Save(string path, ExplorerSession session, bool createCopy = false)
    {
        Validate(session);
        var destination = Path.GetFullPath(path);
        var temporary = Path.Combine(Path.GetDirectoryName(destination)!, ".ec-session-" + Guid.NewGuid().ToString("N") + ".tmp");
        try
        {
            using (var stream = new FileStream(temporary, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            {
                JsonSerializer.Serialize(stream, session, Options);
                if (stream.Length > MaxFileBytes) throw new InvalidDataException(FileTooLargeMessage);
                stream.Flush(flushToDisk: true);
            }
            if (createCopy) return MoveToAvailablePath(temporary, destination);
            // Replace only after the complete snapshot has been written successfully.
            File.Move(temporary, destination, overwrite: true);
            return destination;
        }
        finally { if (File.Exists(temporary)) File.Delete(temporary); }
    }

    private static string MoveToAvailablePath(string temporary, string destination)
    {
        while (true)
        {
            destination = GetAvailablePath(destination);
            try
            {
                File.Move(temporary, destination, overwrite: false);
                return destination;
            }
            catch (IOException) when (File.Exists(destination) || Directory.Exists(destination))
            {
                // Another save took the suggested name; try the next one.
            }
        }
    }

    public static string GetAvailablePath(string path)
    {
        var destination = Path.GetFullPath(path);
        var directory = Path.GetDirectoryName(destination)!;
        var extension = Path.GetExtension(destination);
        var name = Path.GetFileNameWithoutExtension(destination);
        long nextNumber = 1;
        var suffixStart = name.LastIndexOf('(');
        if (suffixStart > 0 && name.EndsWith(')')
            && long.TryParse(name.AsSpan(suffixStart + 1, name.Length - suffixStart - 2),
                NumberStyles.None, CultureInfo.InvariantCulture, out var suffix)
            && suffix is > 0 and < long.MaxValue)
        {
            name = name[..suffixStart];
            nextNumber = suffix + 1;
        }

        var candidate = destination;
        while (true)
        {
            if (!File.Exists(candidate) && !Directory.Exists(candidate))
                return candidate;
            candidate = Path.Combine(directory, name + "(" + nextNumber.ToString(CultureInfo.InvariantCulture) + ")" + extension);
            nextNumber++;
        }
    }

    public static void Validate(ExplorerSession session)
    {
        if (session.Format != FormatName || session.Version != CurrentVersion)
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
        if (session.History == null || session.History.Count > ExplorerSession.HistoryLimit)
            throw new InvalidDataException("The session contains an invalid calculation history.");
        foreach (var job in session.History)
        {
            if (job?.Request is not { } request || request.OperationId == null
                || !CalculationCatalog.All.Any(operation => operation.Id == request.OperationId)
                || request.Equation == null || request.Equation.Length > 100_000
                || request.Arguments == null || request.Arguments.Count > 512
                || request.Arguments.Any(pair => pair.Key.Length > 1024 || pair.Value == null || pair.Value.Length > 100_000)
                || request.TimeoutSeconds is < 0 or > 86_400 || request.MaxItems is < 1 or > 100_000
                || !CalculationStatus.IsKnown(job.Status)
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

#nullable enable
using System.IO;
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
            var data = JsonSerializer.Deserialize<SessionFileData>(stream, Options)
                ?? throw new InvalidDataException("The session file is empty.");
            var session = data.ToSession();
            Validate(session);
            return session;
        }
        catch (JsonException error) { throw new InvalidDataException("This is not a valid Explorer session file.", error); }
    }

    public static string Save(string path, ExplorerSession session)
    {
        Validate(session);
        var destination = Path.GetFullPath(path);
        var temporary = Path.Combine(Path.GetDirectoryName(destination)!, ".ec-session-" + Guid.NewGuid().ToString("N") + ".tmp");
        try
        {
            using (var stream = new FileStream(temporary, FileMode.CreateNew, FileAccess.Write, FileShare.None))
            {
                JsonSerializer.Serialize(stream, SessionFileData.FromSession(session), Options);
                if (stream.Length > MaxFileBytes) throw new InvalidDataException(FileTooLargeMessage);
                stream.Flush(flushToDisk: true);
            }
            // Replace only after the complete snapshot has been written successfully.
            File.Move(temporary, destination, overwrite: true);
            return destination;
        }
        finally { if (File.Exists(temporary)) File.Delete(temporary); }
    }

    public static void Validate(ExplorerSession session)
    {
        if (session.Format != FormatName || session.Version != CurrentVersion)
            throw new InvalidDataException("This session format or version is not supported.");
        if (session.Equation == null || !CurveEquationText.TryParse(session.Equation, out _, out _))
            throw new InvalidDataException("The session contains an invalid curve equation.");
        if (session.History == null || session.History.Count > ExplorerSession.HistoryLimit)
            throw new InvalidDataException("The session contains an invalid calculation history.");
        foreach (var job in session.History)
        {
            if (job?.Request is not { } request || request.OperationId == null
                || !CalculationCatalog.All.Any(operation => operation.Id == request.OperationId)
                || request.Equation == null || request.Equation.Length > CurveEquationText.MaxTextLength
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

}

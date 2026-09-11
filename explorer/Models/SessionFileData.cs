#nullable enable
using System.Text.Json.Serialization;

namespace EllipticCurves.Explorer.Models;

// Only document data is persisted. Editor and presentation settings are local.
internal sealed record SessionFileData
{
    [JsonRequired] public string Format { get; init; } = SessionFile.FormatName;
    [JsonRequired] public int Version { get; init; } = SessionFile.CurrentVersion;
    public required string Equation { get; init; }
    public required List<CalculationSession> History { get; init; }

    internal static SessionFileData FromSession(ExplorerSession session) => new()
    {
        Format = session.Format, Version = session.Version, Equation = session.Equation, History = session.History
    };

    internal ExplorerSession ToSession() => (ExplorerSession.New() with
    {
        Format = Format, Version = Version, Equation = Equation, History = History
    }).DataOnly();
}

#nullable enable
namespace EllipticCurves.Explorer.Computations;

// These values are also stored in session files and displayed in saved reports.
public static class CalculationStatus
{
    public const string Running = "Running";
    public const string Completed = "Completed";
    public const string Cancelled = "Cancelled";
    public const string TimedOut = "Timed out";
    public const string Failed = "Failed";
    public const string Interrupted = "Interrupted";

    public static bool IsKnown(string? status) =>
        status is Running or Completed or Cancelled or TimedOut or Failed or Interrupted;
}

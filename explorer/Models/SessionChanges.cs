#nullable enable
namespace EllipticCurves.Explorer.Models;

// Visual settings are saved with the session, but only edits to its data need
// an unsaved-changes warning. Compare without serializing reports or committing edits.
public static class SessionChanges
{
    public static bool Equal(ExplorerSession left, ExplorerSession right) =>
        left.Equation == right.Equation && left.SliderStep == right.SliderStep && left.Preset == right.Preset
        && left.SliderOffsets.SequenceEqual(right.SliderOffsets)
        && left.History.Count == right.History.Count
        && left.History.Zip(right.History).All(pair => SameCalculation(pair.First, pair.Second));

    private static bool SameCalculation(CalculationSession left, CalculationSession right)
    {
        var a = left.Request;
        var b = right.Request;
        return left.StartedAt == right.StartedAt && left.Status == right.Status && left.Stage == right.Stage
            && left.Elapsed == right.Elapsed && left.Percent == right.Percent && left.Result == right.Result
            && a.OperationId == b.OperationId && a.Equation == b.Equation && a.TimeoutSeconds == b.TimeoutSeconds
            && a.MaxItems == b.MaxItems && a.Arguments.Count == b.Arguments.Count
            && a.Arguments.All(pair => b.Arguments.TryGetValue(pair.Key, out var value) && pair.Value == value);
    }
}

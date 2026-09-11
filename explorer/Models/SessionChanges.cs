#nullable enable
namespace EllipticCurves.Explorer.Models;

// Compare persisted state without serializing potentially large text reports or
// committing an incomplete edit just to decide whether it is safe to leave.
public static class SessionChanges
{
    public static bool Equal(ExplorerSession left, ExplorerSession right) =>
        left.Equation == right.Equation && left.SliderStep == right.SliderStep && left.Preset == right.Preset
        && left.SliderOffsets.SequenceEqual(right.SliderOffsets)
        && left.ShowGrid == right.ShowGrid && left.ShowPoints == right.ShowPoints && left.ComplexView == right.ComplexView
        && left.CoefficientsExpanded == right.CoefficientsExpanded
        && left.EquationScrollOffset == right.EquationScrollOffset && left.TorusScrollOffset == right.TorusScrollOffset
        // Background period mapping supplies O when no point has been selected.
        && PointKey(left.SelectedTorusPoint) == PointKey(right.SelectedTorusPoint)
        && left.Plot == right.Plot && left.TorusCamera == right.TorusCamera
        && left.EquationPanel == right.EquationPanel && left.ResultsPanel == right.ResultsPanel
        && left.SelectedResult == right.SelectedResult && left.History.Count == right.History.Count
        && left.History.Zip(right.History).All(pair => SameCalculation(pair.First, pair.Second));

    private static string PointKey(string? value) => value ?? EllipticCurvePoint.Infinity.ToString();

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

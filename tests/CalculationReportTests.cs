using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public class CalculationReportTests
{
    [Theory]
    [InlineData(nameof(EllipticCurveQ.GetRankLowerBound))]
    [InlineData(nameof(EllipticCurveQ.Regulator))]
    [InlineData(nameof(EllipticCurveQ.Saturate))]
    [InlineData(nameof(EllipticCurveQ.CreateIsogeny))]
    public void EquivalentPointInputsProduceTheSameReportAndRemainEditable(string method)
    {
        var operation = CalculationCatalog.All.First(o => o.Member?.Name == method);
        var parameter = operation.Parameters.Single(p => CalculationInput.IsPoints(p.ValueType));
        string expected = null;
        foreach (var input in new[] { "0.50, -2.500\nO", " (2/4, -10/4) \r\n o ", "0,5; -2,5\nO" })
        {
            var request = new CalculationRequest(operation.Id, "y^2 = x^3 - x",
                new() { [parameter.Key] = input }, MaxItems: 1);
            var job = new CalculationJobViewModel(request, operation.Title, new DateTime(2026, 9, 16));
            var report = job.Report;
            Assert.Contains(string.Join(Environment.NewLine, parameter.Label + ":",
                "  [0]: (1/2, -5/2)", "  [1]: O", "  Items displayed: 2 / 2"), report);
            expected ??= report;
            Assert.Equal(expected, report);
            Assert.Equal(input, job.Request.Arguments[parameter.Key]);

            job.Status = CalculationStatus.Completed;
            job.Result = "Completed result";
            var restored = CalculationJobViewModel.FromSession(job.CaptureSession());
            Assert.Equal(job.Report, restored.Report);
            using var workbench = new WorkbenchViewModel();
            using var repeat = new CalculationFormViewModel(operation, restored.Equation, workbench, restored.Request);
            Assert.Equal(input, repeat.Fields.Single(f => f.Parameter.Key == parameter.Key).Text);
        }
    }

    [Fact]
    public void EmptyAndInvalidPointInputsRemainReadable()
    {
        var operation = CalculationCatalog.All.Single(o => o.Member?.Name == nameof(EllipticCurveQ.GetRankLowerBound));
        var parameter = operation.Parameters.Single(p => p.Key == "points");
        Assert.Contains("Items displayed: 0 / 0", CalculationFormatter.FormatInput(parameter, " \n\t"));
        const string invalid = "(1, 2\ninvalid coordinate";
        var job = new CalculationJobViewModel(new(operation.Id, "y^2 = x^3 - x", new() { ["points"] = invalid }), operation.Title);
        job.Status = CalculationStatus.Failed;
        job.Result = "Invalid point input";
        Assert.Contains("points: " + invalid, job.Report);
        Assert.Contains("Invalid point input", job.Report);
    }
}

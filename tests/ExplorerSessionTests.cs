using System.Text.Json;
using System.Text.Json.Nodes;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerSessionTests
{
    [Fact]
    public async Task SavedTorusSelectionSurvivesDeferredMappingAndACachedView()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var points = new[] { new EllipticCurvePoint(-1, 0), new EllipticCurvePoint(0, 0) };
        using var model = new ComplexTorusViewModel();
        model.RestoreSelection(points[1].ToString());
        model.Update(curve, points, false);
        Assert.Equal(points[1].ToString(), model.SessionSelection);
        model.Update(curve, points, true);
        await model.PendingUpdate;
        Assert.Equal(points[1], model.SelectedPoint.Point);
        model.RestoreSelection(points[0].ToString());
        model.Update(curve, points, true);
        Assert.Equal(points[0], model.SelectedPoint.Point);
    }

    private static ExplorerSession Example() => new()
    {
        Equation = "y^2 + x*y + y = x^3 + 1/3*x^2 - 7/11*x + 2/5",
        SliderStep = "1/7", SliderOffsets = new[] { 2, -3, 4, -5, 6 }, ShowGrid = false, ShowPoints = false,
        ComplexView = true, CoefficientsExpanded = true, EquationScrollOffset = 37, TorusScrollOffset = 21,
        Plot = new(1.25, -2.5, 9.75), TorusCamera = new(-47, 26, 8.4), SelectedTorusPoint = "O",
        EquationPanel = new(false, 350), ResultsPanel = new(true, 420), SelectedResult = 1,
        History = new()
        {
            new(new("Q.TorsionStructure", "y^2 = x^3 - x", new() { ["example"] = "1/3\nπ" }, 35, 17),
                new DateTime(2026, 9, 11, 19, 20, 21, DateTimeKind.Local), "Completed", "Done", TimeSpan.FromSeconds(2.3), 100, "Exact result: ℚ\n1/3"),
            new(new("Q.TorsionStructure", "y^2 = x^3 + x", new()), DateTime.Now,
                "Running", "Computing", TimeSpan.FromSeconds(12), null, "Waiting for the calculation to finish…")
        }
    };

    [Fact]
    public void FileRoundTripRestoresExactCurveSlidersHistoryAndRequestWithoutRunningWork()
    {
        var path = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        try
        {
            var original = Example();
            SessionFile.Save(path, original);
            var loaded = SessionFile.Load(path);
            Assert.Equal(JsonSerializer.Serialize(original), JsonSerializer.Serialize(loaded));
            using var model = new MainViewModel();
            model.RestoreSession(loaded);
            Assert.Equal(original.Equation, model.Equation.Text);
            Assert.Equal(new BigRational(-7, 11), model.Snapshot.Curve.A4);
            Assert.Equal(new BigRational(1, 7), model.Step.ExactValue);
            Assert.Equal(original.SliderOffsets, model.ActiveCoefficients.Select(c => (int)c.SliderOffset));
            var coefficient = model.ActiveCoefficients[3];
            var minimum = coefficient.SliderMinimum;
            coefficient.SliderOffset++;
            Assert.Equal(new BigRational(-7, 11) + new BigRational(1, 7), coefficient.ExactValue);
            Assert.Equal(minimum, coefficient.SliderMinimum);
            using var workbench = new WorkbenchViewModel();
            workbench.RestoreHistory(loaded.History, loaded.SelectedResult);
            Assert.False(workbench.IsBusy);
            Assert.Null(workbench.Active);
            Assert.Same(workbench.Jobs[1], workbench.Selected);
            Assert.Equal("Interrupted", workbench.Selected.Status);
            Assert.Equal(original.History[0].StartedAt, workbench.Jobs[0].StartedAt);
            Assert.Equal(original.History[0].Result, workbench.Jobs[0].Result);
            Assert.Equal(35, workbench.Jobs[0].Request.TimeoutSeconds);
            Assert.Equal(17, workbench.Jobs[0].Request.MaxItems);
            Assert.Equal("1/3\nπ", workbench.Jobs[0].Request.Arguments["example"]);
            Assert.Contains("1/3", workbench.Jobs[0].Report);
        }
        finally { File.Delete(path); }
    }

    [Theory]
    [InlineData("invalid-json")]
    [InlineData("version")]
    [InlineData("missing-version")]
    [InlineData("format")]
    [InlineData("curve")]
    [InlineData("step")]
    [InlineData("plot")]
    [InlineData("panel")]
    [InlineData("history")]
    [InlineData("unknown-operation")]
    [InlineData("null-arguments")]
    [InlineData("selection")]
    public void InvalidFilesAreRejectedBeforeRestoringState(string defect)
    {
        var path = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        try
        {
            var json = JsonSerializer.SerializeToNode(Example());
            switch (defect)
            {
                case "version": json["Version"] = 999; break;
                case "missing-version": json.AsObject().Remove("Version"); break;
                case "format": json["Format"] = "Other app"; break;
                case "curve": json["Equation"] = "garbage"; break;
                case "step": json["SliderStep"] = "0"; break;
                case "plot": json["Plot"]["VerticalSpan"] = 0; break;
                case "panel": json["EquationPanel"]["Width"] = -1; break;
                case "history": json["History"] = null; break;
                case "unknown-operation": json["History"][0]["Request"]["OperationId"] = "missing"; break;
                case "null-arguments": json["History"][0]["Request"]["Arguments"] = null; break;
                case "selection": json["SelectedResult"] = 9; break;
            }
            File.WriteAllText(path, defect == "invalid-json" ? "{broken" : json.ToJsonString());
            Assert.Throws<InvalidDataException>(() => SessionFile.Load(path));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void InvalidSaveDoesNotReplaceExistingSession()
    {
        var path = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        try
        {
            var saved = Example();
            SessionFile.Save(path, saved);
            var original = File.ReadAllBytes(path);
            saved.History.Add(null);
            Assert.Throws<InvalidDataException>(() => SessionFile.Save(path, saved));
            Assert.Equal(original, File.ReadAllBytes(path));
        }
        finally { File.Delete(path); }
    }
}

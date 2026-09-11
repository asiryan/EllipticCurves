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
    public void ChangeComparisonIncludesEquationEditingParametersAndResults()
    {
        var baseline = ExplorerSession.New();
        SessionFile.Validate(baseline);
        Assert.True(SessionChanges.Equal(baseline, ExplorerSession.New()));
        var changes = new[]
        {
            baseline with { Equation = "y^2 = x^3 + x" },
            baseline with { SliderStep = "1/7" },
            baseline with { SliderOffsets = new[] { 1, 0 } },
            baseline with { Preset = null },
            baseline with { History = Example().History }
        };
        Assert.All(changes, changed => Assert.False(SessionChanges.Equal(baseline, changed)));
        var history = Example();
        Assert.False(SessionChanges.Equal(history, history with { History = history.History.Skip(1).ToList() }));
        Assert.False(SessionChanges.Equal(history, history with { History = new() }));
        var updated = history with { History = history.History.ToList() };
        Assert.True(SessionChanges.Equal(history, updated));
        updated.History[0] = updated.History[0] with { Result = "A changed result" };
        Assert.False(SessionChanges.Equal(history, updated));
        updated.History[0] = history.History[0] with
        {
            Request = history.History[0].Request with { Arguments = new() { ["example"] = "Different input" } }
        };
        Assert.False(SessionChanges.Equal(history, updated));
    }

    [Fact]
    public void VisualSettingsDoNotCountAsEditsOrHideDataChanges()
    {
        var baseline = ExplorerSession.New();
        var visualChanges = new[]
        {
            baseline with { ComplexView = true },
            baseline with { ShowGrid = false },
            baseline with { ShowPoints = false },
            baseline with { CoefficientsExpanded = true },
            baseline with { EquationPanel = new(false, 238) },
            baseline with { ResultsPanel = new(true, 400) },
            baseline with { EquationScrollOffset = 10 },
            baseline with { TorusScrollOffset = 20 },
            baseline with { SelectedTorusPoint = "(0, 0)" },
            baseline with { SelectedTorusPoint = null },
            baseline with { FitRealViewWhenShown = false },
            baseline with { Plot = new(200, -300, 70) },
            baseline with { TorusCamera = new(45, 30, 8) }
        };
        Assert.All(visualChanges, changed =>
        {
            Assert.True(SessionChanges.Equal(baseline, changed));
            Assert.True(SessionChanges.Equal(changed, baseline));
            Assert.False(SessionChanges.Equal(baseline, changed with { Equation = "y^2 = x^3 + 7" }));
            Assert.False(SessionChanges.Equal(baseline, changed with { History = Example().History }));
        });
    }

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

    [Fact]
    public async Task ChoosingOriginCancelsDeferredSessionSelection()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var savedPoint = new EllipticCurvePoint(0, 0);
        using var model = new ComplexTorusViewModel();
        model.RestoreSelection(savedPoint.ToString());
        model.Update(curve, Array.Empty<EllipticCurvePoint>(), true);
        await model.PendingUpdate;
        Assert.True(model.SelectedPoint.Point.IsInfinity);

        model.SelectedPoint = model.Points[0];
        Assert.Equal(EllipticCurvePoint.Infinity.ToString(), model.SessionSelection);
        model.Update(curve, new[] { savedPoint }, true);
        await model.PendingUpdate;
        Assert.True(model.SelectedPoint.Point.IsInfinity);
    }

    [Fact]
    public void RestoringSessionWithoutSliderOffsetsRecentersTheSliders()
    {
        using var model = new MainViewModel();
        model.ActiveCoefficients[0].SliderOffset = 7;
        model.FlushUpdate();
        var saved = ExplorerSession.New() with { Equation = model.Equation.Text, SliderOffsets = Array.Empty<int>() };
        SessionFile.Validate(saved);
        model.RestoreSession(saved);
        Assert.All(model.ActiveCoefficients, coefficient => Assert.Equal(0, coefficient.SliderOffset));
    }

    [Fact]
    public void SavingAnExistingLongFileNameDoesNotExceedTheFileSystemNameLimit()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-long-name-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var path = Path.Combine(directory, new string('s', 220) + ".ec");
        try
        {
            File.WriteAllText(path, JsonSerializer.Serialize(ExplorerSession.New()));
            var saved = ExplorerSession.New() with { Equation = "y^2 = x^3 + 7", Preset = null };
            SessionFile.Save(path, saved);
            Assert.Equal(saved.Equation, SessionFile.Load(path).Equation);
            Assert.Equal(new[] { path }, Directory.GetFiles(directory));
        }
        finally { File.Delete(path); Directory.Delete(directory); }
    }

    [Theory]
    [InlineData("session.ec", "session(1).ec")]
    [InlineData("session(1).ec", "session(2).ec")]
    [InlineData("curve.study.EC", "curve.study(1).EC")]
    [InlineData("curve (sample).ec", "curve (sample)(1).ec")]
    public void SavingACopyAddsANumberAndPreservesTheExistingFile(string name, string expected)
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-save-copy-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var original = Path.Combine(directory, name);
        try
        {
            Assert.Equal(original, SessionFile.GetAvailablePath(original));
            SessionFile.Save(original, ExplorerSession.New());
            var bytes = File.ReadAllBytes(original);
            var suggested = SessionFile.GetAvailablePath(original);
            Assert.Equal(Path.Combine(directory, expected), suggested);
            Assert.Equal(new[] { original }, Directory.GetFiles(directory));
            var edited = ExplorerSession.New() with { Equation = "y^2 = x^3 + 7", Preset = null };
            var actual = SessionFile.Save(original, edited, createCopy: true);
            Assert.Equal(suggested, actual);
            Assert.Equal(bytes, File.ReadAllBytes(original));
            Assert.Equal(edited.Equation, SessionFile.Load(actual).Equation);
            Assert.Equal(2, Directory.GetFiles(directory).Length);
        }
        finally
        {
            foreach (var file in Directory.GetFiles(directory)) File.Delete(file);
            Directory.Delete(directory);
        }
    }

    [Fact]
    public void SavingACopyRechecksTheNameSuggestedBeforeTheDialogOpened()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-save-suggestion-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var original = Path.Combine(directory, "session.ec");
        try
        {
            SessionFile.Save(original, ExplorerSession.New());
            var suggested = SessionFile.GetAvailablePath(original);
            Assert.Equal(Path.Combine(directory, "session(1).ec"), suggested);
            File.WriteAllText(suggested, "Created while the dialog was open");
            var actual = SessionFile.Save(suggested, ExplorerSession.New(), createCopy: true);
            Assert.Equal(Path.Combine(directory, "session(2).ec"), actual);
            Assert.Equal("Created while the dialog was open", File.ReadAllText(suggested));
            Assert.Equal("y^2 = x^3 - x", SessionFile.Load(actual).Equation);
        }
        finally
        {
            foreach (var file in Directory.GetFiles(directory)) File.Delete(file);
            Directory.Delete(directory);
        }
    }

    [Fact]
    public void SavingACopySkipsFilesAndDirectoriesButUsesAGapInTheNumbers()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-save-gap-" + Guid.NewGuid());
        var occupiedFolder = Path.Combine(directory, "session(2).ec");
        Directory.CreateDirectory(occupiedFolder);
        try
        {
            foreach (var name in new[] { "session.ec", "session(1).ec", "session(4).ec" })
                File.WriteAllText(Path.Combine(directory, name), "Keep this file");
            var actual = SessionFile.Save(Path.Combine(directory, "session.ec"), ExplorerSession.New(), createCopy: true);
            Assert.Equal(Path.Combine(directory, "session(3).ec"), actual);
            Assert.Equal("y^2 = x^3 - x", SessionFile.Load(actual).Equation);
            Assert.Equal("Keep this file", File.ReadAllText(Path.Combine(directory, "session(4).ec")));
        }
        finally
        {
            foreach (var file in Directory.GetFiles(directory)) File.Delete(file);
            Directory.Delete(occupiedFolder);
            Directory.Delete(directory);
        }
    }

    [Fact]
    public async Task ConcurrentCopiesKeepEverySnapshotAndDoNotOverwriteEachOther()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-save-concurrent-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var original = Path.Combine(directory, "session.ec");
        try
        {
            SessionFile.Save(original, ExplorerSession.New());
            var bytes = File.ReadAllBytes(original);
            var copies = await Task.WhenAll(Enumerable.Range(1, 8).Select(index => Task.Run(() =>
            {
                var equation = $"y^2 = x^3 + {index}";
                var path = SessionFile.Save(original, ExplorerSession.New() with { Equation = equation, Preset = null }, createCopy: true);
                return (path, equation);
            })));
            Assert.Equal(8, copies.Select(copy => copy.path).Distinct().Count());
            Assert.All(copies, copy => Assert.Equal(copy.equation, SessionFile.Load(copy.path).Equation));
            Assert.Equal(bytes, File.ReadAllBytes(original));
            Assert.Equal(9, Directory.GetFiles(directory).Length);
        }
        finally
        {
            foreach (var file in Directory.GetFiles(directory)) File.Delete(file);
            Directory.Delete(directory);
        }
    }

    private static ExplorerSession Example() => new()
    {
        Equation = "y^2 + x*y + y = x^3 + 1/3*x^2 - 7/11*x + 2/5",
        SliderStep = "1/7", SliderOffsets = new[] { 2, -3, 4, -5, 6 }, ShowGrid = false, ShowPoints = false,
        ComplexView = true, CoefficientsExpanded = true, EquationScrollOffset = 37, TorusScrollOffset = 21,
        Plot = new(1.25, -2.5, 9.75), TorusCamera = new(-47, 26, 8.4), SelectedTorusPoint = "O",
        EquationPanel = new(false, 350), ResultsPanel = new(true, 420),
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
            Assert.Equal("EllipticCurves.Explorer.Session", loaded.Format);
            Assert.Equal(1, loaded.Version);
            Assert.DoesNotContain("\"SelectedResult\"", File.ReadAllText(path));
            Assert.Equal(JsonSerializer.Serialize(original), JsonSerializer.Serialize(loaded));
            Assert.True(SessionChanges.Equal(original, loaded));
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
            workbench.RestoreHistory(loaded.History);
            Assert.False(workbench.IsBusy);
            Assert.Null(workbench.Active);
            Assert.Same(workbench.Jobs[0], workbench.Selected);
            Assert.Equal("Interrupted", workbench.Jobs[1].Status);
            Assert.Equal(original.History[0].StartedAt, workbench.Jobs[0].StartedAt);
            Assert.Equal(original.History[0].Result, workbench.Jobs[0].Result);
            Assert.Equal(35, workbench.Jobs[0].Request.TimeoutSeconds);
            Assert.Equal(17, workbench.Jobs[0].Request.MaxItems);
            Assert.Equal("1/3\nπ", workbench.Jobs[0].Request.Arguments["example"]);
            Assert.Contains("1/3", workbench.Jobs[0].Report);
            workbench.RestoreHistory(Array.Empty<CalculationSession>());
            Assert.Empty(workbench.Jobs);
            Assert.Null(workbench.Selected);
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

    [Fact]
    public void LockedDestinationPreservesTheFileAndCleansTheTemporarySnapshot()
    {
        if (!OperatingSystem.IsWindows()) return;
        var directory = Path.Combine(Path.GetTempPath(), "ec-locked-save-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var path = Path.Combine(directory, "session.ec");
        try
        {
            SessionFile.Save(path, ExplorerSession.New());
            var original = File.ReadAllBytes(path);
            using (var locked = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                var error = Record.Exception(() => SessionFile.Save(path, Example()));
                Assert.True(error is IOException or UnauthorizedAccessException, "A locked session file must reject replacement.");
            }
            Assert.Equal(original, File.ReadAllBytes(path));
            Assert.Equal(new[] { path }, Directory.GetFiles(directory));
        }
        finally { File.Delete(path); Directory.Delete(directory); }
    }
}

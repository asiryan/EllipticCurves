using System.Numerics;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerEquationTests
{
    [Theory]
    [InlineData("y^2 = x^3 - 106.16*x - 0.32", "0", "0", "0", "-2654/25", "-8/25")]
    [InlineData("y^2 = x^3 − 106,16x − 0,32", "0", "0", "0", "-2654/25", "-8/25")]
    [InlineData("y^2 + xy + y = x^3 - x", "1", "0", "1", "-1", "0")]
    [InlineData("y^2+8.325*x*y-2/7*y=x^3+3/5*x^2-1e-5*x+1000", "333/40", "3/5", "-2/7", "-1/100000", "1000")]
    [InlineData("x^3 - x = y^2 + xy + y", "1", "0", "1", "-1", "0")]
    [InlineData("y^2 + xy + y - x^3 + x = 0", "1", "0", "1", "-1", "0")]
    [InlineData("0 = -y^2 + x^3 - x", "0", "0", "0", "-1", "0")]
    [InlineData("-2y^2 = -2x^3 + 4x - 2", "0", "0", "0", "-2", "1")]
    [InlineData("y*(y+1) = x*(x*x-1)", "0", "0", "1", "-1", "0")]
    [InlineData("y(y+1) = x(x^2-1)", "0", "0", "1", "-1", "0")]
    [InlineData("y^2 = x^3 + (1/2)x + 1/(2+3)", "0", "0", "0", "1/2", "1/5")]
    [InlineData("y^2 = x^3 - x/7 - .00025", "0", "0", "0", "-1/7", "-1/4000")]
    [InlineData("Y^2 = X^3 + 1E+3·X - 8.325", "0", "0", "0", "1000", "-333/40")]
    [InlineData("y^2 = x^3 + 2x - 3x + x + 1 - 1", "0", "0", "0", "0", "0")]
    [InlineData("(y+1)^2 = x^3", "0", "0", "2", "0", "-1")]
    public void ParsesExactWeierstrassEquations(string input, string a1, string a2, string a3, string a4, string a6)
    {
        Assert.True(CurveEquationText.TryParse(input, out var curve, out var error), error);
        var expected = new[] { a1, a2, a3, a4, a6 }.Select(ParseRational).ToArray();
        Assert.Equal(expected, new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 });
    }

    [Theory]
    [InlineData("")]
    [InlineData("y² = x³ - x")]
    [InlineData("y^2 = x³ - x")]
    [InlineData("y^2 =")]
    [InlineData("y^2 = x^3 -")]
    [InlineData("y^2 = x^3 + 1/")]
    [InlineData("y^2 = x^3 + 1e-")]
    [InlineData("y^2 = x^3 + 1/0")]
    [InlineData("y^2 = x^3 + 1/(1-1)")]
    [InlineData("y^2 = x^3 + 1/x")]
    [InlineData("y^2 = x^3 + z")]
    [InlineData("y^2 = x^3 + sin(x)")]
    [InlineData("y^2 = x^3 = 0")]
    [InlineData("y^2 = x^3 + (x")]
    [InlineData("y^2 = x^4")]
    [InlineData("y^2 = x^3 + x^2y")]
    [InlineData("y^2 = x^30")]
    [InlineData("y^2 = 2x^3 - x")]
    [InlineData("y^2 = x^2")]
    [InlineData("y = x^3")]
    [InlineData("y^2 = x^3 + 1 2")]
    [InlineData("y^2 = x^3 + 1e5000")]
    public void RejectsUnsupportedOrIncompleteExpressions(string input)
    {
        Assert.False(CurveEquationText.TryParse(input, out var curve, out var error));
        Assert.Null(curve);
        Assert.NotEmpty(error);
    }

    [Fact]
    public void GuardsInputLengthAndNesting()
    {
        Assert.False(CurveEquationText.TryParse("y^2=x^3+" + new string('1', CurveEquationText.MaxTextLength), out _, out _));
        Assert.False(CurveEquationText.TryParse("y^2=x^3+" + new string('(', 100) + "1" + new string(')', 100), out _, out _));
        Assert.False(CurveEquationText.TryParse("y^2=x^3" + string.Concat(Enumerable.Repeat("+1", 4096)), out _, out var error));
        Assert.Contains("too many terms", error);
    }

    [Fact]
    public void MaximumSizeCoefficientsRoundTripInAllFivePositions()
    {
        var large = BigInteger.One << (RationalText.MaxValueBits - 1);
        var fraction = new BigRational(large - 1, large - 3);
        var original = new EllipticCurveQ(fraction, -fraction, fraction, -fraction, fraction);
        var formatted = CurveEquationText.Format(original);
        Assert.InRange(formatted.Length, 4097, CurveEquationText.MaxTextLength);
        Assert.True(CurveEquationText.TryParse(formatted, out var parsed, out var error), error);
        Assert.Equal(original, parsed);

        var tiny = new BigRational(1, large);
        formatted = RationalText.Format(tiny);
        Assert.InRange(formatted.Length, 1, RationalText.MaxTextLength);
        Assert.True(RationalText.TryParse(formatted, out var restored));
        Assert.Equal(tiny, restored);
    }

    [Theory]
    [InlineData("1e4096")]
    [InlineData("1e-4096")]
    public async Task ScientificCoefficientsSurviveSlidersCalculationsAndSessionFiles(string coefficient)
    {
        using var model = new MainViewModel();
        model.ShowPoints = false;
        model.Equation.Text = "y^2=x^3+" + coefficient + "*x";
        model.FlushUpdate();
        var original = model.Snapshot.Curve;
        var request = new CalculationRequest("curve.overview", CurveEquationText.Format(original), new());
        Assert.Contains("Discriminant", await CalculationEngine.ExecuteAsync(request));
        model.SimpleCoefficients[0].SliderOffset = 1;
        model.FlushUpdate();
        Assert.True(model.Equation.IsValid);
        Assert.Equal(original.A4 + new BigRational(1, 100), model.Snapshot.Curve.A4);
        Assert.True(CurveEquationText.TryParse(model.Equation.Text, out var parsed, out var error), error);
        Assert.Equal(model.Snapshot.Curve, parsed);
        var path = Path.Combine(Path.GetTempPath(), "ec-scientific-" + Guid.NewGuid() + ".ec");
        try
        {
            SessionFile.Save(path, ExplorerSession.New() with { Equation = model.Equation.Text });
            var saved = SessionFile.Load(path);
            using var reopened = new MainViewModel();
            reopened.ShowPoints = false;
            reopened.RestoreSession(saved);
            Assert.Equal(model.Snapshot.Curve, reopened.Snapshot.Curve);
            Assert.Contains("Discriminant", await CalculationEngine.ExecuteAsync(request with { Equation = saved.Equation }));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    public void SliderOverflowIsInvalidAndCanBeUndoneByReturningToItsAnchor()
    {
        using var model = new MainViewModel();
        model.ShowPoints = false;
        var largest = new BigRational((BigInteger.One << RationalText.MaxValueBits) - 1);
        model.Equation.SetCurve(new EllipticCurveQ(0, 0, 0, largest, 0));
        model.FlushUpdate();
        var snapshot = model.Snapshot;
        model.Step.Text = "1";
        model.SimpleCoefficients[0].SliderOffset = 1;
        model.FlushUpdate();
        Assert.False(model.Equation.IsValid);
        Assert.True(model.HasInputError);
        Assert.Same(snapshot, model.Snapshot);
        Assert.Throws<InvalidDataException>(() => SessionFile.Validate(ExplorerSession.New() with { Equation = model.Equation.Text }));
        model.SimpleCoefficients[0].SliderOffset = 0;
        model.FlushUpdate();
        Assert.True(model.Equation.IsValid);
        Assert.Equal(largest, model.Snapshot.Curve.A4);
    }

    [Theory]
    [InlineData("2")]
    [InlineData("1/2")]
    public void RejectsExplosiveNestedPowers(string constant)
    {
        // Cross the parser's limit without exhausting memory if the guard regresses.
        for (var i = 0; i < 10; i++) constant = "(" + constant + ")^3";
        Assert.False(CurveEquationText.TryParse("y^2=x^3+" + constant, out var curve, out var error));
        Assert.Null(curve);
        Assert.Contains("too large", error);
    }

    [Fact]
    public void GuardsGrowthDuringLeadingCoefficientNormalization()
    {
        Assert.False(CurveEquationText.TryParse("(1e-4000)^2*y^2=(1e-4000)^2*x^3+1e4000", out _, out var error));
        Assert.Contains("too large", error);
    }

    [Fact]
    public void LargeScientificCoefficientsRemainExact()
    {
        Assert.True(CurveEquationText.TryParse("y^2=x^3+1e4000*x+1e-4000", out var curve, out var error), error);
        Assert.Equal(ParseRational("1e4000"), curve.A4);
        Assert.Equal(ParseRational("1e-4000"), curve.A6);
    }

    [Fact]
    public void FormattedEquationsRoundTripExactly()
    {
        var random = new Random(8427);
        for (var i = 0; i < 60; i++)
        {
            var values = Enumerable.Range(0, 5).Select(_ => new BigRational(random.Next(-10000, 10000), random.Next(1, 80))).ToArray();
            var original = new EllipticCurveQ(values[0], values[1], values[2], values[3], values[4]);
            Assert.True(CurveEquationText.TryParse(CurveEquationText.Format(original), out var parsed, out var error), error);
            Assert.Equal(values, new[] { parsed.A1, parsed.A2, parsed.A3, parsed.A4, parsed.A6 });
            Assert.Equal(original.Discriminant, parsed.Discriminant);
        }
    }

    [Fact]
    public async Task FormulaEditsDebounceWithoutRewritingTheTextOrResettingTheViewport()
    {
        using var model = new MainViewModel();
        model.ShowPoints = false;
        var previous = model.Snapshot;
        var resets = 0;
        model.ViewResetRequested += (_, _) => resets++;
        model.CurveResetRequested += (_, _) => resets++;
        const string input = "y^2 = x^3 - 106.16x - 0.32";
        model.Equation.Text = input;
        Assert.Same(previous, model.Snapshot);
        await model.PendingUpdate;
        Assert.Equal(ParseRational("-2654/25"), model.Snapshot.Curve.A4);
        Assert.Equal(ParseRational("-8/25"), model.Snapshot.Curve.A6);
        Assert.Equal(input, model.Equation.Text);
        Assert.Equal(0, resets);
    }

    [Fact]
    public async Task IncompleteFormulaCancelsPendingCurveAndOnlyValidatesOnCommit()
    {
        using var model = new MainViewModel();
        var previous = model.Snapshot;
        model.Equation.Text = "y^2=x^3+8.325x";
        var pending = model.PendingUpdate;
        model.Equation.Text = "y^2=x^3+1/";
        Assert.True(model.HasIncompleteInput);
        Assert.False(model.HasInputError);
        Assert.Empty(model.Equation.Error);
        await pending;
        Assert.Same(previous, model.Snapshot);
        model.Equation.CommitEdit();
        Assert.True(model.HasInputError);
        model.Equation.Text = "y^2=x^3+1/7";
        model.FlushUpdate();
        Assert.False(model.HasInputError);
        Assert.Equal(new BigRational(1, 7), model.Snapshot.Curve.A6);
    }

    [Fact]
    public void SlidersUpdateTheFormulaAndResetCanRecoverFromInvalidInput()
    {
        using var model = new MainViewModel();
        model.Equation.Text = "y^2=x^3+8.325x";
        model.SimpleCoefficients[0].SliderOffset = 1;
        model.FlushUpdate();
        Assert.Equal("y^2 = x^3 + 8.335*x", model.Equation.Text);
        Assert.Equal(ParseRational("8.335"), model.Snapshot.Curve.A4);
        model.Equation.Text = "y^2=x^3-";
        model.Equation.CommitEdit();
        model.SimpleCoefficients[0].ResetCommand.Execute(null);
        model.FlushUpdate();
        Assert.Equal("y^2 = x^3", model.Equation.Text);
        Assert.True(model.Equation.IsValid);
        Assert.True(model.Snapshot.IsSingular);
    }

    [Fact]
    public void PresetCancelsInvalidFormulaAndSynchronizesAllCoefficients()
    {
        using var model = new MainViewModel();
        model.Equation.Text = "y^2=x^3-";
        model.ApplyPreset(CurvePreset.All.Single(preset => preset.Name == "37.a1"));
        Assert.Equal("y^2 + y = x^3 - x", model.Equation.Text);
        Assert.True(model.IsGeneralForm);
        Assert.Equal(BigRational.One, model.Coefficients[2].ExactValue);
        Assert.False(model.HasIncompleteInput);
        Assert.Equal(model.Snapshot.Curve.A3, model.Equation.Curve.A3);
    }

    private static BigRational ParseRational(string text)
    {
        Assert.True(RationalText.TryParse(text, out var value));
        return value;
    }
}

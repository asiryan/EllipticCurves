using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerModelTests
{
    [Fact]
    public void SnapshotUsesExactLibraryInvariants()
    {
        var snapshot = new CurveSnapshot(new EllipticCurveQ(0, 0, 0, -1, 0));
        Assert.Equal("64", snapshot.Discriminant);
        Assert.Equal("1728", snapshot.JInvariant);
        Assert.Equal("2", snapshot.Components);
        Assert.Equal("48", snapshot.C4);
        Assert.Equal("0", snapshot.C6);
        Assert.False(snapshot.IsSingular);
    }

    [Fact]
    public void SingularSnapshotDoesNotRequestUndefinedInvariants()
    {
        var snapshot = new CurveSnapshot(new EllipticCurveQ(0, 0, 0, 0, 0));
        Assert.True(snapshot.IsSingular);
        Assert.Equal("Undefined", snapshot.JInvariant);
        Assert.Equal("—", snapshot.Components);
        Assert.Equal("0", snapshot.Discriminant);
    }

    [Theory]
    [InlineData("1.2", 6, 5)]
    [InlineData("-0.1", -1, 10)]
    [InlineData("3.5", 7, 2)]
    [InlineData("100", 100, 1)]
    [InlineData("-100", -100, 1)]
    [InlineData("8.325", 333, 40)]
    [InlineData("333/40", 333, 40)]
    [InlineData(" -2 / 7 ", -2, 7)]
    [InlineData("1e-5", 1, 100000)]
    [InlineData("-1.25E+3", -1250, 1)]
    [InlineData("100.1", 1001, 10)]
    [InlineData("-1000000", -1000000, 1)]
    [InlineData(".00025", 1, 4000)]
    [InlineData("−2.5", -5, 2)]
    public void DecimalInputRemainsExact(string input, int numerator, int denominator)
    {
        var coefficient = new CoefficientViewModel("a₄", "x", () => { }) { Text = input };
        Assert.True(coefficient.IsValid);
        Assert.Equal(new BigRational(numerator, denominator), coefficient.ExactValue);
    }

    [Theory]
    [InlineData("")]
    [InlineData("-")]
    [InlineData("NaN")]
    [InlineData("1/")]
    [InlineData("1/0")]
    [InlineData("1/2/3")]
    [InlineData("1e-")]
    [InlineData("1.2.3")]
    [InlineData("3,5")]
    [InlineData("8,325")]
    [InlineData("−2,5")]
    [InlineData("Infinity")]
    [InlineData("1e999999999")]
    public void InvalidInputPreservesLastCurveAndCanBeCorrected(string input)
    {
        using var model = new MainViewModel();
        var previous = model.Snapshot;
        var coefficient = model.SimpleCoefficients[0];
        coefficient.Text = input;
        Assert.True(model.HasIncompleteInput);
        Assert.False(model.HasInputError);
        Assert.Empty(coefficient.Error);
        coefficient.CommitEdit();
        Assert.True(model.HasInputError);
        Assert.Same(previous, model.Snapshot);
        coefficient.Text = "-2";
        model.FlushUpdate();
        Assert.False(model.HasInputError);
        Assert.Equal(new BigRational(-2), model.Snapshot.Curve.A4);
    }

    [Fact]
    public async Task UnchangedInputKeepsThePresetSnapshotAndSamples()
    {
        using var model = new MainViewModel();
        await model.PendingSamples;
        var preset = model.SelectedPreset;
        var snapshot = model.Snapshot;
        var samples = model.Samples;
        Assert.NotEmpty(samples);

        // Opening a calculation commits the equation and flushes pending input.
        model.Equation.CommitEdit();
        model.FlushUpdate();
        model.SimpleCoefficients[0].CommitEdit();
        model.Equation.Text = "y^2=x^3-1.0*x";
        model.Equation.CommitEdit();
        model.ShowPoints = true;
        await model.PendingUpdate;
        model.FlushUpdate();

        Assert.Same(preset, model.SelectedPreset);
        Assert.Same(snapshot, model.Snapshot);
        Assert.Same(samples, model.Samples);
        Assert.Equal("y^2=x^3-1.0*x", model.Equation.Text);
        Assert.DoesNotContain("Updating", model.InputStatus);
    }

    [Fact]
    public async Task RevertingAnEditCancelsThePendingReplacement()
    {
        using var model = new MainViewModel();
        await model.PendingSamples;
        var snapshot = model.Snapshot;
        var samples = model.Samples;
        model.Equation.Text = "y^2=x^3+2*x";
        var pending = model.PendingUpdate;
        model.Equation.Text = "y^2=x^3-x";
        await pending;
        model.FlushUpdate();
        Assert.Same(snapshot, model.Snapshot);
        Assert.Same(samples, model.Samples);
        Assert.DoesNotContain("Updating", model.InputStatus);
    }

    [Fact]
    public void UnchangedStepDoesNotRecenterTheSliders()
    {
        using var model = new MainViewModel();
        var coefficient = model.SimpleCoefficients[0];
        coefficient.SliderOffset = 20;
        var value = coefficient.ExactValue;
        var minimum = coefficient.SliderMinimum;
        var maximum = coefficient.SliderMaximum;
        model.Step.CommitEdit();
        model.Step.Text = "1/100";
        model.Step.CommitEdit();
        model.SetStepCommand.Execute("0.01");
        Assert.Equal(20, coefficient.SliderOffset);
        Assert.Equal(value, coefficient.ExactValue);
        Assert.Equal(minimum, coefficient.SliderMinimum);
        Assert.Equal(maximum, coefficient.SliderMaximum);

        model.Step.Text = "0.1";
        Assert.Equal(0, coefficient.SliderOffset);
        coefficient.SliderOffset = 1;
        Assert.Equal(value + new BigRational(1, 10), coefficient.ExactValue);
    }

    [Fact]
    public void StepsAndSlidersPreserveTheExactAnchorWithoutClampingCoefficients()
    {
        using var model = new MainViewModel();
        var coefficient = model.SimpleCoefficients[0];
        coefficient.Text = "8.325";
        coefficient.SliderOffset = 1;
        Assert.Equal("8.335", coefficient.Text);
        coefficient.SliderOffset = 0;
        Assert.Equal(new BigRational(333, 40), coefficient.ExactValue);
        coefficient.Text = "1000.325";
        coefficient.SliderOffset = 2;
        Assert.Equal("1000.345", coefficient.Text);
        coefficient.SliderOffset = 0;
        Assert.Equal("1000.325", coefficient.Text);
        model.Step.Text = "1/7";
        Assert.Equal(0, coefficient.SliderOffset);
        coefficient.SliderOffset = 1;
        Assert.Equal(new BigRational(40013, 40) + new BigRational(1, 7), coefficient.ExactValue);
        coefficient.ResetCommand.Execute(null);
        Assert.Equal(BigRational.Zero, coefficient.ExactValue);
        Assert.Equal("0", coefficient.Text);
        Assert.Equal(0, coefficient.SliderOffset);
    }

    [Theory]
    [InlineData("0")]
    [InlineData("-1")]
    [InlineData("1/")]
    public void InvalidStepDoesNotChangeCoefficientAndDoesNotPreventReset(string input)
    {
        using var model = new MainViewModel();
        var coefficient = model.SimpleCoefficients[0];
        coefficient.Text = "8.325";
        model.Step.Text = input;
        coefficient.SliderOffset = 3;
        Assert.Equal(new BigRational(333, 40), coefficient.ExactValue);
        coefficient.Text = "-";
        coefficient.CommitEdit();
        coefficient.ResetCommand.Execute(null);
        Assert.True(coefficient.IsValid);
        Assert.Empty(coefficient.Error);
        Assert.Equal(BigRational.Zero, coefficient.ExactValue);
    }

    [Fact]
    public void LongDecimalsAndLargeRatiosDoNotPassThroughDecimalOrDouble()
    {
        const string input = "0.123456789012345678901234567890123456789";
        Assert.True(RationalText.TryParse(input, out var value));
        Assert.Equal(new BigRational(System.Numerics.BigInteger.Parse("123456789012345678901234567890123456789"), System.Numerics.BigInteger.Pow(10, 39)), value);
        Assert.Equal(input, RationalText.Format(value));
        var large = System.Numerics.BigInteger.Pow(10, 400);
        Assert.Equal(1.0, CurvePlotData.ToDouble(new BigRational(large + 1, large - 1)), 12);
    }

    [Fact]
    public void PresetsResetBothViewsWhileFitTargetsOnlyTheCurrentView()
    {
        using var model = new MainViewModel();
        model.ShowPoints = false;
        var curveResets = 0;
        var viewResets = 0;
        model.CurveResetRequested += (_, _) => curveResets++;
        model.ViewResetRequested += (_, _) => viewResets++;
        model.ApplyPreset(CurvePreset.All.Single(preset => preset.Name == "48.a3"));
        Assert.Equal((1, 0), (curveResets, viewResets));
        model.FitCommand.Execute(null);
        Assert.Equal((1, 1), (curveResets, viewResets));
        model.ResetCommand.Execute(null);
        Assert.Equal((2, 1), (curveResets, viewResets));
    }

    [Fact]
    public async Task EditingIsDebouncedAndNeverResetsTheView()
    {
        using var model = new MainViewModel();
        model.ShowPoints = false;
        var previous = model.Snapshot;
        var updates = 0;
        var viewResets = 0;
        model.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(model.Snapshot)) updates++; };
        model.ViewResetRequested += (_, _) => viewResets++;
        model.CurveResetRequested += (_, _) => viewResets++;
        model.SimpleCoefficients[0].Text = "8.3";
        var superseded = model.PendingUpdate;
        model.SimpleCoefficients[0].Text = "8.325";
        Assert.Same(previous, model.Snapshot);
        await Task.WhenAll(superseded, model.PendingUpdate);
        Assert.Equal(new BigRational(333, 40), model.Snapshot.Curve.A4);
        Assert.Equal(1, updates);
        Assert.Equal(0, viewResets);
        var plotted = model.Snapshot;
        model.SimpleCoefficients[0].Text = "9";
        var abandoned = model.PendingUpdate;
        model.SimpleCoefficients[0].Text = "1/";
        await abandoned;
        Assert.Same(plotted, model.Snapshot);
    }

    [Fact]
    public void FormulaAutomaticallySelectsCoefficientsForBothForms()
    {
        using var model = new MainViewModel();
        model.Equation.Text = "y^2 + 2/7*xy = x^3 + 1000x";
        model.FlushUpdate();
        Assert.False(model.IsSimpleForm);
        Assert.Equal(5, model.ActiveCoefficients.Count);
        Assert.Equal(new BigRational(2, 7), model.Snapshot.Curve.A1);
        Assert.Equal(new BigRational(1000), model.Snapshot.Curve.A4);
        model.Equation.Text = "y^2 = x^3 + 8.325x";
        model.FlushUpdate();
        Assert.True(model.IsSimpleForm);
        Assert.Equal(2, model.ActiveCoefficients.Count);
        Assert.False(model.HasIncompleteInput);
        Assert.Equal(BigRational.Zero, model.Snapshot.Curve.A1);
        Assert.Equal(new BigRational(333, 40), model.Snapshot.Curve.A4);
    }

    [Fact]
    public void UnrepresentablePlotRetainsExactInvariantsWithoutCrashing()
    {
        Assert.True(RationalText.TryParse("1e400", out var coefficient));
        var snapshot = new CurveSnapshot(new EllipticCurveQ(0, 0, 0, coefficient, 1));
        Assert.Equal(coefficient, snapshot.Curve.A4);
        Assert.True(snapshot.IsPlotUnavailable);
        Assert.False(snapshot.Plot.TryEvaluate(0, out _, out _));
        Assert.NotEmpty(snapshot.Discriminant);
    }

    [Fact]
    public void GeneralModelUsesBothLinearYTerms()
    {
        var data = new CurvePlotData(new EllipticCurveQ(2, 0, 3, -1, 1));
        Assert.True(data.TryEvaluate(1, out var upper, out var lower));
        Assert.Equal((-5 + Math.Sqrt(29)) / 2, upper, 12);
        Assert.Equal((-5 - Math.Sqrt(29)) / 2, lower, 12);
    }

    [Fact]
    public void RealRootsSeparateDisconnectedBranches()
    {
        var data = new CurvePlotData(new EllipticCurveQ(0, 0, 0, -1, 0));
        Assert.Equal(3, data.Roots.Count);
        Assert.Equal(-1, data.Roots[0], 10);
        Assert.Equal(0, data.Roots[1], 10);
        Assert.Equal(1, data.Roots[2], 10);
        Assert.False(data.TryEvaluate(-2, out _, out _));
        Assert.True(data.TryEvaluate(-0.5, out _, out _));
        Assert.False(data.TryEvaluate(0.5, out _, out _));
        Assert.True(data.TryEvaluate(2, out _, out _));
    }

    [Theory]
    [InlineData("1e-12", 1e-6)]
    [InlineData("1e-100", 1e-50)]
    [InlineData("1e-300", 1e-150)]
    public void SmallCurvesPreserveForbiddenIntervalsWithoutUnderflow(string coefficient, double scale)
    {
        Assert.True(RationalText.TryParse(coefficient, out var value));
        var data = new CurvePlotData(new EllipticCurveQ(0, 0, 0, -value, 0));
        Assert.True(data.IsDrawable);
        Assert.Equal(3, data.Roots.Count);
        Assert.False(data.TryEvaluate(-2 * scale, out _, out _));
        Assert.False(data.TryEvaluate(0.5 * scale, out _, out _));
        Assert.True(data.TryEvaluate(-0.5 * scale, out var upper, out var lower));
        Assert.True(upper > 0 && double.IsFinite(upper));
        Assert.Equal(-upper, lower);
        var expected = Math.Sqrt(0.375) * Math.Sqrt(scale) * scale;
        Assert.InRange(upper / expected, 1 - 1e-12, 1 + 1e-12);
        Assert.True(data.TryEvaluate(0, out upper, out lower));
        Assert.Equal(0, upper);
        Assert.Equal(0, lower);
    }

    [Fact]
    public void SingularIsolatedPointIsNotLostBySignChangeRootSearch()
    {
        var data = new CurvePlotData(new EllipticCurveQ(0, 0, 0, -3, -2));
        Assert.Equal(2, data.Roots.Count);
        Assert.Equal(-1, data.Roots[0], 10);
        Assert.Equal(2, data.Roots[1], 10);
        Assert.True(data.TryEvaluate(-1, out var upper, out var lower));
        Assert.Equal(0, upper, 12);
        Assert.Equal(0, lower, 12);
        Assert.False(data.TryEvaluate(-0.9, out _, out _));
        Assert.False(data.TryEvaluate(-1.1, out _, out _));
    }

    [Fact]
    public void NumericalBranchesSatisfyOriginalEquationForGeneralCoefficients()
    {
        var random = new Random(3701);
        for (var example = 0; example < 100; example++)
        {
            var a = Enumerable.Range(0, 5).Select(_ => new BigRational(random.Next(-1000, 1001), 10)).ToArray();
            var data = new CurvePlotData(new EllipticCurveQ(a[0], a[1], a[2], a[3], a[4]));
            var coefficients = a.Select(CurvePlotData.ToDouble).ToArray();
            for (var x = -20.0; x <= 20; x += 0.25)
            {
                if (!data.TryEvaluate(x, out var upper, out var lower)) continue;
                foreach (var y in new[] { upper, lower })
                {
                    var lhs = y * y + coefficients[0] * x * y + coefficients[2] * y;
                    var rhs = x * x * x + coefficients[1] * x * x + coefficients[3] * x + coefficients[4];
                    var scale = 1 + Math.Abs(y * y) + Math.Abs(coefficients[0] * x * y) + Math.Abs(coefficients[2] * y) + Math.Abs(rhs);
                    Assert.True(Math.Abs(lhs - rhs) <= 1e-10 * scale);
                }
            }
        }
    }

    [Fact]
    public async Task RapidEditsOnlyPublishSamplesForLatestCurve()
    {
        using var model = new MainViewModel();
        for (var i = 0; i < 12; i++) model.SimpleCoefficients[1].SetExact(i);
        model.ApplyPreset(CurvePreset.All.Single(preset => preset.Name == "37.a1"));
        await model.PendingSamples;
        Assert.NotEmpty(model.Samples);
        Assert.All(model.Samples, point => Assert.True(model.Snapshot.Curve.IsOnCurve(point)));
        Assert.Contains("|m| ≤ 12", model.SampleStatus);
        model.ApplyPreset(CurvePreset.All.Single(preset => preset.Name == "The cusp"));
        await model.PendingSamples;
        Assert.Empty(model.Samples);
        Assert.True(model.Snapshot.IsSingular);
    }

    [Fact]
    public async Task HidingSamplesCancelsPendingPublication()
    {
        using var model = new MainViewModel();
        var pending = model.PendingSamples;
        model.ShowPoints = false;
        await pending;
        Assert.Empty(model.Samples);
        Assert.Contains("hidden", model.SampleStatus);
    }

    [Fact]
    public async Task ClosingDuringSearchIsSafeAndDisposalIsIdempotent()
    {
        var model = new MainViewModel();
        var pending = model.PendingSamples;
        var snapshot = model.Snapshot;
        model.SimpleCoefficients[0].Text = "8.325";
        var pendingUpdate = model.PendingUpdate;
        model.Dispose();
        model.Dispose();
        await Task.WhenAll(pending, pendingUpdate);
        Assert.Same(snapshot, model.Snapshot);
        Assert.Empty(model.Samples);
    }
}

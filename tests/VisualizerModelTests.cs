using EllipticCurves.Visualizer.Models;
using EllipticCurves.Visualizer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class VisualizerModelTests
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
    [InlineData("3,5", 7, 2)]
    [InlineData("100", 100, 1)]
    [InlineData("-100", -100, 1)]
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
    [InlineData("1.25")]
    [InlineData("100.1")]
    [InlineData("-101")]
    public void InvalidInputPreservesLastCurveAndCanBeCorrected(string input)
    {
        using var model = new MainViewModel();
        var previous = model.Snapshot;
        model.Coefficients[3].Text = input;
        Assert.True(model.HasInputError);
        Assert.Same(previous, model.Snapshot);
        model.Coefficients[3].Value = -2;
        Assert.False(model.HasInputError);
        Assert.Equal(new BigRational(-2), model.Snapshot.Curve.A4);
    }

    [Fact]
    public void SliderUsesTenthsAndClampsAtItsBounds()
    {
        var coefficient = new CoefficientViewModel("a₄", "x", () => { }) { Value = 0.30000000000000004 };
        Assert.Equal(new BigRational(3, 10), coefficient.ExactValue);
        coefficient.Value = 101;
        Assert.Equal(new BigRational(100), coefficient.ExactValue);
        coefficient.Value = -101;
        Assert.Equal(new BigRational(-100), coefficient.ExactValue);
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
        for (var i = 0; i < 12; i++) model.Coefficients[4].Value = i;
        model.ApplyPreset(CurvePreset.All[1]);
        await model.PendingSamples;
        Assert.NotEmpty(model.Samples);
        Assert.All(model.Samples, point => Assert.True(model.Snapshot.Curve.IsOnCurve(point)));
        Assert.Contains("|m| ≤ 12", model.SampleStatus);
        model.ApplyPreset(CurvePreset.All[3]);
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
        model.Dispose();
        model.Dispose();
        await pending;
        Assert.Empty(model.Samples);
    }
}

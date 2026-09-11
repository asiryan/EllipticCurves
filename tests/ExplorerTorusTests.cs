using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerTorusTests
{
    [Fact]
    public void LogarithmCoordinatesRespectBothPeriodsOfAShearedLattice()
    {
        // z = .2 omega1 + .7 omega2 with omega1 = 2, omega2 = 1 + 3i.
        var expected = new TorusCoordinates(0.2, 0.7);
        for (var a = -3; a <= 3; a++)
            for (var b = -3; b <= 3; b++)
            {
                var point = TorusCoordinates.FromLogarithm(1.1 + 2 * a + b, 2.1 + 3 * b, 2, 1, 3);
                Assert.Equal(expected.U, point.U, 12);
                Assert.Equal(expected.V, point.V, 12);
            }
    }

    [Fact]
    public void ClassicCurveMapsItsTwoTorsionToTheHalfPeriodPoints()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var lattice = TorusLattice.Create(curve, default);
        Assert.Equal(0, lattice.TauReal, 10);
        Assert.Equal(1, lattice.TauImaginary, 10);
        var mapped = lattice.MapPoints(new[]
        {
            new EllipticCurvePoint(-1, 0), new EllipticCurvePoint(0, 0), new EllipticCurvePoint(1, 0)
        }, default);
        Assert.Equal(4, mapped.Points.Count);
        Assert.Equal(0, mapped.SkippedCount);
        Assert.Equal(TorusPoint.Origin, mapped.Points[0]);
        AssertCoordinates(mapped.Points[1], 0, 0.5);
        AssertCoordinates(mapped.Points[2], 0.5, 0.5);
        AssertCoordinates(mapped.Points[3], 0.5, 0);
    }

    [Fact]
    public void GeneralEquationPointsKeepTheirExactOriginalCoordinates()
    {
        var curve = new EllipticCurveQ(0, 0, 1, -1, 0);
        var point = new EllipticCurvePoint(0, 0);
        var opposite = curve.Negate(point);
        var lattice = TorusLattice.Create(curve, default);
        var mapped = lattice.MapPoints(new[] { point, opposite, point }, default);
        Assert.Equal(2, mapped.SampleCount);
        Assert.Equal(3, mapped.Points.Count);
        Assert.Equal(point, mapped.Points[1].Point);
        Assert.Equal(opposite, mapped.Points[2].Point);
        var sumU = mapped.Points[1].Coordinates.U + mapped.Points[2].Coordinates.U;
        var sumV = mapped.Points[1].Coordinates.V + mapped.Points[2].Coordinates.V;
        Assert.True(Math.Abs(sumU - Math.Round(sumU)) < 1e-9);
        Assert.True(Math.Abs(sumV - Math.Round(sumV)) < 1e-9);
    }

    [Fact]
    public void OneRealComponentUsesAShearedPeriodBasis()
    {
        var curve = new EllipticCurveQ(0, 0, 0, 0, 1);
        var lattice = TorusLattice.Create(curve, default);
        Assert.Equal(0.5, lattice.TauReal, 10);
        Assert.True(lattice.TauImaginary > 0);
        var mapped = lattice.MapPoints(new[] { new EllipticCurvePoint(0, 1), new EllipticCurvePoint(0, -1) }, default);
        Assert.Equal(3, mapped.Points.Count);
        Assert.All(mapped.Points, point => Assert.Equal(0, point.Coordinates.V, 10));
    }

    [Fact]
    public void InvalidPointsAreReportedWithoutLosingTheOriginOrValidPoints()
    {
        var lattice = TorusLattice.Create(new EllipticCurveQ(0, 0, 0, -1, 0), default);
        var mapped = lattice.MapPoints(new[] { new EllipticCurvePoint(2, 0), new EllipticCurvePoint(1, 0) }, default);
        Assert.Equal(1, mapped.SkippedCount);
        Assert.Equal(2, mapped.Points.Count);
        Assert.Equal(2, mapped.SampleCount);
        Assert.Throws<OperationCanceledException>(() => lattice.MapPoints(new[] { new EllipticCurvePoint(1, 0) }, new CancellationToken(true)));
        Assert.Throws<ArgumentException>(() => TorusLattice.Create(new EllipticCurveQ(0, 0, 0, 0, 0), default));
        Assert.Throws<ArithmeticException>(() => TorusCoordinates.FromLogarithm(0, 0, 0, 0, 1));
    }

    [Fact]
    public async Task InactiveAndSingularViewsDoNotStartPeriodCalculations()
    {
        var calls = 0;
        using var model = new ComplexTorusViewModel((curve, token) => { calls++; return TorusLattice.Create(curve, token); });
        model.Update(new EllipticCurveQ(0, 0, 0, -1, 0), Array.Empty<EllipticCurvePoint>(), false);
        await model.PendingUpdate;
        Assert.Equal(0, calls);
        model.Update(new EllipticCurveQ(0, 0, 0, 0, 0), Array.Empty<EllipticCurvePoint>(), true);
        await model.PendingUpdate;
        Assert.Equal(0, calls);
        Assert.False(model.HasLattice);
        Assert.False(model.IsBusy);
        Assert.Contains("singular", model.Status);
    }

    [Fact]
    public async Task LateCalculationCannotRestoreAReplacedCurve()
    {
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var lattice = TorusLattice.Create(curve, default);
        using var started = new ManualResetEventSlim();
        using var release = new ManualResetEventSlim();
        using var model = new ComplexTorusViewModel((_, _) =>
        {
            started.Set();
            release.Wait(TimeSpan.FromSeconds(5));
            return lattice; // Deliberately finishes even after its cancellation token was cancelled.
        });
        model.Update(curve, Array.Empty<EllipticCurvePoint>(), true);
        var oldRequest = model.PendingUpdate;
        try
        {
            Assert.True(await Task.Run(() => started.Wait(TimeSpan.FromSeconds(5))));
            model.Update(new EllipticCurveQ(0, 0, 0, 0, 0), Array.Empty<EllipticCurvePoint>(), true);
        }
        finally { release.Set(); }
        await oldRequest.WaitAsync(TimeSpan.FromSeconds(5));
        Assert.False(model.HasLattice);
        Assert.Empty(model.Points);
        Assert.False(model.IsBusy);
        Assert.Contains("singular", model.Status);
    }

    private static void AssertCoordinates(TorusPoint point, double u, double v)
    {
        Assert.Equal(u, point.Coordinates.U, 9);
        Assert.Equal(v, point.Coordinates.V, 9);
    }
}

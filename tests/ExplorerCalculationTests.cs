using System.Numerics;
using System.Reflection;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerCalculationTests
{
    private const string Classic = "y^2 = x^3 - x";
    private static CalculationOperation Operation(string method, Type type = null) => CalculationCatalog.All.First(o =>
        o.Member is MethodInfo member && member.Name == method && member.DeclaringType == (type ?? typeof(EllipticCurveQ)));
    private static CalculationRequest Request(CalculationOperation operation, params (string Key, string Value)[] overrides)
    {
        var fields = operation.Parameters.ToDictionary(p => p.Key, p => p.Default);
        foreach (var (key, value) in overrides) fields[key] = value;
        return new(operation.Id, Classic, fields);
    }
    public static IEnumerable<object[]> OfflineOperations() => CalculationCatalog.All.Where(o => o.Context != CalculationContext.Database).Select(o => new object[] { o.Id });

    [Theory, MemberData(nameof(OfflineOperations))]
    public async Task EveryOfflineMenuActionRunsWithRepresentativeInputs(string id)
    {
        var operation = CalculationCatalog.Get(id);
        var request = Request(operation);
        if (operation.Member?.Name == "Regulator") request = request with { Equation = "y^2 + y = x^3 - x" };
        var output = await CalculationEngine.ExecuteAsync(request);
        Assert.False(string.IsNullOrWhiteSpace(output), id);
    }

    [Fact]
    public void EveryPublicCurveAndFieldMethodHasAnExecutableMenuEntry()
    {
        Assert.Equal(CalculationCatalog.All.Count, CalculationCatalog.All.Select(o => o.Id).Distinct().Count());
        foreach (var type in new[] { typeof(EllipticCurveQ), typeof(EllipticCurveFp), typeof(EllipticCurveFq), typeof(FiniteField), typeof(BigRational) })
            foreach (var method in type.GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly)
                .Where(CalculationCatalog.IsMathematicalMethod))
                Assert.Contains(CalculationCatalog.All, o => Equals(o.Member, method));
        foreach (var operation in CalculationCatalog.All)
        {
            Assert.Equal(operation.Parameters.Count, operation.Parameters.Select(p => p.Key).Distinct().Count());
            foreach (var field in operation.Parameters)
            {
                if (operation.Id == "database.import") continue;
                Assert.True(CalculationInput.Validate(field, field.Default) == "", operation.Id + ": " + field.Key);
            }
        }
        Assert.Contains(CalculationCatalog.All, o => o.Id == "Q.TorsionStructure");
        Assert.Contains(CalculationCatalog.All, o => o.Member?.Name == "StableFaltingsHeight");
        Assert.Contains(CalculationCatalog.All, o => o.Member?.Name == "EstimateAnalyticRank");
    }

    [Theory]
    [InlineData("8.325", "333/40")]
    [InlineData("8,325", "333/40")]
    [InlineData("-2/7", "-2/7")]
    [InlineData("1e-5", "1/100000")]
    public void InputsRetainExactRationals(string text, string exact)
        => Assert.Equal(exact, CalculationInput.ParseScalar(typeof(BigRational), text).ToString());

    [Fact]
    public async Task TorsionAndRankIncludeGroupAndProofInformation()
    {
        var torsion = await CalculationEngine.ExecuteAsync(Request(CalculationCatalog.Get("Q.TorsionStructure")));
        Assert.Contains("2", torsion);
        var bounds = await CalculationEngine.ExecuteAsync(Request(Operation("GetRankBounds")));
        Assert.Contains("Lower Bound: 0", bounds);
        Assert.Contains("Upper Bound: 0", bounds);
        Assert.Contains("Is Exact: True", bounds);
        Assert.Contains("Reason:", bounds);
    }

    [Fact]
    public async Task BothRankActionsDefaultToParallelExecutionAndAllowSequentialOverride()
    {
        Assert.Equal(1, new RankComputationOptions().MaxDegreeOfParallelism);
        foreach (var operation in CalculationCatalog.All.Where(o => o.Member?.Name == nameof(EllipticCurveQ.GetRankBounds)))
        {
            var degree = operation.Parameters.Single(p => p.Key.EndsWith(".MaxDegreeOfParallelism"));
            Assert.Equal(CalculationInput.DefaultRankWorkers.ToString(), degree.Default);
            Assert.InRange(int.Parse(degree.Default), 1, 4);
            Assert.Contains("1 runs sequentially", degree.Help);
            using var workbench = new WorkbenchViewModel();
            using var form = new CalculationFormViewModel(operation, "y^2 + y = x^3 - x", workbench);
            var field = form.Fields.Single(f => f.Parameter.Key == degree.Key);
            Assert.Equal(degree.Default, field.Text);
            field.Text = "0";
            Assert.False(form.CanRun);
            field.Text = "1";
            Assert.True(form.CanRun);
            var request = form.CreateRequest();
            var sequential = await CalculationEngine.ExecuteAsync(request);
            Assert.Contains("Exact Rank: 1", sequential);
            using var repeat = new CalculationFormViewModel(operation, request.Equation, workbench, request);
            Assert.Equal("1", repeat.Fields.Single(f => f.Parameter.Key == degree.Key).Text);
            var defaults = await CalculationEngine.ExecuteAsync(new(operation.Id, request.Equation, new()));
            Assert.Contains("Exact Rank: 1", defaults);
            var invalid = request with { Arguments = new(request.Arguments) { [degree.Key] = "0" } };
            await Assert.ThrowsAsync<FormatException>(() => CalculationEngine.ExecuteAsync(invalid));
        }
    }

    [Theory]
    [InlineData("FaltingsHeight")]
    [InlineData("StableFaltingsHeight")]
    [InlineData("CanonicalHeight")]
    [InlineData("GetPeriods")]
    public async Task CertifiedCalculationsIncludeExactEnclosures(string method)
    {
        var text = await CalculationEngine.ExecuteAsync(Request(Operation(method)));
        Assert.Contains("certified enclosure", text);
        Assert.Contains("Lower Bound:", text);
        Assert.Contains("Upper Bound:", text);
    }

    [Fact]
    public async Task OutParametersAndMapsAreUsable()
    {
        var request = Request(CalculationCatalog.Get("map.change.forward"), ("u", "2"));
        Assert.Contains("Mapped point: (0, 0)", await CalculationEngine.ExecuteAsync(request));
        var inverse = await CalculationEngine.ExecuteAsync(Request(CalculationCatalog.Get("map.change.inverse"), ("u", "2")));
        Assert.Contains("U: 1/2", inverse);
        Assert.Contains("Mapped point: O", await CalculationEngine.ExecuteAsync(Request(CalculationCatalog.Get("map.two"))));
        var isomorphism = await CalculationEngine.ExecuteAsync(Request(Operation("TryGetIsomorphism")));
        Assert.Contains("Succeeded: True", isomorphism);
        Assert.Contains("isomorphism:", isomorphism);
        var failure = Request(Operation("TryGetIsomorphism"), ("target", "y^2 = x^3 + 1"));
        Assert.Contains("isomorphism: unavailable", await CalculationEngine.ExecuteAsync(failure));
    }

    [Fact]
    public async Task PrimeCurveUsesEnteredCoefficientsAndSupportsInvalidMembership()
    {
        var count = await CalculationEngine.ExecuteAsync(Request(Operation("CountPoints", typeof(EllipticCurveFp))));
        Assert.Contains("Result: 8", count);
        var predicate = Request(Operation("IsOnCurve", typeof(EllipticCurveFp)), ("point.x", "0"), ("point.y", "1"));
        Assert.Contains("Result: False", await CalculationEngine.ExecuteAsync(predicate));
        var nonintegral = Request(Operation("CountPoints", typeof(EllipticCurveFp))) with { Equation = "y^2 = x^3 - x/16" };
        Assert.Contains("Result: 8", await CalculationEngine.ExecuteAsync(nonintegral));
    }

    [Fact]
    public async Task ExtensionElementsAndCurvesAcceptPolynomialCoordinates()
    {
        var pow = Request(Operation("Pow", typeof(FiniteField)), ("a", "0; 1"), ("exponent", "2"));
        Assert.Contains("[3]", await CalculationEngine.ExecuteAsync(pow));
        Assert.Contains("Result: 32", await CalculationEngine.ExecuteAsync(Request(Operation("CountPoints", typeof(EllipticCurveFq)))));
        var predicate = Request(Operation("IsOnCurve", typeof(EllipticCurveFq)), ("point.x", "0"), ("point.y", "1"));
        Assert.Contains("Result: False", await CalculationEngine.ExecuteAsync(predicate));
    }

    [Fact]
    public async Task RejectsOffCurveGroupInputsAndHandlesSingularOverview()
    {
        await Assert.ThrowsAsync<ArgumentException>(() => CalculationEngine.ExecuteAsync(Request(Operation("Double"), ("P.y", "1"))));
        var singular = Request(CalculationCatalog.Get("curve.overview")) with { Equation = "y^2 = x^3" };
        Assert.Contains("Singular: True", await CalculationEngine.ExecuteAsync(singular));
    }

    [Fact]
    public async Task TruncatedEnumerationCannotBeMistakenForACompleteList()
    {
        var request = Request(Operation("Points", typeof(EllipticCurveFp))) with { MaxItems = 2 };
        var text = await CalculationEngine.ExecuteAsync(request);
        Assert.Contains("OUTPUT TRUNCATED", text);
        Assert.Contains("[0]: O", text);
        Assert.DoesNotContain("[2]:", text);
    }

    [Fact]
    public async Task StoredDatabaseIncludesMetadataAndDistinguishesDecimalApproximations()
    {
        var json = File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "Fixtures", "lmfdb-37a1.json"));
        var text = await CalculationEngine.ExecuteAsync(Request(CalculationCatalog.Get("database.import"), ("json", json)));
        Assert.Contains("37.a1", text);
        Assert.Contains("Faltings Height", text);
        Assert.Contains("not a certified error interval", text);
        Assert.Contains("Isogeny Matrix", text);
    }

    [Fact]
    public void FormPreservesOptionDefaultsAndIgnoresCoordinatesAtInfinity()
    {
        using var workbench = new WorkbenchViewModel();
        using var form = new CalculationFormViewModel(Operation("CanonicalHeight"), Classic, workbench);
        Assert.Equal("256", form.Fields.Single(f => f.Parameter.Key == "options.PrecisionBits").Text);
        form.Fields.Single(f => f.Parameter.Key == "point.x").Text = "-";
        Assert.False(form.CanRun);
        form.Fields.Single(f => f.Parameter.Key == "point.infinity").Checked = true;
        Assert.True(form.CanRun);
        form.Timeout = "-1";
        Assert.False(form.CanRun);
        form.Timeout = "0";
        Assert.True(form.CanRun);
        form.MaxItems = "100001";
        Assert.False(form.CanRun);
    }

    [Fact]
    public void SearchFindsFunctionsAndReportsEmptyResults()
    {
        var menu = new ExplorerMenuViewModel { Search = "Faltings" };
        Assert.Contains(menu.Operations, o => o.Member?.Name == "FaltingsHeight");
        menu.Search = "not-a-real-operation";
        Assert.Empty(menu.Operations);
    }
}

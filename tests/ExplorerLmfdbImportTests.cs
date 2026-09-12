using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ExplorerLmfdbImportTests
{
    private static readonly LmfdbCurveFormula Formula = new("37.a1", 37, "y^2 + y = x^3 - x");

    [Theory]
    [InlineData("37", 37, 37)]
    [InlineData("11-100", 11, 100)]
    [InlineData("1 – 500 000", 1, 500000)]
    [InlineData("500000", 500000, 500000)]
    public void ParsesConductorAndInclusiveRange(string text, int minimum, int maximum)
    {
        Assert.True(LmfdbConductorRange.TryParse(text, out var range));
        Assert.Equal(new(minimum, maximum), range);
    }

    [Theory]
    [InlineData("")][InlineData("0")][InlineData("-37")][InlineData("100-11")]
    [InlineData("500001")][InlineData("1-500001")][InlineData("37.a1")]
    [InlineData("1-2-3")][InlineData("1e3")][InlineData("37,38")]
    public void RejectsInvalidConductorInput(string text) => Assert.False(LmfdbConductorRange.TryParse(text, out _));

    [Fact]
    public async Task LoadsOnlyIdentifiersAndExactEquationsFromRealApiFixture()
    {
        var handler = new LmfdbTestHttpHandler(await File.ReadAllTextAsync(Path.Combine(AppContext.BaseDirectory, "Fixtures", "lmfdb-formulas-37.json")));
        using var http = new HttpClient(handler);
        var curves = await new LmfdbCurveSearch(http).SearchAsync(new(37, 37));
        Assert.Equal(new[] { "37.a1", "37.b1", "37.b2", "37.b3" }, curves.Select(row => row.Label));
        Assert.Equal(Formula, curves[0]);
        Assert.All(curves, row => Assert.True(CurveEquationText.TryParse(row.Equation, out _, out _)));
        var url = Assert.Single(handler.Requests);
        Assert.Contains("_fields=lmfdb_label,ainvs", url);
        Assert.Contains("conductor=i37", url);
        Assert.DoesNotContain("rank", url);
        Assert.DoesNotContain("torsion", url);
        Assert.DoesNotContain("/data/", url);
    }

    [Fact]
    public async Task UsesRangeAndSeekCursorWithoutTheApiOffsetLimit()
    {
        var handler = new LmfdbTestHttpHandler("{\"data\":[]}");
        using var captured = new HttpClient(handler);
        await new LmfdbCurveSearch(captured).SearchAsync(new(11, 500000), Formula);
        var url = Uri.UnescapeDataString(Assert.Single(handler.Requests));
        Assert.Contains("py{\"$gte\":11,\"$lte\":500000}", url);
        Assert.Contains("$or=py[{\"conductor\":{\"$gt\":37}},{\"conductor\":37,\"lmfdb_label\":{\"$gt\":\"37.a1\"}}]", url);
        Assert.DoesNotContain("_offset", url);
    }

    [Fact]
    public async Task CoefficientsNeverPassThroughFloatingPoint()
    {
        const string coefficient = "9007199254740993";
        using var http = new HttpClient(new LmfdbTestHttpHandler("{\"data\":[{\"lmfdb_label\":\"37.a1\",\"ainvs\":[0,0,1," + coefficient + ",0]}]}"));
        var result = Assert.Single(await new LmfdbCurveSearch(http).SearchAsync(new(37, 37)));
        Assert.Contains(coefficient, result.Equation);
    }

    [Theory]
    [InlineData("{\"data\":[{\"lmfdb_label\":\"38.a1\",\"ainvs\":[0,0,1,-1,0]}]}")]
    [InlineData("{\"data\":[{\"lmfdb_label\":\"37.a1\",\"ainvs\":[0,0,0,0,0]}]}")]
    [InlineData("{\"data\":[{\"lmfdb_label\":\"37.a1\",\"ainvs\":[0,0,-1,0]}]}")]
    [InlineData("{\"data\":[{\"lmfdb_label\":\"37.a1\",\"ainvs\":[0,0,1,1.5,0]}]}")]
    [InlineData("{\"data\":null}")]
    public async Task RejectsMalformedOrUnrelatedFormulas(string json)
    {
        using var http = new HttpClient(new LmfdbTestHttpHandler(json));
        await Assert.ThrowsAsync<FormatException>(() => new LmfdbCurveSearch(http).SearchAsync(new(37, 37)));
    }

    [Fact]
    public async Task RejectsRepeatedPages()
    {
        using var http = new HttpClient(new LmfdbTestHttpHandler("{\"data\":[{\"lmfdb_label\":\"37.a1\",\"ainvs\":[0,0,1,-1,0]}]}"));
        await Assert.ThrowsAsync<FormatException>(() => new LmfdbCurveSearch(http).SearchAsync(new(37, 37), Formula));
    }

    [Fact]
    public async Task CaptchaPageLeavesSearchAvailableWithoutImportableResults()
    {
        var challenge = new HttpResponseMessage(System.Net.HttpStatusCode.OK)
            { Content = new StringContent("<html>CAPTCHA</html>", System.Text.Encoding.UTF8, "text/html") };
        using var http = new HttpClient(new ResponseHandler(challenge));
        using var model = new LmfdbImportViewModel(new LmfdbCurveSearch(http).SearchAsync);
        await model.SearchAsync();
        Assert.Contains("CAPTCHA", model.Error);
        Assert.True(model.CanSearch);
        Assert.False(model.CanImport);
    }

    private sealed class ResponseHandler(HttpResponseMessage response) : HttpMessageHandler
    {
        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            => Task.FromResult(response);
    }

    [Fact]
    public async Task PagingRetainsTheSelectedPageOnNetworkFailureAndCanRetry()
    {
        var firstPage = Enumerable.Range(0, 100).Select(i => Formula with { Label = "37.a" + (i + 1) }).ToArray();
        var calls = new List<LmfdbCurveFormula>();
        var fail = false;
        using var model = new LmfdbImportViewModel((range, cursor, token) =>
        {
            calls.Add(cursor);
            if (fail) throw new HttpRequestException("Offline");
            return Task.FromResult<IReadOnlyList<LmfdbCurveFormula>>(cursor == null ? firstPage : new[] { Formula with { Label = "37.b1" } });
        });
        await model.SearchAsync();
        Assert.True(model.CanNext);
        Assert.True(model.CanImport);
        fail = true;
        await model.NextAsync();
        Assert.Same(firstPage, model.Curves);
        Assert.False(model.CanPrevious);
        Assert.NotEmpty(model.Error);
        fail = false;
        await model.NextAsync();
        Assert.True(model.CanPrevious);
        Assert.False(model.CanNext);
        Assert.Equal("37.b1", model.Selected.Label);
        await model.PreviousAsync();
        Assert.Same(firstPage, model.Curves);
        Assert.Null(calls[^1]);
    }

    [Fact]
    public async Task ChangingQueryDiscardsLateResultsEvenIfServerIgnoresCancellation()
    {
        var response = new TaskCompletionSource<IReadOnlyList<LmfdbCurveFormula>>();
        using var model = new LmfdbImportViewModel((range, cursor, token) => response.Task);
        var pending = model.SearchAsync();
        Assert.True(model.IsBusy);
        model.Query = "11";
        response.SetResult(new[] { Formula });
        await pending;
        Assert.Empty(model.Curves);
        Assert.False(model.CanImport);
        Assert.False(model.IsBusy);
    }

    [Fact]
    public async Task StopAndCloseDoNotAcceptLateResults()
    {
        foreach (var close in new[] { false, true })
        {
            var response = new TaskCompletionSource<IReadOnlyList<LmfdbCurveFormula>>();
            using var model = new LmfdbImportViewModel((range, cursor, token) => response.Task);
            var pending = model.SearchAsync();
            if (close) model.Dispose(); else model.Cancel();
            response.SetResult(new[] { Formula });
            await pending;
            Assert.Empty(model.Curves);
            Assert.False(model.CanImport);
            Assert.False(model.IsBusy);
        }
    }

    [Fact]
    public async Task EmptyResultAndInvalidQueryCannotBeImported()
    {
        var calls = 0;
        using var model = new LmfdbImportViewModel((range, cursor, token) =>
        {
            calls++;
            return Task.FromResult<IReadOnlyList<LmfdbCurveFormula>>(Array.Empty<LmfdbCurveFormula>());
        });
        await model.SearchAsync();
        Assert.Contains("No curves found", model.Status);
        Assert.False(model.CanImport);
        model.Query = "500001";
        await model.SearchAsync();
        Assert.Equal(1, calls);
        Assert.NotEmpty(model.InputError);
    }

    [Fact]
    public void ImportIsDiscoverableInTheLmfdbCategoryAndSearch()
    {
        var model = new ExplorerMenuViewModel { Group = "LMFDB · internet" };
        Assert.Contains(ExplorerMenuViewModel.ImportCurve, model.Operations);
        model.Search = "conductor";
        Assert.Contains(ExplorerMenuViewModel.ImportCurve, model.Operations);
        model.Search = "Faltings";
        Assert.DoesNotContain(ExplorerMenuViewModel.ImportCurve, model.Operations);
        model.Search = "";
        model.Group = "Prime fields";
        Assert.DoesNotContain(ExplorerMenuViewModel.ImportCurve, model.Operations);
    }
}

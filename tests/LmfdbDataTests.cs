using System.Net;
using System.Numerics;
using System.Text.Json;
using System.Text.Json.Nodes;
using EllipticCurves;
using Xunit;
using static EllipticCurves.Tests.ExtendedReferenceTests;

namespace EllipticCurves.Tests;

public class LmfdbDataTests
{
    internal static string Fixture(string label = "37a1") => File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "Fixtures", "lmfdb-" + label + ".json"));
    private static JsonObject Snapshot() => JsonNode.Parse(Fixture()).AsObject();
    private static JsonObject Curve(JsonObject snapshot) => snapshot["data"][0][0].AsObject();
    private static JsonObject Mw(JsonObject snapshot) => snapshot["data"][2][0].AsObject();

    [Theory]
    [InlineData("37a1"), InlineData("48a3"), InlineData("389a1"), InlineData("5077a1"), InlineData("11a1")]
    public void StoredDataMatchesNativeArithmetic(string label)
    {
        var data = LmfdbEllipticCurve.FromStoredDataJson(Fixture(label)); var e = data.GlobalMinimalModel;
        Assert.Equal(e.Conductor, data.Conductor); Assert.Equal(e.TorsionStructure, data.TorsionStructure);
        var local = e.LocalData; Assert.Equal(local.Count, data.LocalData.Count);
        for (int i = 0; i < local.Count; i++)
        {
            var a = local[i]; var b = data.LocalData[i];
            Assert.Equal(a.Prime, b.Prime); Assert.Equal(a.KodairaSymbol, b.KodairaSymbol);
            Assert.Equal(a.TamagawaNumber, b.TamagawaNumber); Assert.Equal(a.ReductionType, b.ReductionType);
            Assert.Equal(a.ConductorValuation, b.ConductorValuation); Assert.Equal(a.DiscriminantValuation, b.DiscriminantValuation);
            Assert.Equal(a.JDenominatorValuation, b.JDenominatorValuation); Assert.Equal(a.RootNumber, b.RootNumber);
        }
        Assert.Equal(e.TamagawaProduct, data.TamagawaProduct);
        var periods = e.GetPeriods(); StoredNear(periods.RealPeriod, data.RealPeriod); StoredNear(periods.Area, data.PeriodArea);
        for (int i = 0; i < data.MordellWeilGenerators.Count; i++) StoredNear(e.CanonicalHeight(data.MordellWeilGenerators[i]), data.GeneratorHeights[i]);
        StoredNear(e.Regulator(data.MordellWeilGenerators), data.Regulator);
        Assert.Equal(data.Rank, data.RankBounds.LowerBound); Assert.Equal(data.Rank, data.RankBounds.UpperBound);
        Assert.All(data.TorsionGenerators, p => Assert.Contains(p, e.TorsionPoints));
        var changed = e.ChangeModel(new BigRational(2, 3), 5, -2, 7).Target;
        Assert.All(data.GetGeneratorsOnModel(changed), p => Assert.True(changed.IsOnCurve(p)));
    }
    // Stored decimal accuracy is not certified by LMFDB. These checks only compare
    // the supplied approximations with a generous tolerance below the native target.
    private static void StoredNear(RealEnclosure native, LmfdbRealValue stored)
    {
        Assert.NotNull(stored); var error = new BigRational(1, BigInteger.Pow(10, 20)); var value = stored.AsRational();
        Assert.True(native.LowerBound <= value + error && native.UpperBound >= value - error, $"Stored {value} versus native {native}.");
    }
    [Fact]
    public void MissingFieldsRemainNullAndStoredPrecisionIsNotLost()
    {
        var snapshot = Snapshot(); var mw = Mw(snapshot);
        foreach (var key in new[] { "gens", "heights", "area", "rank_bounds", "real_period" }) mw.Remove(key);
        var data = LmfdbEllipticCurve.FromStoredDataJson(snapshot.ToJsonString());
        Assert.Null(data.MordellWeilGenerators); Assert.Null(data.GeneratorHeights); Assert.Null(data.PeriodArea); Assert.Null(data.RankBounds); Assert.Null(data.RealPeriod);
        Assert.Empty(data.TorsionGenerators); Assert.Equal(97, data.Regulator.StoredPrecisionBits);
        Curve(snapshot)["regulator"] = JsonValue.Create("1.234e-30");
        data = LmfdbEllipticCurve.FromStoredDataJson(snapshot.ToJsonString());
        Assert.Equal(new BigRational(1234, BigInteger.Pow(10, 33)), data.Regulator.AsRational());
        Assert.Null(data.Regulator.StoredPrecisionBits);
    }
    [Fact]
    public void GeneratorCoordinatesAreWeightedProjective()
    {
        var snapshot = Snapshot(); var e = new EllipticCurveQ(0, 0, 1, -1, 0); var p = e.Multiply(new EllipticCurvePoint(0, 0), 5);
        Assert.True(p.X.Den > 1); var z = BigInteger.Parse("2");
        // 5P=(1/4,-5/8), so the stored triple is [1,-5,2].
        Assert.Equal(new BigRational(1, 4), p.X);
        Mw(snapshot)["gens"] = new JsonArray(new JsonArray(JsonValue.Create((p.X * new BigRational(z * z)).Num.ToString()), JsonValue.Create((p.Y * new BigRational(z * z * z)).Num.ToString()), JsonValue.Create(z.ToString())));
        var data = LmfdbEllipticCurve.FromStoredDataJson(snapshot.ToJsonString()); Assert.Equal(p, data.MordellWeilGenerators[0]);
    }
    [Fact]
    public void MalformedAndIncompleteSnapshotsAreRejected()
    {
        var s = Snapshot(); s["totals"][5] = 2; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); Mw(s)["lmfdb_label"] = "11.a1"; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); Mw(s)["gens"] = new JsonArray(new JsonArray(0, 1, 1)); Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); Mw(s)["tamagawa_product"] = 9; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["data"][5] = new JsonArray(); s["totals"][5] = 0; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
    }
    [Fact]
    public async Task AsyncLookupPaginatesAndVerifiesIsomorphism()
    {
        // A quadratic twist with exactly the same j-invariant must be rejected.
        var wrong = new { ainvs = new[] { 0, 0, 0, -4, 2 }, lmfdb_label = "wrong" };
        var first = JsonSerializer.Serialize(new { data = Enumerable.Repeat(wrong, 100).ToArray() });
        var second = "{\"data\":[{\"ainvs\":[0,0,1,-1,0],\"lmfdb_label\":\"37.a1\"}]}";
        using var handler = new LmfdbTestHttpHandler(first, second, Fixture()); using var client = new HttpClient(handler);
        var e = new EllipticCurveQ(0, 0, 1, -1, 0).ChangeModel(2, 3, -2, 5).Target;
        var data = await LmfdbEllipticCurve.FetchAsync(e, client);
        Assert.Equal("37.a1", data.Label); Assert.Equal(3, handler.Requests.Count);
        Assert.Contains("_offset=100", handler.Requests[1]); Assert.Contains("/data/37.a1?", handler.Requests[2]);
        Assert.NotNull(data.RealPeriod); _ = data.LocalData; _ = data.Regulator; Assert.Equal(3, handler.Requests.Count);
        // The caller can continue to use its client after fetching.
        Assert.Equal(HttpStatusCode.OK, (await client.GetAsync("https://www.lmfdb.org/test")).StatusCode);
    }
    [Fact]
    public async Task LookupReportsMissingCurveHttpFailureCaptchaAndCancellation()
    {
        var e = new EllipticCurveQ(0, 0, 1, -1, 0);
        using var missing = new HttpClient(new LmfdbTestHttpHandler("{\"data\":[]}"));
        await Assert.ThrowsAsync<InvalidOperationException>(() => LmfdbEllipticCurve.FetchAsync(e, missing));
        using var http = new HttpClient(new LmfdbTestHttpHandler("failure") { Status = HttpStatusCode.ServiceUnavailable });
        await Assert.ThrowsAsync<HttpRequestException>(() => LmfdbEllipticCurve.FetchAsync(e, http));
        using var captcha = new HttpClient(new LmfdbTestHttpHandler("<html>CAPTCHA challenge</html>"));
        await Assert.ThrowsAsync<InvalidOperationException>(() => LmfdbEllipticCurve.FetchAsync(e, captcha));
        using var handler = new LmfdbTestHttpHandler(Fixture()); using var canceled = new HttpClient(handler);
        await Assert.ThrowsAnyAsync<OperationCanceledException>(() => LmfdbEllipticCurve.FetchAsync(e, canceled, new CancellationToken(true)));
        Assert.Empty(handler.Requests);
    }
}

using System.Numerics;
using System.Text.Json.Nodes;
using EllipticCurves;
using Xunit;

namespace EllipticCurves.Tests;

public class LmfdbMetadataTests
{
    private static JsonObject Snapshot() => JsonNode.Parse(LmfdbDataTests.Fixture()).AsObject();
    [Fact]
    public void OptionalMetadataRemainsOptional()
    {
        var snapshot = Snapshot(); var row = snapshot["data"][0][0].AsObject(); var cls = snapshot["data"][1][0].AsObject();
        foreach (string key in new[] { "Clabel", "cm", "degree", "manin_constant", "sha", "torsion", "faltings_height", "stable_faltings_height", "isogeny_degrees", "class_size" }) row.Remove(key);
        foreach (string key in new[] { "anlist", "aplist", "class_size", "isogeny_matrix" }) cls.Remove(key);
        var mw = snapshot["data"][2][0].AsObject(); foreach (string key in new[] { "sha_an", "special_value", "xcoord_integral_points" }) mw.Remove(key);
        var data = LmfdbEllipticCurve.FromStoredDataJson(snapshot.ToJsonString());
        Assert.Null(data.CmDiscriminant); Assert.Null(data.ModularDegree); Assert.Null(data.ManinConstant); Assert.Null(data.BsdShaOrder);
        Assert.Null(data.TorsionOrder); Assert.Null(data.FaltingsHeight); Assert.Null(data.StableFaltingsHeight); Assert.Null(data.IsogenyDegrees);
        Assert.Null(data.IsogenyClassSize); Assert.Null(data.FourierCoefficients); Assert.Null(data.PrimeFourierCoefficients); Assert.Null(data.IsogenyMatrix);
        Assert.Null(data.AnalyticShaOrder); Assert.Null(data.LeadingLValue); Assert.Null(data.IntegralPointXCoordinates);
        Assert.NotNull(data.RealPeriod);
    }

    [Fact]
    public void ExtendedMetadataIsValidatedAndImmutable()
    {
        var s = Snapshot(); s["data"][1][0]["lmfdb_iso"] = "11.a";
        Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["data"][0][0].AsObject().Remove("lmfdb_iso"); s["data"][1][0]["lmfdb_iso"] = "11.a";
        Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["totals"][1] = 2; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["data"][1][0]["isogeny_matrix"] = new JsonArray(new JsonArray(1, 2));
        Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["data"][1][0]["anlist"][0] = 1; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        s = Snapshot(); s["data"][0][0]["degree"] = -1; Assert.Throws<FormatException>(() => LmfdbEllipticCurve.FromStoredDataJson(s.ToJsonString()));
        var data = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture());
        Assert.Throws<NotSupportedException>(() => ((IList<BigInteger>)data.FourierCoefficients)[1] = 7);
        Assert.Throws<NotSupportedException>(() => ((IList<int>)data.IsogenyMatrix[0])[0] = 7);
        Assert.Equal("37a1", data.CremonaLabel); Assert.Equal("37.a", data.IsogenyClassLabel);
        Assert.Equal(2, data.ModularDegree); Assert.Equal(1, data.BsdShaOrder); Assert.Equal(0, data.CmDiscriminant);
    }

    [Theory]
    [InlineData("389a1"), InlineData("5077a1")]
    public void StoredLeadingValueUsesTheRankFactorial(string label)
    {
        var data = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture(label));
        var result = data.GlobalMinimalModel.EstimateAnalyticRank(new AnalyticRankOptions { CertifyLowRanks = false });
        double factorial = Enumerable.Range(1, data.Rank).Aggregate(1.0, (a, b) => a * b);
        Assert.NotNull(result.EstimatedRank);
        Assert.True(Math.Abs(result.Derivatives[data.Rank] / factorial - data.LeadingLValue.Approximation) < 1e-8);
    }
}

using Xunit;

namespace EllipticCurves.Tests;

public sealed class LmfdbCancellationTests
{
    [Fact]
    public void GeneratorMappingHonorsCancellationAndKeepsExactCoordinates()
    {
        var data = LmfdbEllipticCurve.FromStoredDataJson(LmfdbDataTests.Fixture());
        var model = data.GlobalMinimalModel.ChangeModel(2, 1, 1, 1).Target;
        var token = new CancellationToken(true);
        var error = Assert.ThrowsAny<OperationCanceledException>(() => data.GetGeneratorsOnModel(model, token));
        Assert.Equal(token, error.CancellationToken);

        var mapped = data.GetGeneratorsOnModel(model, default);
        Assert.Equal(data.GetGeneratorsOnModel(model), mapped);
        Assert.Equal(data.MordellWeilGenerators.Count, mapped.Count);
        Assert.All(mapped, point => Assert.True(model.IsOnCurve(point)));
    }
}

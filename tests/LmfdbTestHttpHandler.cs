using System.Net;

namespace EllipticCurves.Tests;

internal sealed class LmfdbTestHttpHandler(params string[] responses) : HttpMessageHandler
{
    public readonly List<string> Requests = new();
    public HttpStatusCode Status = HttpStatusCode.OK;
    protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken token)
    {
        token.ThrowIfCancellationRequested(); Requests.Add(request.RequestUri.ToString());
        return Task.FromResult(new HttpResponseMessage(Status) { Content = new StringContent(responses[Math.Min(Requests.Count - 1, responses.Length - 1)]) });
    }
}

#nullable enable
using System.Globalization;
using System.IO;
using System.Net.Http;
using System.Numerics;
using System.Text.Json;
using System.Text.RegularExpressions;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Computations;

public sealed class LmfdbCurveSearch(HttpClient? client = null)
{
    public const int PageSize = 100;
    private const int MaxResponseBytes = 2_000_000;
    private static readonly HttpClient SharedClient = new() { Timeout = TimeSpan.FromSeconds(30) };
    private static readonly Regex LabelPattern = new(@"^([1-9][0-9]{0,5})\.[a-z]+[1-9][0-9]*$", RegexOptions.CultureInvariant);
    private readonly HttpClient http = client ?? SharedClient;

    public async Task<IReadOnlyList<LmfdbCurveFormula>> SearchAsync(LmfdbConductorRange range,
        LmfdbCurveFormula? after = null, CancellationToken cancellationToken = default)
    {
        if (range.Minimum < 1 || range.Maximum < range.Minimum || range.Maximum > LmfdbConductorRange.Limit)
            throw new ArgumentOutOfRangeException(nameof(range));
        var filter = range.Minimum == range.Maximum ? "i" + range.Minimum.ToString(CultureInfo.InvariantCulture)
            : "py" + JsonSerializer.Serialize(new Dictionary<string, int> { ["$gte"] = range.Minimum, ["$lte"] = range.Maximum });
        var url = "https://www.lmfdb.org/api/ec_curvedata/?_format=json&_fields=lmfdb_label,ainvs"
            + "&_sort=conductor,lmfdb_label&conductor=" + Uri.EscapeDataString(filter);
        if (after != null)
        {
            if (ReadConductor(after.Label) != after.Conductor || after.Conductor < range.Minimum || after.Conductor > range.Maximum)
                throw new ArgumentException("The page cursor does not belong to this conductor range.", nameof(after));
            // Seek after the last (conductor, label), rather than using API offsets,
            // which LMFDB caps at 10,000. Each request still loads at most 100 formulas.
            var cursor = "py" + JsonSerializer.Serialize(new object[]
            {
                new { conductor = new Dictionary<string, int> { ["$gt"] = after.Conductor } },
                new { conductor = after.Conductor, lmfdb_label = new Dictionary<string, string> { ["$gt"] = after.Label } }
            });
            url += "&%24or=" + Uri.EscapeDataString(cursor);
        }
        using var timeout = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
        timeout.CancelAfter(TimeSpan.FromSeconds(30));
        var token = timeout.Token;
        using var request = new HttpRequestMessage(HttpMethod.Get, url);
        request.Headers.Accept.ParseAdd("application/json");
        request.Headers.UserAgent.ParseAdd("EllipticCurves.Explorer/" + typeof(LmfdbCurveSearch).Assembly.GetName().Version!.ToString(3));
        using var response = await http.SendAsync(request, HttpCompletionOption.ResponseHeadersRead, token).ConfigureAwait(false);
        response.EnsureSuccessStatusCode();
        if (response.Content.Headers.ContentType?.MediaType?.Contains("html", StringComparison.OrdinalIgnoreCase) == true)
            throw new FormatException("LMFDB returned a web page instead of formulas. It may be asking for a CAPTCHA; try again later.");
        if (response.Content.Headers.ContentLength > MaxResponseBytes) throw new FormatException("The LMFDB response is too large.");
        using var stream = await response.Content.ReadAsStreamAsync(token).ConfigureAwait(false);
        using var buffer = new MemoryStream();
        var bytes = new byte[8192];
        int read;
        while ((read = await stream.ReadAsync(bytes.AsMemory(), token).ConfigureAwait(false)) != 0)
        {
            if (buffer.Length + read > MaxResponseBytes) throw new FormatException("The LMFDB response is too large.");
            buffer.Write(bytes, 0, read);
        }
        token.ThrowIfCancellationRequested();
        using var json = JsonDocument.Parse(buffer.ToArray());
        if (!json.RootElement.TryGetProperty("data", out var data) || data.ValueKind != JsonValueKind.Array || data.GetArrayLength() > PageSize)
            throw new FormatException("LMFDB returned an invalid curve list.");
        var result = new List<LmfdbCurveFormula>();
        var previous = after;
        foreach (var row in data.EnumerateArray())
        {
            token.ThrowIfCancellationRequested();
            var label = row.GetProperty("lmfdb_label").GetString() ?? "";
            var conductor = ReadConductor(label);
            if (conductor < range.Minimum || conductor > range.Maximum
                || (previous != null && (conductor < previous.Conductor
                    || (conductor == previous.Conductor && string.CompareOrdinal(label, previous.Label) <= 0))))
                throw new FormatException("LMFDB returned formulas outside the requested page.");
            var ainvs = row.GetProperty("ainvs");
            if (ainvs.ValueKind != JsonValueKind.Array || ainvs.GetArrayLength() != 5)
                throw new FormatException("LMFDB did not return five curve coefficients.");
            var values = ainvs.EnumerateArray().Select(ReadInteger).Select(value => new BigRational(value)).ToArray();
            var curve = new EllipticCurveQ(values[0], values[1], values[2], values[3], values[4]);
            if (curve.Discriminant.IsZero) throw new FormatException("LMFDB returned a singular curve.");
            var equation = CurveEquationText.Format(curve);
            if (!CurveEquationText.TryParse(equation, out _, out _)) throw new FormatException("The returned formula exceeds the editor's limits.");
            previous = new(label, conductor, equation);
            result.Add(previous);
        }
        return result.AsReadOnly();
    }

    private static int ReadConductor(string label)
    {
        if (label.Length > 64) throw new FormatException("LMFDB returned an invalid curve label.");
        var match = LabelPattern.Match(label);
        if (!match.Success) throw new FormatException("LMFDB returned an invalid curve label.");
        return int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
    }

    private static BigInteger ReadInteger(JsonElement value)
    {
        var text = value.ValueKind == JsonValueKind.String ? value.GetString()! : value.GetRawText();
        if (text.Length > 4096 || !BigInteger.TryParse(text, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out var integer))
            throw new FormatException("LMFDB returned an invalid integer coefficient.");
        return integer;
    }
}

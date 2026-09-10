using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Net.Http;
using System.Numerics;
using System.Text.Json;
using System.Threading;
using System.Threading.Tasks;

namespace EllipticCurves
{
    /// <summary>
    /// Cached LMFDB data for a curve over Q. Candidate lookup verifies exact Q-isomorphism;
    /// a second request loads the curve's table snapshot. Properties never make network calls.
    /// Missing optional database fields remain null, and stored real values are approximations.
    /// </summary>
    public sealed partial class LmfdbEllipticCurve
    {
        private const string Origin = "https://www.lmfdb.org";
        private static readonly HttpClient SharedClient = new HttpClient { Timeout = TimeSpan.FromSeconds(30) };

        /// <summary>Synchronously fetches the matching curve and its data. Prefer FetchAsync in asynchronous applications.</summary>
        public LmfdbEllipticCurve(EllipticCurveQ ellipticCurve)
            : this(FetchSnapshotAsync(ellipticCurve, SharedClient, CancellationToken.None).GetAwaiter().GetResult()) { }

        /// <summary>Fetches and caches data. A supplied HttpClient remains owned by the caller.</summary>
        public static async Task<LmfdbEllipticCurve> FetchAsync(EllipticCurveQ ellipticCurve, HttpClient httpClient = null, CancellationToken cancellationToken = default)
            => new LmfdbEllipticCurve(await FetchSnapshotAsync(ellipticCurve, httpClient ?? SharedClient, cancellationToken).ConfigureAwait(false));

        /// <summary>Reads an offline snapshot from /EllipticCurve/Q/data/{label}?_format=json.</summary>
        public static LmfdbEllipticCurve FromStoredDataJson(string json)
        {
            using var document = JsonDocument.Parse(json ?? throw new ArgumentNullException(nameof(json)));
            return new LmfdbEllipticCurve(document.RootElement);
        }

        /// <summary>Algebraic rank recorded by LMFDB.</summary>
        public int Rank { get; }
        /// <summary>Analytic rank, when recorded.</summary>
        public int? AnalyticRank { get; }
        /// <summary>The conductor.</summary>
        public BigInteger Conductor { get; }
        /// <summary>The curve label, for example 37.a1.</summary>
        public string Label { get; }
        /// <summary>The curve's human-facing page.</summary>
        public string Url => Origin + "/EllipticCurve/Q/" + Uri.EscapeDataString(Label) + "/";
        /// <summary>Torsion structure, for example Z/2Z x Z/4Z.</summary>
        public string TorsionStructure { get; }
        /// <summary>The model used for all stored points and invariants.</summary>
        public EllipticCurveQ GlobalMinimalModel { get; }
        /// <summary>Local data at all bad primes, or null when the table is absent.</summary>
        public IReadOnlyList<LocalReductionData> LocalData { get; }
        /// <summary>The product of local Tamagawa numbers, when recorded.</summary>
        public BigInteger? TamagawaProduct { get; }
        /// <summary>Stored free generators, on GlobalMinimalModel, or null if absent.</summary>
        public IReadOnlyList<EllipticCurvePoint> MordellWeilGenerators { get; }
        /// <summary>Stored torsion generators, on GlobalMinimalModel, or null if absent.</summary>
        public IReadOnlyList<EllipticCurvePoint> TorsionGenerators { get; }
        /// <summary>Canonical heights in the same order as MordellWeilGenerators, when recorded.</summary>
        public IReadOnlyList<LmfdbRealValue> GeneratorHeights { get; }
        /// <summary>Stored rank bounds, when present.</summary>
        public LmfdbRankBounds RankBounds { get; }
        /// <summary>Stored regulator approximation.</summary>
        public LmfdbRealValue Regulator { get; }
        /// <summary>Stored BSD real period, including the number of real components.</summary>
        public LmfdbRealValue RealPeriod { get; }
        /// <summary>Stored area of the complex period lattice.</summary>
        public LmfdbRealValue PeriodArea { get; }

        /// <summary>Maps the stored free generators to an exactly Q-isomorphic model; null if generators are absent.</summary>
        public IReadOnlyList<EllipticCurvePoint> GetGeneratorsOnModel(EllipticCurveQ model)
        {
            if (model == null) throw new ArgumentNullException(nameof(model));
            if (!GlobalMinimalModel.TryGetIsomorphism(model, out var map)) throw new ArgumentException("Models are not Q-isomorphic.", nameof(model));
            return MordellWeilGenerators == null ? null : Array.AsReadOnly(MordellWeilGenerators.Select(map.Map).ToArray());
        }

        private LmfdbEllipticCurve(JsonElement root)
        {
            var curves = Table(root, "ec_curvedata");
            if (!curves.HasValue || curves.Value.GetArrayLength() != 1) throw new FormatException("LMFDB: expected exactly one curve record.");
            var row = curves.Value[0];
            Label = row.GetProperty("lmfdb_label").GetString();
            if (string.IsNullOrWhiteSpace(Label)) throw new FormatException("LMFDB: missing label.");
            GlobalMinimalModel = Model(row.GetProperty("ainvs"));
            if (GlobalMinimalModel.Discriminant.IsZero) throw new FormatException("LMFDB: singular curve.");
            Conductor = Integer(row.GetProperty("conductor"));
            Rank = checked((int)Integer(row.GetProperty("rank")));
            AnalyticRank = Optional(row, "analytic_rank") is JsonElement ar ? checked((int)Integer(ar)) : (int?)null;
            if (Conductor <= 0 || Rank < 0 || AnalyticRank < 0) throw new FormatException("LMFDB: invalid conductor or rank.");
            var torsion = row.GetProperty("torsion_structure");
            TorsionStructure = torsion.ValueKind == JsonValueKind.String ? torsion.GetString() :
                torsion.GetArrayLength() == 0 ? "Z/1Z" : string.Join(" x ", torsion.EnumerateArray().Select(x => "Z/" + Integer(x).ToString(CultureInfo.InvariantCulture) + "Z"));
            Regulator = OptionalReal(row, "regulator");

            var mw = Table(root, "ec_mwbsd");
            if (mw.HasValue && mw.Value.GetArrayLength() > 1) throw new FormatException("LMFDB: multiple Mordell-Weil records.");
            if (mw.HasValue && mw.Value.GetArrayLength() == 1)
            {
                var m = mw.Value[0]; CheckLabel(m);
                MordellWeilGenerators = Points(m, "gens");
                TorsionGenerators = Points(m, "torsion_generators");
                if (Optional(m, "heights") is JsonElement heights)
                    GeneratorHeights = Array.AsReadOnly(heights.EnumerateArray().Select(Real).ToArray());
                if (GeneratorHeights != null && MordellWeilGenerators != null && GeneratorHeights.Count != MordellWeilGenerators.Count)
                    throw new FormatException("LMFDB: generator and height counts differ.");
                if (Optional(m, "rank_bounds") is JsonElement bounds)
                {
                    if (bounds.GetArrayLength() != 2) throw new FormatException("LMFDB: invalid rank bounds.");
                    RankBounds = new LmfdbRankBounds(checked((int)Integer(bounds[0])), checked((int)Integer(bounds[1])));
                }
                TamagawaProduct = Optional(m, "tamagawa_product") is JsonElement cp ? Integer(cp) : (BigInteger?)null;
                Regulator = OptionalReal(m, "regulator") ?? Regulator;
                RealPeriod = OptionalReal(m, "real_period");
                PeriodArea = OptionalReal(m, "area");
            }
            var local = Table(root, "ec_localdata");
            if (local.HasValue)
            {
                var values = new List<LocalReductionData>();
                foreach (var l in local.Value.EnumerateArray())
                {
                    CheckLabel(l);
                    var prime = Integer(l.GetProperty("prime"));
                    int n = checked((int)Integer(l.GetProperty("discriminant_valuation")));
                    int f = checked((int)Integer(l.GetProperty("conductor_valuation")));
                    int j = checked((int)Integer(l.GetProperty("j_denominator_valuation")));
                    int c = checked((int)Integer(l.GetProperty("tamagawa_number")));
                    int w = checked((int)Integer(l.GetProperty("root_number")));
                    int reduction = checked((int)Integer(l.GetProperty("reduction_type")));
                    if (prime < 2 || n < 0 || f < 0 || j < 0 || c < 1 || (w != 1 && w != -1) || Math.Abs(reduction) > 1)
                        throw new FormatException("LMFDB: invalid local data.");
                    var type = n == 0 ? ReductionType.Good : reduction == 1 ? ReductionType.SplitMultiplicative :
                        reduction == -1 ? ReductionType.NonSplitMultiplicative : ReductionType.Additive;
                    values.Add(new LocalReductionData(prime, n, f, j, Kodaira(checked((int)Integer(l.GetProperty("kodaira_symbol")))), type, c, w));
                }
                values.Sort((a, b) => a.Prime.CompareTo(b.Prime));
                if (values.Select(x => x.Prime).Distinct().Count() != values.Count) throw new FormatException("LMFDB: duplicate local prime.");
                if (Optional(row, "bad_primes") is JsonElement bad && !bad.EnumerateArray().Select(Integer).OrderBy(x => x).SequenceEqual(values.Select(x => x.Prime)))
                    throw new FormatException("LMFDB: local data does not cover the recorded bad primes.");
                if (TamagawaProduct.HasValue && values.Aggregate(BigInteger.One, (a, x) => a * x.TamagawaNumber) != TamagawaProduct.Value)
                    throw new FormatException("LMFDB: inconsistent Tamagawa product.");
                LocalData = values.AsReadOnly();
            }
            ReadExtendedData(root, row);
        }

        private void CheckLabel(JsonElement row)
        {
            if (row.GetProperty("lmfdb_label").GetString() != Label) throw new FormatException("LMFDB: mismatched record labels.");
        }
        private IReadOnlyList<EllipticCurvePoint> Points(JsonElement row, string field)
        {
            if (!(Optional(row, field) is JsonElement array)) return null;
            var points = new List<EllipticCurvePoint>();
            foreach (var p in array.EnumerateArray())
            {
                if (p.GetArrayLength() != 3) throw new FormatException("LMFDB: expected weighted projective coordinates.");
                var z = Integer(p[2]);
                if (z.IsZero) throw new FormatException("LMFDB: generator at infinity.");
                var point = new EllipticCurvePoint(new BigRational(Integer(p[0]), z * z), new BigRational(Integer(p[1]), z * z * z));
                if (!GlobalMinimalModel.IsOnCurve(point)) throw new FormatException("LMFDB: generator is not on the stored model.");
                points.Add(point);
            }
            return points.AsReadOnly();
        }
        private static JsonElement? Table(JsonElement root, string name)
        {
            var tables = root.GetProperty("tables"); var data = root.GetProperty("data");
            if (tables.GetArrayLength() != data.GetArrayLength()) throw new FormatException("LMFDB: inconsistent table snapshot.");
            JsonElement? found = null;
            for (int i = 0; i < tables.GetArrayLength(); i++) if (tables[i].GetString() == name)
            {
                if (found.HasValue || data[i].ValueKind != JsonValueKind.Array) throw new FormatException("LMFDB: invalid table snapshot.");
                found = data[i];
                if (root.TryGetProperty("totals", out var totals) && Integer(totals[i]) != data[i].GetArrayLength())
                    throw new FormatException("LMFDB: incomplete " + name + " table.");
            }
            return found;
        }
        private static JsonElement? Optional(JsonElement row, string name)
            => row.TryGetProperty(name, out var value) && value.ValueKind != JsonValueKind.Null ? value : (JsonElement?)null;
        private static LmfdbRealValue OptionalReal(JsonElement row, string name) => Optional(row, name) is JsonElement x ? Real(x) : null;
        private static LmfdbRealValue Real(JsonElement x)
        {
            if (x.ValueKind == JsonValueKind.Object)
                return new LmfdbRealValue(x.GetProperty("data").GetString(), Optional(x, "prec") is JsonElement p ? checked((int)Integer(p)) : (int?)null);
            return new LmfdbRealValue(x.ValueKind == JsonValueKind.String ? x.GetString() : x.GetRawText(), null);
        }
        private static BigInteger Integer(JsonElement x) => BigInteger.Parse(x.ValueKind == JsonValueKind.String ? x.GetString() : x.GetRawText(), CultureInfo.InvariantCulture);
        private static EllipticCurveQ Model(JsonElement x)
        {
            if (x.ValueKind == JsonValueKind.String)
            {
                using var doc = JsonDocument.Parse(x.GetString());
                return Model(doc.RootElement);
            }
            if (x.GetArrayLength() != 5) throw new FormatException("LMFDB: expected five a-invariants.");
            return new EllipticCurveQ(new BigRational(Integer(x[0])), new BigRational(Integer(x[1])), new BigRational(Integer(x[2])), new BigRational(Integer(x[3])), new BigRational(Integer(x[4])));
        }
        private static string Kodaira(int code)
        {
            if (code == 0 || code == int.MinValue) throw new FormatException("LMFDB: invalid Kodaira code.");
            var n = Math.Abs(code);
            string symbol = n == 1 ? "I0" : n == 2 ? "II" : n == 3 ? "III" : n == 4 ? "IV" : "I" + (n - 4).ToString(CultureInfo.InvariantCulture);
            return symbol + (code < 0 ? "*" : "");
        }

        private static async Task<JsonElement> FetchSnapshotAsync(EllipticCurveQ curve, HttpClient client, CancellationToken token)
        {
            if (curve == null) throw new ArgumentNullException(nameof(curve));
            if (curve.Discriminant.IsZero) throw new ArgumentException("The curve is singular.", nameof(curve));
            token.ThrowIfCancellationRequested();
            var j = curve.JInvariant; string label = null;
            for (int offset = 0; offset <= 10000; offset += 100)
            {
                string url = Origin + "/api/ec_curvedata/?jinv=li" + j.Num.ToString(CultureInfo.InvariantCulture) + "," + j.Den.ToString(CultureInfo.InvariantCulture) +
                    "&_format=json&_fields=ainvs,lmfdb_label&_sort=lmfdb_label&_offset=" + offset.ToString(CultureInfo.InvariantCulture);
                var candidates = (await ReadJsonAsync(client, url, token).ConfigureAwait(false)).GetProperty("data");
                foreach (var candidate in candidates.EnumerateArray())
                {
                    token.ThrowIfCancellationRequested();
                    if (curve.TryGetIsomorphism(Model(candidate.GetProperty("ainvs")), out _))
                    { label = candidate.GetProperty("lmfdb_label").GetString(); break; }
                }
                if (label != null) break;
                if (candidates.GetArrayLength() < 100) throw new InvalidOperationException("LMFDB: no Q-isomorphic curve was found.");
                if (offset == 10000) throw new InvalidOperationException("LMFDB: candidate pagination limit reached before finding a Q-isomorphic curve.");
            }
            var snapshot = await ReadJsonAsync(client, Origin + "/EllipticCurve/Q/data/" + Uri.EscapeDataString(label) + "?_format=json&_limit=10000", token).ConfigureAwait(false);
            var parsed = new LmfdbEllipticCurve(snapshot);
            if (parsed.Label != label || !curve.TryGetIsomorphism(parsed.GlobalMinimalModel, out _)) throw new FormatException("LMFDB: fetched curve does not match the request.");
            return snapshot;
        }
        private static async Task<JsonElement> ReadJsonAsync(HttpClient client, string url, CancellationToken token)
        {
            using var request = new HttpRequestMessage(HttpMethod.Get, url);
            request.Headers.UserAgent.ParseAdd("EllipticCurves/" + typeof(LmfdbEllipticCurve).Assembly.GetName().Version.ToString(3));
            request.Headers.Accept.ParseAdd("application/json");
            using var response = await client.SendAsync(request, HttpCompletionOption.ResponseContentRead, token).ConfigureAwait(false);
            response.EnsureSuccessStatusCode();
            var json = await response.Content.ReadAsStringAsync().ConfigureAwait(false);
            token.ThrowIfCancellationRequested();
            if (json.IndexOf("captcha", StringComparison.OrdinalIgnoreCase) >= 0) throw new InvalidOperationException("LMFDB returned a CAPTCHA challenge.");
            using var document = JsonDocument.Parse(json);
            return document.RootElement.Clone();
        }
    }
}

using System;
using System.Globalization;
using System.Net.Http;
using System.Numerics;
using System.Text;
using System.Text.Json;

namespace EllipticCurves
{
    /// <summary>
    /// Thin helper that fetches basic invariants for an elliptic curve over ℚ from the LMFDB API,
    /// matching by the exact rational j-invariant and verifying Q–isomorphism via (c4,c6,Δ).
    /// 
    /// Design:
    ///  • Single network call to <c>/api/ec_curvedata</c> (no secondary requests).
    ///  • Minimal local cache of the selected record (a-invariants, conductor, rank, etc.).
    ///  • No retries/backoff: the caller controls lifetime by constructing a new instance.
    /// 
    /// Safety:
    ///  • Throws <see cref="InvalidOperationException"/> when the record is not initialized
    ///    or if no Q–isomorphic curve is found for the given j-invariant.
    /// </summary>
    public sealed partial class EllipticCurveLMFDB
    {
        /// <summary>
        /// Construct and immediately fetch/cached the matching LMFDB record
        /// for <paramref name="ellipticCurve"/> (matching by j–invariant and Q–isomorphism).
        /// </summary>
        public EllipticCurveLMFDB(EllipticCurveQ ellipticCurve) => GetLmfdbRecord(ellipticCurve);

        #region Properties

        /// <summary>
        /// Algebraic rank (Mordell–Weil rank) from LMFDB curvedata.
        /// </summary>
        public int Rank
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();
                return _lmfdbCache.Value.rank;
            }
        }

        /// <summary>
        /// Analytic rank (ord_{s=1} L(E,s)) when present in LMFDB; otherwise <c>null</c>.
        /// </summary>
        public int? AnalyticRank
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();
                return _lmfdbCache.Value.analyticRank;
            }
        }

        /// <summary>
        /// Conductor N of the (minimal) elliptic curve over ℚ.
        /// </summary>
        public BigInteger Conductor
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();
                return _lmfdbCache.Value.conductor;
            }
        }

        /// <summary>
        /// LMFDB label (e.g. <c>"48.a3"</c>) of the matched minimal model.
        /// </summary>
        public string Label
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();
                return _lmfdbCache.Value.label;
            }
        }

        /// <summary>
        /// Canonical human-facing LMFDB URL of the curve's page (or null if label is missing).
        /// </summary>
        public string Url
        {
            get
            {
                if (string.IsNullOrEmpty(Label))
                {
                    return null;
                }
                return $"https://www.lmfdb.org/EllipticCurve/Q/{Label}/";
            }
        }

        /// <summary>
        /// Torsion structure in a compact textual form, e.g. <c>"Z/2Z"</c> or <c>"Z/2Z x Z/4Z"</c>.
        /// </summary>
        public string TorsionStructure
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();
                return _lmfdbCache.Value.torsionStructure;
            }
        }

        /// <summary>
        /// Minimal (global) integral Weierstrass model reconstructed from the a-invariants in LMFDB.
        /// </summary>
        public EllipticCurveQ GlobalMinimalModel
        {
            get
            {
                if (!_lmfdbCache.HasValue)
                    throw new InvalidOperationException();

                var a = _lmfdbCache.Value.ainvs;
                return new EllipticCurveQ(
                    new BigRational(a[0]),
                    new BigRational(a[1]),
                    new BigRational(a[2]),
                    new BigRational(a[3]),
                    new BigRational(a[4]));
            }
        }

        #endregion

        #region Private methods

        // --- Local cache of a single curvedata record ---
        // Stored fields:
        //  • ainvs: minimal integral a-invariants [a1,a2,a3,a4,a6]
        //  • conductor: N
        //  • rank: algebraic rank
        //  • analyticRank: analytic rank if present (nullable)
        //  • label: LMFDB label (used for links)
        //  • torsionStructure: textual "Z/nZ" or "Z/aZ x Z/bZ" form
        private (BigInteger[] ainvs,
            BigInteger conductor,
            int rank,
            int? analyticRank,
            string label,
            string torsionStructure)? _lmfdbCache;

        /// <summary>
        /// Fetch a list of candidates by exact rational j-invariant, then select the unique
        /// curve Q–isomorphic to <paramref name="ellipticCurve"/> by verifying (c4,c6,Δ) scaling.
        /// Caches the selected record; subsequent property reads are served from the cache.
        /// </summary>
        private void GetLmfdbRecord(EllipticCurveQ ellipticCurve)
        {
            if (_lmfdbCache.HasValue) return;

            // 1) Build the API URL: jinv is a numeric[] [num,den], queried as "li<num>,<den>"
            var jinv = ellipticCurve.JInvariant; // BigRational (exact)
            string jnum = jinv.Num.ToString(CultureInfo.InvariantCulture);
            string jden = jinv.Den.ToString(CultureInfo.InvariantCulture);
            string url =
                $"https://www.lmfdb.org/api/ec_curvedata/?jinv=li{jnum},{jden}" +
                "&_format=json" +
                "&_fields=ainvs,conductor,rank,analytic_rank,iso_nlabel,torsion_structure,lmfdb_label";

            // 2) Perform HTTP GET (single-shot; the caller controls object lifetime)
            using var httpClient = new HttpClient();
            httpClient.Timeout = new TimeSpan(0, 0, 10);
            // Use a common browser UA to avoid any strict filters in frontends
            httpClient.DefaultRequestHeaders.Add(
                "User-Agent",
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0.0.0 Safari/537.36");

            using var response = httpClient.GetAsync(url).GetAwaiter().GetResult();
            var json = response.Content.ReadAsStringAsync().GetAwaiter().GetResult();

            if (json.Contains("captcha"))
                throw new InvalidOperationException("LMFDB: responded with a CAPTCHA challenge; wait and try again later or reduce request rate.");

            using var doc = JsonDocument.Parse(json);

            if (!doc.RootElement.TryGetProperty("data", out var data) || data.ValueKind != JsonValueKind.Array)
                throw new InvalidOperationException("LMFDB: unexpected response format.");

            if (data.GetArrayLength() == 0)
                throw new InvalidOperationException("LMFDB: no Q-isomorphic curve found for this j-invariant.");

            // Local invariants for Q–isomorphism validation
            var c4E = ellipticCurve.C4;
            var c6E = ellipticCurve.C6;
            var dE = ellipticCurve.Discriminant;

            // 3) Iterate candidates and select the one Q–isomorphic to the input curve
            for (int i = 0; i < data.GetArrayLength(); i++)
            {
                var row = data[i];
                var ainvs = ParseAinvs(row.GetProperty("ainvs"));
                var a1 = ainvs[0]; var a2 = ainvs[1]; 
                var a3 = ainvs[2]; var a4 = ainvs[3]; 
                var a6 = ainvs[4];

                var (c4C, c6C, dC) = InternalMath.InvariantsIntFromAinvs(a1, a2, a3, a4, a6);

                if (InternalMath.IsQIsomorphic(c4E, c6E, dE, c4C, c6C, dC, out _))
                {
                    var conductor = ReadBigInteger(row.GetProperty("conductor"));
                    var rank = row.GetProperty("rank").GetInt32();
                    var lmfdb_label = row.GetProperty("lmfdb_label").GetString();
                    var torsionStruct = FormatTorsionStructure(row.GetProperty("torsion_structure"));

                    // analytic_rank may be absent or null for some entries
                    int? analyticRank = null;
                    if (row.TryGetProperty("analytic_rank", out var arEl) && arEl.ValueKind == JsonValueKind.Number)
                        analyticRank = arEl.GetInt32();

                    _lmfdbCache = (ainvs, conductor, rank, analyticRank, lmfdb_label, torsionStruct);
                    return;
                }
            }

            throw new InvalidOperationException("LMFDB: unexpected error occured.");
        }

        /// <summary>
        /// Accepts either an integer array (e.g. [2,4]) or a ready string (e.g. "Z/2Z x Z/4Z")
        /// and returns a normalized textual form "Z/nZ" or "Z/aZ x Z/bZ".
        /// </summary>
        private static string FormatTorsionStructure(JsonElement el)
        {
            if (el.ValueKind == JsonValueKind.Array)
            {
                int n = el.GetArrayLength();
                if (n == 0) return "Z/1Z";
                var sb = new StringBuilder();

                for (int i = 0; i < n; i++)
                {
                    if (i > 0) sb.Append(" x ");
                    int k;
                    if (el[i].ValueKind == JsonValueKind.Number) k = el[i].GetInt32();
                    else if (el[i].ValueKind == JsonValueKind.String) k = int.Parse(el[i].GetString(), CultureInfo.InvariantCulture);
                    else k = 1;
                    sb.Append("Z/").Append(k).Append("Z");
                }

                return sb.ToString();
            }
            if (el.ValueKind == JsonValueKind.String)
            {
                var s = el.GetString() ?? string.Empty;
                return string.IsNullOrWhiteSpace(s) ? "Z/1Z" : s;
            }
            return "Z/1Z";
        }

        // ---- Helpers (parsing + invariants + Q-isomorphism check) ----

        /// <summary>
        /// Parse a-invariants from LMFDB payload. Supports either a JSON array of numbers
        /// or a string form "[a1,a2,a3,a4,a6]".
        /// </summary>
        private static BigInteger[] ParseAinvs(JsonElement el)
        {
            // Common case: an array of numbers
            if (el.ValueKind == JsonValueKind.Array)
            {
                var a = new BigInteger[5];
                for (int i = 0; i < 5; i++)
                    a[i] = ReadBigInteger(el[i]);
                return a;
            }

            // Fallback: legacy string representation "[a1,a2,a3,a4,a6]"
            if (el.ValueKind == JsonValueKind.String)
            {
                var s = el.GetString() ?? throw new FormatException("LMFDB: ainvs is null string.");
                s = s.Trim();
                if (s.StartsWith("[")) s = s.Substring(1);
                if (s.EndsWith("]")) s = s.Substring(0, s.Length - 1);

                var rawParts = s.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (rawParts.Length != 5) throw new FormatException("LMFDB: bad ainvs format.");

                var a = new BigInteger[5];
                for (int i = 0; i < 5; i++)
                    a[i] = BigInteger.Parse(rawParts[i].Trim(), CultureInfo.InvariantCulture);
                return a;
            }

            throw new FormatException("LMFDB: unexpected ainvs type");
        }

        /// <summary>
        /// Parse a BigInteger from either a numeric JSON token or a string token.
        /// </summary>
        private static BigInteger ReadBigInteger(JsonElement el)
        {
            return el.ValueKind switch
            {
                JsonValueKind.Number => BigInteger.Parse(el.GetRawText(), CultureInfo.InvariantCulture),
                JsonValueKind.String => BigInteger.Parse(el.GetString()!, CultureInfo.InvariantCulture),
                _ => throw new FormatException("LMFDB: bad integer field.")
            };
        }

        #endregion
    }
}

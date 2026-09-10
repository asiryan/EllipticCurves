using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text.Json;

namespace EllipticCurves
{
    public sealed partial class LmfdbEllipticCurve
    {
        /// <summary>Stored Cremona label, when present.</summary>
        public string CremonaLabel { get; private set; }
        /// <summary>LMFDB isogeny-class label, when present.</summary>
        public string IsogenyClassLabel { get; private set; }
        /// <summary>Stored size of the rational isogeny class.</summary>
        public int? IsogenyClassSize { get; private set; }
        /// <summary>Stored geometric CM discriminant; zero means non-CM, null means missing.</summary>
        public int? CmDiscriminant { get; private set; }
        /// <summary>Stored degrees of cyclic rational isogenies from this curve, including one.</summary>
        public IReadOnlyList<int> IsogenyDegrees { get; private set; }
        /// <summary>Stored matrix of minimal isogeny degrees in the database's class order; not computed maps.</summary>
        public IReadOnlyList<IReadOnlyList<int>> IsogenyMatrix { get; private set; }
        /// <summary>Stored Fourier coefficients indexed by n, including the database's a_0=0 sentinel.</summary>
        public IReadOnlyList<BigInteger> FourierCoefficients { get; private set; }
        /// <summary>Stored a_p for consecutive primes 2,3,5,... in that order.</summary>
        public IReadOnlyList<BigInteger> PrimeFourierCoefficients { get; private set; }
        /// <summary>Stored modular degree.</summary>
        public BigInteger? ModularDegree { get; private set; }
        /// <summary>Stored Manin constant.</summary>
        public BigInteger? ManinConstant { get; private set; }
        /// <summary>Stored order of rational torsion.</summary>
        public BigInteger? TorsionOrder { get; private set; }
        /// <summary>Stored Faltings height approximation.</summary>
        public LmfdbRealValue FaltingsHeight { get; private set; }
        /// <summary>Stored stable Faltings height approximation.</summary>
        public LmfdbRealValue StableFaltingsHeight { get; private set; }
        /// <summary>Stored rounded analytic order of Sha predicted by BSD; not a native proof of the group order.</summary>
        public BigInteger? BsdShaOrder { get; private set; }
        /// <summary>Stored unrounded analytic Sha approximation; not a certificate.</summary>
        public LmfdbRealValue AnalyticShaOrder { get; private set; }
        /// <summary>Stored leading Taylor coefficient L^(r)(E,1)/r!, not the unscaled derivative.</summary>
        public LmfdbRealValue LeadingLValue { get; private set; }
        /// <summary>Stored integral-point x-coordinates on GlobalMinimalModel; no native completeness claim.</summary>
        public IReadOnlyList<BigInteger> IntegralPointXCoordinates { get; private set; }

        private void ReadExtendedData(JsonElement root, JsonElement row)
        {
            BigInteger? PositiveInteger(JsonElement record, string name)
            {
                if (!(Optional(record, name) is JsonElement value)) return null;
                var n = Integer(value);
                if (n <= 0) throw new FormatException("LMFDB: invalid positive integer in " + name + ".");
                return n;
            }
            IReadOnlyList<BigInteger> Integers(JsonElement record, string name)
            {
                if (!(Optional(record, name) is JsonElement value)) return null;
                if (value.ValueKind != JsonValueKind.Array) throw new FormatException("LMFDB: expected an array in " + name + ".");
                return Array.AsReadOnly(value.EnumerateArray().Select(Integer).ToArray());
            }
            CremonaLabel = Optional(row, "Clabel")?.GetString();
            IsogenyClassLabel = Optional(row, "lmfdb_iso")?.GetString();
            var expectedClassLabel = Label.TrimEnd("0123456789".ToCharArray());
            if (IsogenyClassLabel != null && IsogenyClassLabel != expectedClassLabel) throw new FormatException("LMFDB: mismatched isogeny-class label.");
            var size = PositiveInteger(row, "class_size"); IsogenyClassSize = size.HasValue ? checked((int)size.Value) : (int?)null;
            CmDiscriminant = Optional(row, "cm") is JsonElement cm ? checked((int)Integer(cm)) : (int?)null;
            if (CmDiscriminant > 0) throw new FormatException("LMFDB: invalid CM discriminant.");
            if (Integers(row, "isogeny_degrees") is IReadOnlyList<BigInteger> degrees)
            {
                if (degrees.Any(d => d <= 0)) throw new FormatException("LMFDB: invalid isogeny degree.");
                IsogenyDegrees = Array.AsReadOnly(degrees.Select(d => checked((int)d)).ToArray());
            }
            ModularDegree = PositiveInteger(row, "degree"); ManinConstant = PositiveInteger(row, "manin_constant");
            TorsionOrder = PositiveInteger(row, "torsion"); BsdShaOrder = PositiveInteger(row, "sha");
            FaltingsHeight = OptionalReal(row, "faltings_height"); StableFaltingsHeight = OptionalReal(row, "stable_faltings_height");

            var classes = Table(root, "ec_classdata");
            if (classes.HasValue && classes.Value.GetArrayLength() > 1) throw new FormatException("LMFDB: multiple isogeny-class records.");
            if (classes.HasValue && classes.Value.GetArrayLength() == 1)
            {
                var record = classes.Value[0]; var label = record.GetProperty("lmfdb_iso").GetString();
                if (label != expectedClassLabel) throw new FormatException("LMFDB: mismatched isogeny-class label.");
                IsogenyClassLabel = label;
                var classSize = PositiveInteger(record, "class_size");
                if (classSize.HasValue)
                {
                    if (IsogenyClassSize.HasValue && classSize != IsogenyClassSize.Value) throw new FormatException("LMFDB: inconsistent isogeny-class size.");
                    IsogenyClassSize = checked((int)classSize.Value);
                }
                FourierCoefficients = Integers(record, "anlist"); PrimeFourierCoefficients = Integers(record, "aplist");
                if (FourierCoefficients != null && ((FourierCoefficients.Count > 0 && FourierCoefficients[0] != 0) ||
                    (FourierCoefficients.Count > 1 && FourierCoefficients[1] != 1))) throw new FormatException("LMFDB: invalid coefficient indexing.");
                if (Optional(record, "isogeny_matrix") is JsonElement matrix)
                {
                    var rows = matrix.EnumerateArray().Select(r => r.EnumerateArray().Select(x => checked((int)Integer(x))).ToArray()).ToArray();
                    int n = rows.Length;
                    if (rows.Any(r => r.Length != n) || (IsogenyClassSize.HasValue && n != IsogenyClassSize.Value))
                        throw new FormatException("LMFDB: invalid isogeny-matrix dimensions.");
                    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                        if (rows[i][j] < 1 || rows[i][j] != rows[j][i] || (i == j && rows[i][j] != 1))
                            throw new FormatException("LMFDB: invalid isogeny matrix.");
                    IsogenyMatrix = Array.AsReadOnly(rows.Select(r => (IReadOnlyList<int>)Array.AsReadOnly(r)).ToArray());
                }
            }
            var mw = Table(root, "ec_mwbsd");
            if (mw.HasValue && mw.Value.GetArrayLength() == 1)
            {
                var record = mw.Value[0];
                AnalyticShaOrder = OptionalReal(record, "sha_an"); LeadingLValue = OptionalReal(record, "special_value");
                IntegralPointXCoordinates = Integers(record, "xcoord_integral_points");
            }
        }
    }
}

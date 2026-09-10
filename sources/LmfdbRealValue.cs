using System;
using System.Globalization;
using System.Numerics;
using System.Text.RegularExpressions;

namespace EllipticCurves
{
    /// <summary>A stored LMFDB decimal approximation, not a certified error interval.</summary>
    public sealed class LmfdbRealValue
    {
        private static readonly Regex DecimalPattern = new Regex(@"\A[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\z", RegexOptions.CultureInvariant);

        /// <summary>The decimal text as supplied by the database.</summary>
        public string DecimalValue { get; }
        /// <summary>Stored precision metadata, when supplied; this is not an error bound.</summary>
        public int? StoredPrecisionBits { get; }
        /// <summary>A double approximation for display.</summary>
        public double Approximation => double.Parse(DecimalValue, CultureInfo.InvariantCulture);
        internal LmfdbRealValue(string value, int? precision)
        {
            DecimalValue = value; StoredPrecisionBits = precision;
            AsRational();
        }
        /// <summary>Converts the stored decimal to an exact rational, without asserting mathematical accuracy.</summary>
        public BigRational AsRational()
        {
            var value = DecimalValue?.Trim() ?? string.Empty;
            if (!DecimalPattern.IsMatch(value)) throw new FormatException("LMFDB: invalid decimal.");
            var parts = value.Split('e', 'E');
            int exponent = 0;
            if (parts.Length == 2 && !int.TryParse(parts[1], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out exponent))
                throw new FormatException("LMFDB: decimal exponent exceeds supported range.");
            if (exponent < -10000 || exponent > 10000) throw new FormatException("LMFDB: decimal exponent exceeds supported range.");
            var mantissa = parts[0]; int dot = mantissa.IndexOf('.');
            if (dot >= 0) { exponent -= mantissa.Length - dot - 1; mantissa = mantissa.Remove(dot, 1); }
            var numerator = BigInteger.Parse(mantissa, CultureInfo.InvariantCulture);
            return exponent >= 0 ? new BigRational(numerator * BigInteger.Pow(10, exponent)) : new BigRational(numerator, BigInteger.Pow(10, -exponent));
        }
        /// <inheritdoc/>
        public override string ToString() => DecimalValue;
    }
}

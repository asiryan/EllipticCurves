using System;

namespace EllipticCurves
{
    /// <summary>Acceptance limits for exact rational text parsing, snapshotted at the start of each call.
    /// Early size estimates and exact integer-size checks enforce these limits; they are not elapsed-time or peak-allocation bounds.</summary>
    public sealed class RationalParsingOptions
    {
        /// <summary>Maximum input length, including whitespace. Must be positive; defaults to 20,000 characters.</summary>
        public int MaxTextLength { get; set; } = 20000;
        /// <summary>Maximum absolute decimal exponent in either fraction component. Must be nonnegative; defaults to 4096.</summary>
        public int MaxExponent { get; set; } = 4096;
        /// <summary>Maximum bit length of the final normalized numerator and denominator. Must be positive; defaults to 32,768 bits.</summary>
        public int MaxValueBits { get; set; } = 32768;
        /// <summary>Maximum accepted bit length of component integers, decimal powers and final fraction products, before normalization.
        /// Must be positive; defaults to 131,072 bits. This separate budget permits large components that cancel in a fraction.</summary>
        public int MaxIntermediateBits { get; set; } = 131072;

        internal RationalParsingOptions Snapshot()
        {
            var limits = new RationalParsingOptions
            {
                MaxTextLength = MaxTextLength, MaxExponent = MaxExponent,
                MaxValueBits = MaxValueBits, MaxIntermediateBits = MaxIntermediateBits
            };
            if (limits.MaxTextLength < 1) throw new ArgumentOutOfRangeException(nameof(MaxTextLength));
            if (limits.MaxExponent < 0) throw new ArgumentOutOfRangeException(nameof(MaxExponent));
            if (limits.MaxValueBits < 1) throw new ArgumentOutOfRangeException(nameof(MaxValueBits));
            if (limits.MaxIntermediateBits < 1) throw new ArgumentOutOfRangeException(nameof(MaxIntermediateBits));
            return limits;
        }
    }
}

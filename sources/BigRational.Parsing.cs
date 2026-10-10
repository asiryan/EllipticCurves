using System;
using System.Globalization;
using System.Numerics;

namespace EllipticCurves
{
    public readonly partial struct BigRational
    {
        /// <summary>Parse exact invariant text: a signed ASCII integer or decimal, optional e/E exponent,
        /// or two such numbers separated by /. Leading/trailing whitespace and Unicode minus are accepted.
        /// No floating-point conversion is used. Default resource limits apply.</summary>
        public static BigRational Parse(string text) => Parse(text, null);

        /// <summary>Parse exact invariant text with the supplied limits; null options select defaults.
        /// Null text throws ArgumentNullException. Malformed text, a zero denominator or an exceeded limit throws FormatException.
        /// Invalid limits throw ArgumentOutOfRangeException.</summary>
        public static BigRational Parse(string text, RationalParsingOptions options)
        {
            var limits = (options ?? new RationalParsingOptions()).Snapshot();
            if (text == null) throw new ArgumentNullException(nameof(text));
            if (!TryParseCore(text, limits, out var value))
                throw new FormatException("Expected an exact rational number within the configured parsing limits.");
            return value;
        }

        /// <summary>Try to parse exact invariant text with default resource limits.
        /// Null, malformed or excessive input and a zero denominator return false with value set to zero.</summary>
        public static bool TryParse(string text, out BigRational value) => TryParse(text, null, out value);

        /// <summary>Try to parse exact invariant text with the supplied limits; null options select defaults.
        /// A failed parse clears value to zero. Invalid limits throw ArgumentOutOfRangeException.</summary>
        public static bool TryParse(string text, RationalParsingOptions options, out BigRational value)
            => TryParseCore(text, (options ?? new RationalParsingOptions()).Snapshot(), out value);

        private static bool TryParseCore(string text, RationalParsingOptions limits, out BigRational value)
        {
            value = Zero;
            if (text == null || text.Length > limits.MaxTextLength) return false;
            int slash = -1;
            for (int i = 0; i < text.Length; i++)
                if (text[i] == '/')
                {
                    if (slash >= 0) return false;
                    slash = i;
                }
            if (!TryParseNumber(text, 0, slash < 0 ? text.Length : slash, limits, out var numerator)) return false;
            var candidate = numerator;
            if (slash >= 0)
            {
                if (!TryParseNumber(text, slash + 1, text.Length, limits, out var denominator) || denominator.IsZero) return false;
                if (numerator.IsZero) return true;
                // Cancel common factors before multiplying; N/N can be small even
                // when each component exceeds the final-value budget.
                var first = BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator.Num), BigInteger.Abs(denominator.Num));
                var second = BigInteger.GreatestCommonDivisor(numerator.Den, denominator.Den);
                var a = numerator.Num / first; var b = denominator.Den / second;
                var c = numerator.Den / second; var d = denominator.Num / first;
                int productLimit = Math.Min(limits.MaxValueBits, limits.MaxIntermediateBits);
                if (!TryParseProduct(a, b, productLimit, out var num) || !TryParseProduct(c, d, productLimit, out var den)) return false;
                candidate = new BigRational(num, den);
            }
            if (!ParseWithinBits(candidate.Num, limits.MaxValueBits) || !ParseWithinBits(candidate.Den, limits.MaxValueBits)) return false;
            value = candidate;
            return true;
        }

        private static bool TryParseNumber(string text, int start, int end, RationalParsingOptions limits, out BigRational value)
        {
            value = Zero;
            while (start < end && char.IsWhiteSpace(text[start])) start++;
            while (end > start && char.IsWhiteSpace(text[end - 1])) end--;
            if (start == end) return false;
            bool negative = false;
            if (text[start] == '+' || ParseIsMinus(text[start])) { negative = ParseIsMinus(text[start]); start++; }
            int integerStart = start;
            while (start < end && ParseIsDigit(text[start])) start++;
            int integerEnd = start, fractionStart = start;
            if (start < end && text[start] == '.')
            {
                fractionStart = ++start;
                while (start < end && ParseIsDigit(text[start])) start++;
            }
            int fractionEnd = start;
            int integerDigits = integerEnd - integerStart, fractionDigits = fractionEnd - fractionStart;
            int digitCount = integerDigits + fractionDigits;
            if (digitCount == 0) return false;
            int exponent = 0;
            if (start < end && (text[start] == 'e' || text[start] == 'E'))
            {
                start++;
                bool negativeExponent = false;
                if (start < end && (text[start] == '+' || ParseIsMinus(text[start])))
                { negativeExponent = ParseIsMinus(text[start]); start++; }
                int exponentStart = start;
                while (start < end && ParseIsDigit(text[start]))
                {
                    int digit = text[start++] - '0';
                    if (exponent > limits.MaxExponent / 10 ||
                        (exponent == limits.MaxExponent / 10 && digit > limits.MaxExponent % 10)) return false;
                    exponent = exponent * 10 + digit;
                }
                if (start == exponentStart) return false;
                if (negativeExponent) exponent = -exponent;
            }
            if (start != end) return false;

            char Digit(int index) => text[index < integerDigits ? integerStart + index : fractionStart + index - integerDigits];
            int first = 0, last = digitCount - 1;
            while (first <= last && Digit(first) == '0') first++;
            if (first > last) return true;
            while (Digit(last) == '0') last--;
            long scale = (long)fractionDigits - exponent - (digitCount - 1 - last);
            int significantDigits = last - first + 1;
            long power = Math.Abs(scale);
            // 10^k >= 2^(3k). These necessary bounds reject excessive powers
            // and digit strings before allocating BigIntegers, without logarithms.
            if (3L * (significantDigits - 1) >= limits.MaxIntermediateBits || 3L * power >= limits.MaxIntermediateBits) return false;
            var digits = new char[significantDigits];
            for (int i = 0; i < significantDigits; i++) digits[i] = Digit(first + i);
            if (!BigInteger.TryParse(new string(digits), NumberStyles.None, CultureInfo.InvariantCulture, out var numerator) ||
                !ParseWithinBits(numerator, limits.MaxIntermediateBits)) return false;
            var tenPower = BigInteger.Pow(10, (int)power);
            if (!ParseWithinBits(tenPower, limits.MaxIntermediateBits)) return false;
            BigInteger denominator = BigInteger.One;
            if (scale >= 0) denominator = tenPower;
            else if (!TryParseProduct(numerator, tenPower, limits.MaxIntermediateBits, out numerator)) return false;
            if (negative) numerator = -numerator;
            value = new BigRational(numerator, denominator);
            return true;
        }

        private static bool TryParseProduct(BigInteger a, BigInteger b, int maxBits, out BigInteger result)
        {
            result = BigInteger.Zero;
            if (a.IsZero || b.IsZero) return true;
            if ((long)RealArithmetic.BitLength(a) + RealArithmetic.BitLength(b) - 1 > maxBits) return false;
            var product = a * b;
            if (!ParseWithinBits(product, maxBits)) return false;
            result = product;
            return true;
        }

        private static bool ParseWithinBits(BigInteger value, int maxBits) => RealArithmetic.BitLength(value) <= maxBits;
        private static bool ParseIsDigit(char value) => value >= '0' && value <= '9';
        private static bool ParseIsMinus(char value) => value == '-' || value == '\u2212';
    }
}

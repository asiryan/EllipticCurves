using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>An immutable finite-field element in a specified polynomial presentation.
    /// A default-initialized value is invalid, not the zero of an unspecified field.</summary>
    public readonly struct FiniteFieldElement : IEquatable<FiniteFieldElement>
    {
        private readonly BigInteger[] coefficients;
        internal BigInteger[] Coordinates => coefficients ?? Array.Empty<BigInteger>();
        /// <summary>Defining field; null only for an invalid default value.</summary>
        public FiniteField Field { get; }
        /// <summary>Canonical ascending coefficients, with trailing zeros removed. Zero has an empty list.</summary>
        public IReadOnlyList<BigInteger> Coefficients => Array.AsReadOnly(Coordinates);
        /// <summary>Whether this is a valid zero element.</summary>
        public bool IsZero => Field != null && coefficients.Length == 0;
        /// <summary>Whether this is a valid multiplicative identity.</summary>
        public bool IsOne => Field != null && coefficients.Length == 1 && coefficients[0].IsOne;
        internal FiniteFieldElement(FiniteField field, BigInteger[] ownedCoefficients)
        { Field = field; coefficients = ownedCoefficients; }
        private FiniteField RequireField() => Field ?? throw new InvalidOperationException("A default field element is invalid.");
        /// <summary>Multiplicative inverse of a nonzero element.</summary>
        public FiniteFieldElement Inverse(CancellationToken cancellationToken = default) => RequireField().Inverse(this, cancellationToken);
        /// <summary>Signed integer power; a^0 is one, including 0^0.</summary>
        public FiniteFieldElement Pow(BigInteger exponent, CancellationToken cancellationToken = default) => RequireField().Pow(this, exponent, cancellationToken);
        /// <summary>Add elements in the same field presentation.</summary>
        public static FiniteFieldElement operator +(FiniteFieldElement a, FiniteFieldElement b) => a.RequireField().Add(a, b);
        /// <summary>Subtract elements in the same field presentation.</summary>
        public static FiniteFieldElement operator -(FiniteFieldElement a, FiniteFieldElement b) => a.RequireField().Subtract(a, b);
        /// <summary>Additive inverse.</summary>
        public static FiniteFieldElement operator -(FiniteFieldElement a) => a.RequireField().Negate(a);
        /// <summary>Multiply elements in the same field presentation.</summary>
        public static FiniteFieldElement operator *(FiniteFieldElement a, FiniteFieldElement b) => a.RequireField().Multiply(a, b);
        /// <summary>Divide by a nonzero element in the same field presentation.</summary>
        public static FiniteFieldElement operator /(FiniteFieldElement a, FiniteFieldElement b) => a.RequireField().Divide(a, b);
        /// <summary>Add a prime-subfield integer.</summary>
        public static FiniteFieldElement operator +(FiniteFieldElement a, BigInteger b) => a + a.RequireField().CreateElement(b);
        /// <summary>Add a prime-subfield integer.</summary>
        public static FiniteFieldElement operator +(BigInteger a, FiniteFieldElement b) => b + a;
        /// <summary>Subtract a prime-subfield integer.</summary>
        public static FiniteFieldElement operator -(FiniteFieldElement a, BigInteger b) => a - a.RequireField().CreateElement(b);
        /// <summary>Subtract an element from a prime-subfield integer.</summary>
        public static FiniteFieldElement operator -(BigInteger a, FiniteFieldElement b) => b.RequireField().CreateElement(a) - b;
        /// <summary>Multiply by a prime-subfield integer.</summary>
        public static FiniteFieldElement operator *(FiniteFieldElement a, BigInteger b) => a * a.RequireField().CreateElement(b);
        /// <summary>Multiply by a prime-subfield integer.</summary>
        public static FiniteFieldElement operator *(BigInteger a, FiniteFieldElement b) => b * a;
        /// <inheritdoc/>
        public bool Equals(FiniteFieldElement other) => Equals(Field, other.Field) && Coordinates.SequenceEqual(other.Coordinates);
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is FiniteFieldElement element && Equals(element);
        /// <inheritdoc/>
        public override int GetHashCode()
        { int hash = Field?.GetHashCode() ?? 0; foreach (var c in Coordinates) hash = unchecked(hash * 397 ^ c.GetHashCode()); return hash; }
        /// <summary>Equality includes the field presentation.</summary>
        public static bool operator ==(FiniteFieldElement a, FiniteFieldElement b) => a.Equals(b);
        /// <summary>Inequality includes the field presentation.</summary>
        public static bool operator !=(FiniteFieldElement a, FiniteFieldElement b) => !a.Equals(b);
        /// <inheritdoc/>
        public override string ToString() => Field == null ? "<invalid>" : $"[{string.Join(", ", Coordinates)}] in {Field}";
    }
}

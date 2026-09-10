using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    /// <summary>The finite field F_p[t]/(f), with proved prime p and proved irreducible f.
    /// Elements in different polynomial presentations are not implicitly identified.</summary>
    public sealed class FiniteField : IEquatable<FiniteField>
    {
        private readonly BigInteger[] modulus;
        private readonly int hashCode;
        /// <summary>Prime characteristic.</summary>
        public BigInteger Characteristic { get; }
        /// <summary>Extension degree over the prime field, at least one.</summary>
        public int Degree => modulus.Length - 1;
        /// <summary>Number of elements, p^Degree.</summary>
        public BigInteger Order { get; }
        /// <summary>Monic defining polynomial, with coefficients in ascending power order.</summary>
        public IReadOnlyList<BigInteger> Modulus { get; }
        /// <summary>Additive identity.</summary>
        public FiniteFieldElement Zero { get; }
        /// <summary>Multiplicative identity.</summary>
        public FiniteFieldElement One { get; }
        /// <summary>Residue class of t. It need not generate the multiplicative group.</summary>
        public FiniteFieldElement Generator { get; }

        /// <summary>Create a field from a polynomial in ascending coefficient order.
        /// Coefficients are reduced modulo p and the polynomial is made monic before the exact Rabin test.
        /// Primality and irreducibility computations are cancellable; there is no fixed time bound.</summary>
        public FiniteField(BigInteger prime, IReadOnlyList<BigInteger> definingPolynomial, CancellationToken cancellationToken = default)
        {
            if (definingPolynomial == null) throw new ArgumentNullException(nameof(definingPolynomial));
            cancellationToken.ThrowIfCancellationRequested();
            if (!NativeNumberTheory.IsPrime(prime, cancellationToken)) throw new ArgumentOutOfRangeException(nameof(prime));
            var f = FiniteFieldPolynomial.Normalize(definingPolynomial, prime, cancellationToken);
            if (f.Length < 2) throw new ArgumentException("The defining polynomial must have positive degree modulo p.", nameof(definingPolynomial));
            var inverse = BigInteger.ModPow(f[f.Length - 1], prime - 2, prime);
            for (int i = 0; i < f.Length; i++) f[i] = NativeNumberTheory.Mod(f[i] * inverse, prime);
            if (!FiniteFieldPolynomial.IsIrreducible(f, prime, cancellationToken))
                throw new ArgumentException("The defining polynomial is reducible over the prime field.", nameof(definingPolynomial));
            Characteristic = prime; modulus = f; Modulus = Array.AsReadOnly(f); Order = BigInteger.Pow(prime, Degree);
            int hash = prime.GetHashCode();
            foreach (var coefficient in f) hash = unchecked(hash * 397 ^ coefficient.GetHashCode());
            hashCode = hash;
            Zero = new FiniteFieldElement(this, Array.Empty<BigInteger>());
            One = new FiniteFieldElement(this, new BigInteger[] { 1 });
            Generator = CreateElement(0, 1);
        }

        /// <summary>Create an element from ascending polynomial coefficients, reducing modulo p and f.</summary>
        public FiniteFieldElement CreateElement(params BigInteger[] coefficients)
        {
            if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));
            return new FiniteFieldElement(this, FiniteFieldPolynomial.Remainder(
                FiniteFieldPolynomial.Normalize(coefficients, Characteristic, default), modulus, Characteristic, default));
        }

        internal void Validate(FiniteFieldElement a)
        {
            if (!Equals(a.Field)) throw new ArgumentException("The element must belong to this polynomial presentation of the field.", nameof(a));
        }
        /// <summary>Add two elements of this field.</summary>
        public FiniteFieldElement Add(FiniteFieldElement a, FiniteFieldElement b)
        { Validate(a); Validate(b); return new FiniteFieldElement(this, FiniteFieldPolynomial.Add(a.Coordinates, b.Coordinates, Characteristic)); }
        /// <summary>Subtract two elements of this field.</summary>
        public FiniteFieldElement Subtract(FiniteFieldElement a, FiniteFieldElement b)
        { Validate(a); Validate(b); return new FiniteFieldElement(this, FiniteFieldPolynomial.Add(a.Coordinates, b.Coordinates, Characteristic, true)); }
        /// <summary>Negate an element of this field.</summary>
        public FiniteFieldElement Negate(FiniteFieldElement a) => Subtract(Zero, a);
        /// <summary>Multiply two elements of this field.</summary>
        public FiniteFieldElement Multiply(FiniteFieldElement a, FiniteFieldElement b, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(a); Validate(b);
            return new FiniteFieldElement(this, FiniteFieldPolynomial.Multiply(a.Coordinates, b.Coordinates, modulus, Characteristic, cancellationToken));
        }
        /// <summary>Invert a nonzero field element. Zero throws DivideByZeroException.</summary>
        public FiniteFieldElement Inverse(FiniteFieldElement a, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(a);
            if (a.IsZero) throw new DivideByZeroException();
            return Pow(a, Order - 2, cancellationToken);
        }
        /// <summary>Divide by a nonzero element.</summary>
        public FiniteFieldElement Divide(FiniteFieldElement a, FiniteFieldElement b, CancellationToken cancellationToken = default)
        { cancellationToken.ThrowIfCancellationRequested(); Validate(a); Validate(b); return Multiply(a, Inverse(b, cancellationToken), cancellationToken); }
        /// <summary>Raise to a signed integer power. Negative exponents require nonzero input; a^0 is one, including 0^0.</summary>
        public FiniteFieldElement Pow(FiniteFieldElement a, BigInteger exponent, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested(); Validate(a);
            if (exponent.Sign < 0) { a = Inverse(a, cancellationToken); exponent = -exponent; }
            return new FiniteFieldElement(this, FiniteFieldPolynomial.Power(a.Coordinates, exponent, modulus, Characteristic, cancellationToken));
        }

        /// <summary>Enumerate every element in increasing base-p coefficient encoding, starting at zero.
        /// Throws before enumeration if the cardinality exceeds maxElements.</summary>
        public IEnumerable<FiniteFieldElement> Elements(long maxElements = 1000000, CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (maxElements < 0) throw new ArgumentOutOfRangeException(nameof(maxElements));
            if (Order > maxElements) throw new ArithmeticException("Field enumeration exceeds the element limit.");
            for (BigInteger index = 0; index < Order; index++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                var value = index; var coefficients = new List<BigInteger>();
                while (value > 0) { coefficients.Add(value % Characteristic); value /= Characteristic; }
                yield return new FiniteFieldElement(this, coefficients.ToArray());
            }
        }
        /// <inheritdoc/>
        public bool Equals(FiniteField other) => ReferenceEquals(this, other) ||
            (other != null && Characteristic == other.Characteristic && modulus.SequenceEqual(other.modulus));
        /// <inheritdoc/>
        public override bool Equals(object obj) => obj is FiniteField field && Equals(field);
        /// <inheritdoc/>
        public override int GetHashCode() => hashCode;
        /// <inheritdoc/>
        public override string ToString() => $"GF({Characteristic}^{Degree}), modulus=[{string.Join(", ", modulus)}]";
    }
}

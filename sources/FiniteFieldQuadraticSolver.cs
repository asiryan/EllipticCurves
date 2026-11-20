using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Exact quadratic solving for a fixed field presentation. Small odd fields use
    // a square table; larger ones use Tonelli-Shanks without an O(q) allocation.
    // In characteristic two, squaring and z -> z^2+z are F_2-linear maps.
    // Their echelon forms are prepared once, including in even extension degrees.
    internal sealed class FiniteFieldQuadraticSolver
    {
        internal const int MaximumTableOrder = 1048576;
        private readonly FiniteField field;
        private readonly int[] squareRoots;
        private readonly BinaryLinearMap squaring;
        private readonly BinaryLinearMap artinSchreier;
        private readonly BigInteger binaryModulus;
        private FiniteFieldElement nonSquarePower;
        private BigInteger oddPart;
        private int twoAdicOrder;

        internal FiniteFieldQuadraticSolver(FiniteField field, CancellationToken token)
        {
            this.field = field;
            token.ThrowIfCancellationRequested();
            if (field.Characteristic == 2)
            {
                for (int i = field.Modulus.Count - 1; i >= 0; i--)
                {
                    token.ThrowIfCancellationRequested();
                    binaryModulus = (binaryModulus << 1) + field.Modulus[i];
                }
                squaring = new BinaryLinearMap(field, false, token);
                artinSchreier = new BinaryLinearMap(field, true, token);
            }
            else if (field.Order <= MaximumTableOrder)
            {
                squareRoots = new int[(int)field.Order];
                for (int code = 0; code < squareRoots.Length; code++)
                {
                    token.ThrowIfCancellationRequested();
                    var value = Decode(code);
                    int square = Encode(field.Multiply(value, value, token));
                    // Zero means "not a square"; storing code+1 includes sqrt(0).
                    if (squareRoots[square] == 0) squareRoots[square] = code + 1;
                }
            }
            token.ThrowIfCancellationRequested();
        }

        internal bool TrySquareRoot(FiniteFieldElement value, CancellationToken token, out FiniteFieldElement root)
        {
            token.ThrowIfCancellationRequested();
            if (squaring != null) return squaring.TrySolve(value, token, out root);
            if (squareRoots != null)
            {
                int code = squareRoots[Encode(value)];
                root = code == 0 ? field.Zero : Decode(code - 1);
                return code != 0;
            }
            root = field.Zero;
            if (value.IsZero) return true;
            if (field.Order % 4 == 3)
            {
                root = field.Pow(value, (field.Order + 1) / 4, token);
                if (field.Multiply(root, root, token) == value) return true;
                root = field.Zero;
                return false;
            }
            if (!field.Pow(value, (field.Order - 1) / 2, token).IsOne) return false;
            PrepareTonelliShanks(token);
            root = field.Pow(value, (oddPart + 1) / 2, token);
            var remainder = field.Pow(value, oddPart, token);
            var c = nonSquarePower;
            int order = twoAdicOrder;
            while (!remainder.IsOne)
            {
                token.ThrowIfCancellationRequested();
                int i = 0; var power = remainder;
                while (!power.IsOne && i < order)
                {
                    power = field.Multiply(power, power, token);
                    i++;
                }
                if (i == order) throw new ArithmeticException("Finite-field square-root iteration failed.");
                var factor = field.Pow(c, BigInteger.One << (order - i - 1), token);
                root = field.Multiply(root, factor, token);
                c = field.Multiply(factor, factor, token);
                remainder = field.Multiply(remainder, c, token);
                order = i;
            }
            token.ThrowIfCancellationRequested();
            return true;
        }

        internal bool TryArtinSchreierRoot(FiniteFieldElement value, CancellationToken token, out FiniteFieldElement root)
            => artinSchreier.TrySolve(value, token, out root);

        internal bool IsSquare(FiniteFieldElement value, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (squaring != null || value.IsZero) return true;
            return squareRoots != null ? squareRoots[Encode(value)] != 0
                : field.Pow(value, (field.Order - 1) / 2, token).IsOne;
        }

        internal bool HasArtinSchreierRoot(FiniteFieldElement value, CancellationToken token)
            => artinSchreier.HasImage(value, token);

        // Polynomial extended Euclid over F_2, with polynomials packed into bits.
        // u = value*g1 and v = value*g2 modulo the defining polynomial throughout.
        // Irreducibility and value != 0 guarantee gcd(value, modulus) = 1.
        // Squaring the resulting inverse uses bit interleaving and exact reduction,
        // avoiding an exponentiation and all dense polynomial multiplications.
        internal FiniteFieldElement InverseSquare(FiniteFieldElement value, CancellationToken token)
        {
            token.ThrowIfCancellationRequested();
            if (field.Characteristic != 2) throw new InvalidOperationException("Binary polynomial inversion requires characteristic two.");
            field.Validate(value);
            if (value.IsZero) throw new DivideByZeroException();
            BigInteger u = 0;
            for (int i = value.Coordinates.Length - 1; i >= 0; i--)
            {
                token.ThrowIfCancellationRequested();
                u = (u << 1) + value.Coordinates[i];
            }
            var v = binaryModulus;
            BigInteger g1 = 1, g2 = 0;
            while (!u.IsOne)
            {
                token.ThrowIfCancellationRequested();
                if (u.IsZero) throw new ArithmeticException("A nonzero finite-field element must be invertible.");
                int shift = RealArithmetic.BitLength(u) - RealArithmetic.BitLength(v);
                if (shift < 0)
                {
                    var swap = u; u = v; v = swap;
                    swap = g1; g1 = g2; g2 = swap;
                    shift = -shift;
                }
                u ^= v << shift;
                g1 ^= g2 << shift;
            }
            BigInteger square = 0;
            for (int bit = 0; !g1.IsZero; bit++, g1 >>= 1)
            {
                token.ThrowIfCancellationRequested();
                if (!g1.IsEven) square |= BigInteger.One << (2 * bit);
            }
            int degree;
            while ((degree = RealArithmetic.BitLength(square) - 1) >= field.Degree)
            {
                token.ThrowIfCancellationRequested();
                square ^= binaryModulus << (degree - field.Degree);
            }
            int length = RealArithmetic.BitLength(square);
            var coefficients = new BigInteger[length];
            for (int i = 0; i < length; i++)
            {
                token.ThrowIfCancellationRequested();
                coefficients[i] = square.IsEven ? 0 : 1;
                square >>= 1;
            }
            token.ThrowIfCancellationRequested();
            return new FiniteFieldElement(field, coefficients);
        }

        private void PrepareTonelliShanks(CancellationToken token)
        {
            if (nonSquarePower.Field != null) return;
            var nextOddPart = field.Order - 1;
            int nextTwoAdicOrder = 0;
            while (nextOddPart.IsEven) { token.ThrowIfCancellationRequested(); nextOddPart >>= 1; nextTwoAdicOrder++; }
            // In even degree every nonzero prime-subfield element is a square.
            // Start at t (encoding p), without constructing or testing p constants.
            var firstCode = field.Degree % 2 == 0 ? field.Characteristic : BigInteger.One;
            for (var code = firstCode; code < field.Order; code++)
            {
                token.ThrowIfCancellationRequested();
                var coefficients = new List<BigInteger>();
                for (var remaining = code; remaining > 0; remaining /= field.Characteristic)
                {
                    token.ThrowIfCancellationRequested();
                    coefficients.Add(remaining % field.Characteristic);
                }
                var candidate = new FiniteFieldElement(field, coefficients.ToArray());
                if (!field.Pow(candidate, (field.Order - 1) / 2, token).IsOne)
                {
                    var nextNonSquarePower = field.Pow(candidate, nextOddPart, token);
                    oddPart = nextOddPart;
                    twoAdicOrder = nextTwoAdicOrder;
                    nonSquarePower = nextNonSquarePower;
                    return;
                }
            }
            throw new ArithmeticException("No finite-field quadratic nonresidue was found.");
        }

        private int Encode(FiniteFieldElement value)
        {
            int code = 0, prime = (int)field.Characteristic;
            for (int i = value.Coordinates.Length - 1; i >= 0; i--) code = code * prime + (int)value.Coordinates[i];
            return code;
        }

        private FiniteFieldElement Decode(int code)
        {
            if (code == 0) return field.Zero;
            int length = 0, remaining = code, prime = (int)field.Characteristic;
            while (remaining > 0) { length++; remaining /= prime; }
            var coefficients = new BigInteger[length];
            for (int i = 0; i < length; i++) { coefficients[i] = code % prime; code /= prime; }
            return new FiniteFieldElement(field, coefficients);
        }

        private sealed class BinaryLinearMap
        {
            private readonly FiniteField field;
            private readonly BigInteger[] images;
            private readonly BigInteger[] preimages;

            internal BinaryLinearMap(FiniteField field, bool addIdentity, CancellationToken token)
            {
                this.field = field;
                images = new BigInteger[field.Degree];
                preimages = new BigInteger[field.Degree];
                for (int i = 0; i < field.Degree; i++)
                {
                    token.ThrowIfCancellationRequested();
                    var input = Decode(BigInteger.One << i);
                    var square = field.Multiply(input, input, token);
                    var image = Encode(addIdentity ? square + input : square);
                    var preimage = BigInteger.One << i;
                    for (int bit = field.Degree - 1; bit >= 0; bit--)
                    {
                        if ((image & (BigInteger.One << bit)).IsZero) continue;
                        if (images[bit].IsZero) { images[bit] = image; preimages[bit] = preimage; break; }
                        image ^= images[bit]; preimage ^= preimages[bit];
                    }
                }
            }

            internal bool TrySolve(FiniteFieldElement value, CancellationToken token, out FiniteFieldElement root)
            {
                var image = Encode(value); BigInteger preimage = 0;
                for (int bit = images.Length - 1; bit >= 0; bit--)
                {
                    token.ThrowIfCancellationRequested();
                    if ((image & (BigInteger.One << bit)).IsZero) continue;
                    if (images[bit].IsZero) { root = field.Zero; return false; }
                    image ^= images[bit]; preimage ^= preimages[bit];
                }
                root = Decode(preimage);
                token.ThrowIfCancellationRequested();
                return true;
            }

            internal bool HasImage(FiniteFieldElement value, CancellationToken token)
            {
                var image = Encode(value);
                for (int bit = images.Length - 1; bit >= 0; bit--)
                {
                    token.ThrowIfCancellationRequested();
                    if ((image & (BigInteger.One << bit)).IsZero) continue;
                    if (images[bit].IsZero) return false;
                    image ^= images[bit];
                }
                return true;
            }

            private static BigInteger Encode(FiniteFieldElement value)
            {
                BigInteger bits = 0;
                for (int i = value.Coordinates.Length - 1; i >= 0; i--) bits = (bits << 1) + value.Coordinates[i];
                return bits;
            }

            private FiniteFieldElement Decode(BigInteger bits)
            {
                if (bits.IsZero) return field.Zero;
                int length = 0; var remaining = bits;
                while (remaining > 0) { length++; remaining >>= 1; }
                var coefficients = new BigInteger[length];
                for (int i = 0; i < length; i++) { coefficients[i] = bits.IsEven ? 0 : 1; bits >>= 1; }
                return new FiniteFieldElement(field, coefficients);
            }
        }
    }
}

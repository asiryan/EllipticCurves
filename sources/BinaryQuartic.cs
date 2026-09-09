using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    // Binary form a*u^4+b*u^3*v+c*u^2*v^2+d*u*v^3+e*v^4.
    // The formulas are classical invariants/covariants (Cremona III.3.6).
    internal sealed class BinaryQuartic
    {
        internal readonly BigInteger A, B, C, D, E;
        internal BigInteger I => 12 * A * E - 3 * B * D + C * C;
        internal BigInteger J => 72 * A * C * E + 9 * B * C * D - 27 * A * D * D - 27 * B * B * E - 2 * C * C * C;
        internal BigInteger P => 3 * B * B - 8 * A * C;
        internal BigInteger R => B * B * B + 8 * A * A * D - 4 * A * B * C;
        internal BigRational[] Polynomial => new BigRational[] { new BigRational(E), new BigRational(D), new BigRational(C), new BigRational(B), new BigRational(A) };
        internal BinaryQuartic(BigInteger a, BigInteger b, BigInteger c, BigInteger d, BigInteger e)
        { A = a; B = b; C = c; D = d; E = e; }
        internal BigInteger Evaluate(BigInteger u, BigInteger v)
            => (((A * u + B * v) * u + C * v * v) * u + D * v * v * v) * u + E * BigInteger.Pow(v, 4);
        internal BigInteger Evaluate(BigInteger x) => (((A * x + B) * x + C) * x + D) * x + E;
        internal BinaryQuartic Reverse() => new BinaryQuartic(E, D, C, B, A);
        internal BinaryQuartic Scale(BigInteger n) => new BinaryQuartic(n * A, n * B, n * C, n * D, n * E);
        internal bool HasRealPoint(DescentBudget budget)
            => A >= 0 || E >= 0 || DescentPolynomial.HasRealRoot(Polynomial, budget);
        internal bool HasRationalRoot(DescentBudget budget)
            => A.IsZero || DescentPolynomial.HasRationalRoot(Polynomial, budget);

        internal bool Equivalent(BinaryQuartic other, DescentBudget budget)
        {
            budget.Step();
            if (I != other.I || J != other.J) throw new ArgumentException("Quartic invariants must agree before testing equivalence.");
            // Cremona-Fisher (2009), Lemma 10 and the specialization following
            // Theorem 12. The older product-seminvariant criterion is incorrect
            // for a reducible resolvent, even if its auxiliary quartic has a root.
            var a = A; var p = P;
            if (R.IsZero)
            {
                var h = Hessian;
                // A nonzero sextic cannot vanish at all seven affine integers.
                int x = 0;
                while (Sextic(x, 1).IsZero) { budget.Step(); x++; }
                a = Evaluate(x); p = h.Evaluate(x);
            }
            var hh = other.Hessian;
            var test = new BinaryQuartic(a * hh.A - p * other.A, a * hh.B - p * other.B,
                a * hh.C - p * other.C, a * hh.D - p * other.D, a * hh.E - p * other.E);
            return test.HasRationalRoot(budget);
        }

        internal bool TryPoint(DescentBudget budget, out BigInteger u, out BigInteger v, out BigInteger y)
        {
            u = v = y = 0;
            int bound = budget.Options.SearchBound;
            if (bound == 0) return false;
            bool Check(BigInteger x, BigInteger z, out BigInteger squareRoot)
            {
                squareRoot = 0;
                var value = Evaluate(x, z);
                if (value < 0) return false;
                squareRoot = InternalMath.IntegerSqrt(value);
                return squareRoot * squareRoot == value;
            }
            if (!budget.PointStep()) return false;
            if (Check(1, 0, out y)) { u = 1; return true; }
            for (int denominator = 1; denominator <= bound; denominator++)
            for (long numerator = -(long)bound; numerator <= bound; numerator++)
            {
                if (!budget.PointStep()) return false;
                if (BigInteger.GreatestCommonDivisor(numerator, denominator) != 1) continue;
                if (Check(numerator, denominator, out y)) { u = numerator; v = denominator; return true; }
            }
            return false;
        }

        private BinaryQuartic Hessian => new BinaryQuartic(P, 4 * (B * C - 6 * A * D),
            2 * (2 * C * C - 24 * A * E - 3 * B * D), 4 * (C * D - 6 * B * E), 3 * D * D - 8 * C * E);

        private BigInteger Sextic(BigInteger u, BigInteger v)
            => R * BigInteger.Pow(u, 6)
                + 2 * (16 * A * A * E + 2 * A * B * D - 4 * A * C * C + B * B * C) * BigInteger.Pow(u, 5) * v
                + 5 * (8 * A * B * E + B * B * D - 4 * A * C * D) * BigInteger.Pow(u, 4) * v * v
                + 20 * (B * B * E - A * D * D) * BigInteger.Pow(u, 3) * BigInteger.Pow(v, 3)
                - 5 * (8 * A * D * E + B * D * D - 4 * B * C * E) * u * u * BigInteger.Pow(v, 4)
                - 2 * (16 * A * E * E + 2 * B * D * E - 4 * C * C * E + C * D * D) * u * BigInteger.Pow(v, 5)
                - (D * D * D + 8 * B * E * E - 4 * C * D * E) * BigInteger.Pow(v, 6);

        internal (BigRational x, BigRational y) MapPoint(BigInteger u, BigInteger v, BigInteger y)
        {
            if (y.IsZero) throw new ArgumentException("A branch point maps to infinity.");
            var g4 = Hessian.Evaluate(u, v);
            var g6 = Sextic(u, v);
            return (new BigRational(3 * g4, 4 * y * y), new BigRational(27 * g6, 8 * y * y * y));
        }
    }
}

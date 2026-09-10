using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace EllipticCurves
{
    // Division polynomials on y²=x³+Ax+B, represented after removing 2y
    // from every even-indexed psi. Recurrences are identities over Z[A,B,x].
    internal sealed class DivisionPolynomials
    {
        private readonly DescentBudget budget;
        private readonly BigRational[] cubic;
        private readonly Dictionary<int, BigRational[]> psi = new Dictionary<int, BigRational[]>();
        internal DivisionPolynomials(BigRational a, BigRational b, DescentBudget budget)
        {
            this.budget = budget; cubic = new BigRational[] { b, a, 0, 1 };
            psi[0] = new BigRational[] { 0 }; psi[1] = new BigRational[] { 1 }; psi[2] = new BigRational[] { 1 };
            psi[3] = new BigRational[] { -a * a, 12 * b, 6 * a, 0, 3 };
            psi[4] = Scale(new BigRational[] { -a * a * a - 8 * b * b, -4 * a * b, -5 * a * a, 20 * b, 5 * a, 0, 1 }, 2);
        }
        private static BigRational[] Trim(BigRational[] a) { int n = a.Length; while (n > 1 && a[n - 1].IsZero) n--; return a.Take(n).ToArray(); }
        internal BigRational[] Scale(BigRational[] a, BigRational s) => a.Select(x => x * s).ToArray();
        internal BigRational[] Sub(BigRational[] a, BigRational[] b)
        { var c = new BigRational[Math.Max(a.Length, b.Length)]; for (int i = 0; i < c.Length; i++) c[i] = (i < a.Length ? a[i] : 0) - (i < b.Length ? b[i] : 0); return Trim(c); }
        private BigRational[] Mul(BigRational[] a, BigRational[] b)
        {
            var c = new BigRational[a.Length + b.Length - 1];
            for (int i = 0; i < a.Length; i++) for (int j = 0; j < b.Length; j++) { budget.Step(); c[i + j] += a[i] * b[j]; }
            return Trim(c);
        }
        private BigRational[] Square(BigRational[] a) => Mul(a, a);
        private BigRational[] Cube(BigRational[] a) => Mul(Square(a), a);
        private BigRational[] Psi(int n)
        {
            if (psi.TryGetValue(n, out var result)) return result;
            budget.Step(); int m = n / 2;
            if (n % 2 == 0) result = Mul(Psi(m), Sub(Mul(Psi(m + 2), Square(Psi(m - 1))), Mul(Psi(m - 2), Square(Psi(m + 1)))));
            else
            {
                var left = Mul(Psi(m + 2), Cube(Psi(m))); var right = Mul(Psi(m - 1), Cube(Psi(m + 1)));
                var factor = Scale(Square(cubic), 16);
                if (m % 2 == 0) left = Mul(left, factor); else right = Mul(right, factor);
                result = Sub(left, right);
            }
            psi[n] = result; return result;
        }
        internal (BigRational[] numerator, BigRational[] denominator) Multiplication(int n)
        {
            var denominator = Square(Psi(n)); if (n % 2 == 0) denominator = Mul(denominator, Scale(cubic, 4));
            var cross = Mul(Psi(n - 1), Psi(n + 1)); if (n % 2 != 0) cross = Mul(cross, Scale(cubic, 4));
            return (Sub(Mul(new BigRational[] { 0, 1 }, denominator), cross), denominator);
        }
        internal IEnumerable<BigRational> RationalRoots(BigRational[] f)
        {
            f = Trim(f); if (f.Length <= 1) yield break;
            BigInteger denominator = 1; foreach (var a in f) denominator = InternalMath.Lcm(denominator, a.Den);
            var coefficients = f.Select(a => a.Num * (denominator / a.Den)).ToArray(); BigInteger gcd = 0;
            foreach (var a in coefficients) gcd = BigInteger.GreatestCommonDivisor(gcd, a);
            var leading = BigInteger.Abs(coefficients[coefficients.Length - 1] / gcd);
            foreach (var interval in DescentPolynomial.RealRoots(f, new BigRational(1, 4 * leading), budget))
                for (var n = DescentPolynomial.Ceiling(interval.Lower * new BigRational(leading)); n <= DescentPolynomial.Floor(interval.Upper * new BigRational(leading)); n++)
                { budget.Step(); var x = new BigRational(n, leading); if (DescentPolynomial.Evaluate(f, x).IsZero) yield return x; }
        }
    }
}

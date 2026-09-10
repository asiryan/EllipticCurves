using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Small-degree polynomial arithmetic over Q. All root decisions use Sturm sequences,
    // including repeated roots of the quartic equivalence polynomial.
    internal static class DescentPolynomial
    {
        internal static BigInteger Floor(BigRational x)
        {
            var q = BigInteger.DivRem(x.Num, x.Den, out var r);
            return r.Sign < 0 ? q - 1 : q;
        }
        internal static BigInteger Ceiling(BigRational x) => -Floor(-x);
        internal static BigRational Evaluate(BigRational[] f, BigRational x)
        {
            var y = BigRational.Zero;
            for (int i = f.Length - 1; i >= 0; i--) y = y * x + f[i];
            return y;
        }
        private static BigRational[] Trim(BigRational[] f)
        {
            int n = f.Length;
            while (n > 0 && f[n - 1].IsZero) n--;
            return f.Take(n).ToArray();
        }
        private static BigRational[] Derivative(BigRational[] f)
            => f.Skip(1).Select((a, i) => a * (i + 1)).ToArray();
        private static BigRational[] Quotient(BigRational[] f, BigRational[] g, CancellationToken token)
        {
            var r = (BigRational[])f.Clone();
            var result = new BigRational[f.Length - g.Length + 1];
            for (int i = r.Length - 1; i >= g.Length - 1; i--)
            {
                token.ThrowIfCancellationRequested();
                var q = result[i - g.Length + 1] = r[i] / g[g.Length - 1];
                for (int j = 0; j < g.Length; j++) r[i - g.Length + 1 + j] -= q * g[j];
            }
            if (Trim(r).Length != 0) throw new InvalidOperationException("Inexact polynomial division.");
            return Trim(result);
        }
        private static BigRational[] Remainder(BigRational[] f, BigRational[] g, CancellationToken token)
        {
            var r = (BigRational[])f.Clone();
            for (int i = r.Length - 1; i >= g.Length - 1; i--)
            {
                token.ThrowIfCancellationRequested();
                var q = r[i] / g[g.Length - 1];
                for (int j = 0; j < g.Length; j++) r[i - g.Length + 1 + j] -= q * g[j];
            }
            return Trim(r);
        }
        private static List<BigRational[]> Sturm(BigRational[] f, CancellationToken token)
        {
            var sequence = new List<BigRational[]> { Trim(f) };
            if (sequence[0].Length <= 1) return sequence;
            sequence.Add(Trim(Derivative(sequence[0])));
            while (sequence[sequence.Count - 1].Length > 1)
            {
                var r = Remainder(sequence[sequence.Count - 2], sequence[sequence.Count - 1], token);
                if (r.Length == 0) return Sturm(Quotient(sequence[0], sequence[sequence.Count - 1], token), token);
                // Positive normalization limits coefficient growth and preserves variations.
                var scale = r[r.Length - 1].Abs();
                sequence.Add(r.Select(x => -x / scale).ToArray());
            }
            return sequence;
        }
        private static int Variations(List<BigRational[]> sequence, BigRational x)
        {
            int previous = 0, count = 0;
            foreach (var f in sequence)
            {
                int sign = Evaluate(f, x).Sign;
                if (sign == 0) continue;
                if (previous != 0 && previous != sign) count++;
                previous = sign;
            }
            return count;
        }
        internal static List<RationalInterval> RealRoots(BigRational[] f, BigRational width, DescentBudget budget)
        {
            f = Trim(f);
            var roots = new List<RationalInterval>();
            if (f.Length <= 1) return roots;
            budget.Token.ThrowIfCancellationRequested();
            var sequence = Sturm(f, budget.Token);
            BigInteger bound = 1;
            for (int i = 0; i < f.Length - 1; i++)
                bound = BigInteger.Max(bound, 1 + Ceiling((f[i] / f[f.Length - 1]).Abs()));
            var left = new BigRational(-bound); var right = new BigRational(bound);
            var stack = new Stack<(BigRational l, BigRational r, int vl, int vr)>();
            stack.Push((left, right, Variations(sequence, left), Variations(sequence, right)));
            while (stack.Count != 0)
            {
                budget.Step();
                var s = stack.Pop(); int count = s.vl - s.vr;
                if (count == 0) continue;
                if (count == 1 && Evaluate(f, s.r).IsZero)
                { roots.Add(RationalInterval.Exact(s.r)); continue; }
                if (count == 1 && s.r - s.l <= width)
                { roots.Add(new RationalInterval(s.l, s.r)); continue; }
                var mid = (s.l + s.r) / 2;
                int vm = Variations(sequence, mid);
                stack.Push((mid, s.r, vm, s.vr));
                stack.Push((s.l, mid, s.vl, vm));
            }
            return roots;
        }
        internal static bool HasRealRoot(BigRational[] f, DescentBudget budget)
            => RealRoots(f, new BigRational(BigInteger.One << 16), budget).Count != 0;

        internal static bool HasRationalRoot(BigRational[] f, DescentBudget budget)
        {
            f = Trim(f);
            if (f.Length <= 1) return false;
            if (f[0].IsZero) return true;
            BigInteger denominator = 1;
            foreach (var a in f) denominator = InternalMath.Lcm(denominator, a.Den);
            var coefficients = f.Select(a => a.Num * (denominator / a.Den)).ToArray();
            BigInteger gcd = 0;
            foreach (var a in coefficients) gcd = BigInteger.GreatestCommonDivisor(gcd, a);
            var leading = BigInteger.Abs(coefficients[coefficients.Length - 1] / gcd);
            // Rational root theorem: leading*x is integral. Root isolation avoids
            // factoring the often very large constant of an equivalence polynomial.
            foreach (var root in RealRoots(f, new BigRational(1, 4 * leading), budget))
                for (var n = Ceiling(root.Lower * new BigRational(leading)); n <= Floor(root.Upper * new BigRational(leading)); n++)
                {
                    budget.Step();
                    if (Evaluate(f, new BigRational(n, leading)).IsZero) return true;
                }
            return false;
        }
    }
}

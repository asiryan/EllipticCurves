using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    internal static class QuarticLocalSolubility
    {
        internal static bool Everywhere(BinaryQuartic q, IEnumerable<BigInteger> badPrimes, DescentBudget budget)
        {
            if (!q.HasRealPoint(budget)) return false;
            foreach (var p in badPrimes)
                if (!AtPrime(q, p, budget)) return false;
            return true;
        }

        // Callers supply proved primes and nonsingular binary forms. Odd primes of
        // good reduction need no test: the smooth genus-one curve over F_p has a point.
        internal static bool AtPrime(BinaryQuartic q, BigInteger p, DescentBudget budget)
        {
            if (p < 2) throw new ArgumentOutOfRangeException(nameof(p));
            return InChart(q, p, false, budget) || InChart(q.Reverse(), p, true, budget);
        }

        private static int Valuation(BigInteger n, BigInteger p)
        {
            if (n.IsZero) return int.MaxValue;
            int v = 0;
            while (n % p == 0) { n /= p; v++; }
            return v;
        }
        private static bool InChart(BinaryQuartic q, BigInteger p, bool atInfinity, DescentBudget budget)
        {
            // DFS over p-adic balls, with one pending sibling frame per level.
            // The second chart is u=1, v in p Z_p; it includes the point at infinity.
            var stack = new Stack<(BigInteger x, BigInteger modulus, int depth, BigInteger digit)>();
            stack.Push((0, atInfinity ? p : BigInteger.One, atInfinity ? 1 : 0, -1));
            while (stack.Count != 0)
            {
                budget.Step();
                var ball = stack.Pop();
                if (ball.digit < 0)
                {
                    int decision = Decide(q, ball.x, p, ball.depth);
                    if (decision > 0) return true;
                    if (decision < 0) continue;
                    ball.digit = 0;
                }
                if (ball.digit == p) continue;
                stack.Push((ball.x, ball.modulus, ball.depth, ball.digit + 1));
                stack.Push((ball.x + ball.digit * ball.modulus, ball.modulus * p, checked(ball.depth + 1), -1));
            }
            return false;
        }

        private static int Decide(BinaryQuartic q, BigInteger x, BigInteger p, int depth)
        {
            var value = q.Evaluate(x);
            if (value.IsZero) return 1;
            int valuation = Valuation(value, p);
            var unit = value / BigInteger.Pow(p, valuation);
            bool even = valuation % 2 == 0;
            if (even && (p == 2 ? NativeNumberTheory.Mod(unit, 8) == 1
                : BigInteger.ModPow(NativeNumberTheory.Mod(unit, p), (p - 1) / 2, p) == 1)) return 1;

            var derivative = ((4 * q.A * x + 3 * q.B) * x + 2 * q.C) * x + q.D;
            int dv = Valuation(derivative, p);
            // Hensel produces an exact zero, hence a point with y=0. The zero is
            // congruent to x mod p, so also belongs to the infinity chart if used.
            if (valuation > 2L * dv) return 1;

            // Taylor expansion on x+p^depth Z_p bounds the valuation of every
            // change in f. Once valuation and square-class unit are fixed, the
            // entire ball is nonsquare. Otherwise subdivide; no finite sieve is
            // promoted to a positive local-solubility result.
            long variation = Math.Min((long)dv + depth,
                Math.Min((long)Valuation(6 * q.A * x * x + 3 * q.B * x + q.C, p) + 2L * depth,
                Math.Min((long)Valuation(4 * q.A * x + q.B, p) + 3L * depth,
                    (long)Valuation(q.A, p) + 4L * depth)));
            if (variation <= valuation) return 0;
            if (!even || p != 2) return -1;
            if (variation >= valuation + 3L) return -1;
            if (variation >= valuation + 2L && NativeNumberTheory.Mod(unit, 4) == 3) return -1;
            return 0;
        }
    }
}

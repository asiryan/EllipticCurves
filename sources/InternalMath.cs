using System;
using System.Collections.Generic;
using System.Numerics;

namespace EllipticCurves
{
    /// <summary>
    /// Internal numeric helpers used across the library.
    /// 
    /// Scope:
    /// • Invariant computations for integral Weierstrass models (c4, c6, Δ).
    /// • Q–isomorphism checks via scaling of invariants (u^4, u^6, u^12).
    /// • Exact k-th roots over ℚ and ℤ (Newton / integer root tests).
    /// • Fast torsion tests (order divisibility and small Mazur fallback).
    /// • Small finite-field helpers: Legendre symbol, modular exponentiation,
    ///   and counting points on short Weierstrass curves over 𝔽_p.
    /// • Basic integer factorization (Pollard–Rho + Miller–Rabin) to build
    ///   square divisors of |Δ| and enumerations of divisors.
    /// 
    /// Notes:
    /// • All methods are deterministic and allocation-light.
    /// • Big-integer operations can be expensive for huge inputs—these are
    ///   intended for typical arithmetic of elliptic curves over ℚ.
    /// • If you target an older C# language version, replace collection
    ///   expressions like `int[] small = [2,3,...]` with classic initializers
    ///   `new int[] { 2, 3, ... }`.
    /// </summary>
    internal static partial class InternalMath
    {
        /// <summary>
        /// Compute integral invariants (c4, c6, Δ) from integral a-invariants [a1,a2,a3,a4,a6].
        /// Assumes all inputs are integers (integral model). No normalization/scaling is applied here.
        /// </summary>
        public static (BigInteger c4, BigInteger c6, BigInteger Delta) InvariantsIntFromAinvs(BigInteger a1, 
            BigInteger a2, BigInteger a3, BigInteger a4, BigInteger a6)
        {
            BigInteger b2 = a1 * a1 + 4 * a2;
            BigInteger b4 = 2 * a4 + a1 * a3;
            BigInteger b6 = a3 * a3 + 4 * a6;
            BigInteger b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
            BigInteger Delta = -(b2 * b2) * b8 - 8 * b4 * b4 * b4 - 27 * b6 * b6 + 9 * b2 * b4 * b6;
            BigInteger c4 = b2 * b2 - 24 * b4;
            BigInteger c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6;
            return (c4, c6, Delta);
        }

        /// <summary>
        /// Build an integral Weierstrass model of curve C by the scaling (X,Y) = (L^2 x, L^3 y),
        /// where L = lcm(denominators of a_i). Returns the integer invariants (c4,c6,Δ).
        /// </summary>
        public static (BigInteger c4, BigInteger c6, BigInteger Delta) IntegralInvariants(EllipticCurveQ C)
        {
            // L = lcm of denominators of a1..a6
            BigInteger L = InternalMath.Lcm(C.A1.Den,
                           InternalMath.Lcm(C.A2.Den,
                           InternalMath.Lcm(C.A3.Den,
                           InternalMath.Lcm(C.A4.Den, C.A6.Den))));

            // a_i' = L^i * a_i  (become integers)
            BigInteger L1 = L;
            BigInteger L2 = BigInteger.Pow(L, 2);
            BigInteger L3 = BigInteger.Pow(L, 3);
            BigInteger L4 = BigInteger.Pow(L, 4);
            BigInteger L6 = BigInteger.Pow(L, 6);

            BigInteger al1 = C.A1.Num * L1 / C.A1.Den;
            BigInteger al2 = C.A2.Num * L2 / C.A2.Den;
            BigInteger al3 = C.A3.Num * L3 / C.A3.Den;
            BigInteger al4 = C.A4.Num * L4 / C.A4.Den;
            BigInteger al6 = C.A6.Num * L6 / C.A6.Den;

            return InternalMath.InvariantsIntFromAinvs(al1, al2, al3, al4, al6);
        }

        /// <summary>
        /// Check whether two elliptic curves are Q–isomorphic by testing the scaling
        /// relations on invariants: c4_E = u^4 c4_C, c6_E = u^6 c6_C, Δ_E = u^12 Δ_C for some u ∈ ℚ.
        /// Outputs the scaling factor u if successful.
        /// </summary>
        public static bool IsQIsomorphic(BigRational c4E, BigRational c6E, BigRational dE, 
            BigRational rc4C, BigRational rc6C, BigRational rdC, out BigRational u)
        {
            u = default;

            bool c4zero = c4E.IsZero || rc4C.IsZero;
            bool c6zero = c6E.IsZero || rc6C.IsZero;

            // Prefer using c4 when available; otherwise fall back to c6 test.
            if (!c4zero)
            {
                var ratio4 = c4E / rc4C;                  // should equal u^4
                if (!TryRationalKthRoot(ratio4, 4, out var u4)) return false;
                u = u4;

                if (!c6zero && BigRational.Pow(u, 6) != (c6E / rc6C)) return false;
                if (BigRational.Pow(u, 12) != (dE / rdC)) return false;
                return true;
            }
            else if (!c6zero)
            {
                var ratio6 = c6E / rc6C;                  // should equal u^6
                if (!TryRationalKthRoot(ratio6, 6, out var u6)) return false;
                u = u6;

                if (BigRational.Pow(u, 12) != (dE / rdC)) return false;
                return true;
            }
            else
            {
                // For nonsingular curves, c4 and c6 cannot both be zero.
                return false;
            }
        }

        /// <summary>
        /// Try to compute a rational k-th root of r ∈ ℚ (k ≥ 2). Returns true iff
        /// r = root^k for some rational root in canonical form.
        /// Implementation: require both |Num| and Den to be perfect k-th powers in ℤ.
        /// </summary>
        public static bool TryRationalKthRoot(BigRational r, int k, out BigRational root)
        {
            root = default;
            if (r.Sign < 0 && (k % 2 == 0)) return false; // even root of a negative number is not rational real

            var an = BigInteger.Abs(r.Num);
            var ad = r.Den;

            if (!TryIntegerKthRoot(an, k, out var rn)) return false;
            if (!TryIntegerKthRoot(ad, k, out var rd)) return false;

            if (r.Sign < 0) rn = BigInteger.Negate(rn);
            root = new BigRational(rn, rd);
            return BigRational.Pow(root, k) == r; // exact verification
        }

        /// <summary>
        /// Try to compute exact integer k-th root: return true iff n = rt^k for some integer rt (rt ≥ 0).
        /// Otherwise, returns false and rt = floor(n^(1/k)).
        /// </summary>
        public static bool TryIntegerKthRoot(BigInteger n, int k, out BigInteger rt)
        {
            rt = IntegerKthRoot(n, k);
            return BigInteger.Pow(rt, k) == n;
        }

        /// <summary>
        /// floor(n^(1/k)) for n ≥ 0 using Newton's method generalized to k-th roots.
        /// For small n it exits fast; for large n it converges monotonically.
        /// </summary>
        public static BigInteger IntegerKthRoot(BigInteger n, int k)
        {
            if (n.IsZero) return BigInteger.Zero;
            if (n.IsOne) return BigInteger.One;

            // Newton's decreasing iteration must start above the root, including when k does not divide 8.
            int bits = n.ToByteArray().Length * 8;
            BigInteger x = BigInteger.One << Math.Max(1, (bits + k - 1) / k);

            while (true)
            {
                // x_{t+1} = ((k-1)*x + n/x^{k-1}) / k
                var xk_1 = BigInteger.Pow(x, k - 1);
                if (xk_1.IsZero) break;
                var next = ((k - 1) * x + n / xk_1) / k;
                if (next >= x) return x; // converged
                x = next;
            }
            return x;
        }

        /// <summary>
        /// Check whether a point is torsion by first testing divisibility by candidate orders
        /// (coming from gcd of #E(𝔽_p) over a few good primes), and if none matches,
        /// falling back to a bounded exact check up to 12 (Mazur’s bound).
        /// </summary>
        public static bool IsTorsionWithCandidates(EllipticCurveQ E, EllipticCurvePoint P, List<int> orders)
        {
            // Fast path: only test divisors of the gcd of #E(F_p) for selected primes.
            for (int i = 0; i < orders.Count; i++)
            {
                int n = orders[i];
                if (n <= 1) continue;
                if (OrderDivides(E, P, n)) return true;
            }
            // Safe fallback: exact check up to Mazur bound (cheap for the few points we test).
            var Q = P;
            for (int k = 1; k <= 12; k++)
            {
                if (Q.IsInfinity) return true; // order(P) | k
                Q = E.Add(Q, P);
            }
            return false;
        }

        /// <summary>
        /// Return true iff n*P = O (the point at infinity), i.e. the order of P divides n.
        /// Computed with double-and-add (no precomputation).
        /// </summary>
        public static bool OrderDivides(EllipticCurveQ E, EllipticCurvePoint P, int n)
        {
            var Q = E.Multiply(P, new BigInteger(n));
            return Q.IsInfinity;
        }

        /// <summary>
        /// Least common multiple for nonzero BigIntegers. If either argument is 0,
        /// this returns 0 by the standard convention.
        /// </summary>
        public static BigInteger Lcm(BigInteger a, BigInteger b)
        {
            if (a.IsZero || b.IsZero) return BigInteger.Zero;
            var gcd = BigInteger.GreatestCommonDivisor(a, b);
            var lcm = a / gcd * b;
            return BigInteger.Abs(lcm);
        }

        /// <summary>
        /// Evaluate the cubic X^3 + A X + B at integer X.
        /// Used for 2-torsion and Lutz–Nagell searches on short integral models.
        /// </summary>
        public static BigInteger EvalCubic(BigInteger A, BigInteger B, BigInteger x)
        {
            return x * x * x + A * x + B;
        }

        /// <summary>
        /// floor(sqrt(n)) for n ≥ 0 via a monotone Newton iteration specialized to k=2.
        /// </summary>
        public static BigInteger IntegerSqrt(BigInteger n)
        {
            if (n <= 1) return n;
            BigInteger x0 = n, x1 = (n >> 1) + 1;
            while (x1 < x0) { x0 = x1; x1 = (x1 + n / x1) >> 1; }
            return x0;
        }

        /// <summary>
        /// Count points on a short Weierstrass curve Y^2 = X^3 + A X + B over 𝔽_p (p odd and small).
        /// Uses Legendre symbol to count square/non-square values of RHS for each X.
        /// </summary>
        public static int CountPointsFpShort(BigInteger A, BigInteger B, int p)
        {
            int cnt = 1; // point at infinity
            BigInteger Amod = A % p;
            if (Amod < 0) Amod += p;
            BigInteger Bmod = B % p;
            if (Bmod < 0) Bmod += p;
            for (int x = 0; x < p; x++)
            {
                BigInteger xVal = x;
                BigInteger rhs = xVal * xVal % p;                  // x^2
                rhs = (rhs * xVal + Amod * xVal) % p;              // x^3 + A x
                rhs = (rhs + Bmod) % p;                            // + B
                if (rhs < 0) rhs += p;

                int chi = Legendre((int)rhs, p); // 0, 1, -1
                if (chi == 0) cnt += 1;          // one solution y=0
                else if (chi == 1) cnt += 2;     // two square roots mod p
            }
            return cnt;
        }

        /// <summary>
        /// Legendre symbol (a|p) for odd prime p: returns 0 if a ≡ 0 (mod p), 1 if a is a QR, −1 otherwise.
        /// Euler’s criterion: a^((p−1)/2) ≡ (a|p) (mod p).
        /// </summary>
        public static int Legendre(int a, int p)
        {
            if (a % p == 0) return 0;
            int e = (p - 1) / 2;
            int r = ModPow(a, e, p);
            if (r == 1) return 1;
            if (r == p - 1) return -1;
            return 0;
        }

        /// <summary>
        /// Fast modular exponentiation (a^e mod m) with 32-bit ints.
        /// </summary>
        public static int ModPow(int a, int e, int m)
        {
            long res = 1, b = ((a % m) + m) % m;
            while (e > 0)
            {
                if ((e & 1) == 1) res = res * b % m;
                b = b * b % m;
                e >>= 1;
            }
            return (int)res;
        }

        /// <summary>
        /// Return certified prime factors with exponents for |n|, using the shared native factorizer.
        /// Zero retains the legacy empty-result convention; divisor enumeration handles it separately.
        /// </summary>
        public static Dictionary<BigInteger, int> FactorAbs(BigInteger n)
        {
            return n.IsZero ? new Dictionary<BigInteger, int>() : NativeNumberTheory.Factor(n, default);
        }

        /// <summary>
        /// Legacy helper, not used by FactorAbs: splits composite n into probable prime factors using
        /// Miller–Rabin (probable prime) and Pollard–Rho to find nontrivial divisors.
        /// </summary>
        public static void FactorRec(BigInteger n, Dictionary<BigInteger, int> res)
        {
            if (n == 1) return;
            if (IsProbablePrime(n))
            {
                int e;
                if (res.TryGetValue(n, out e)) res[n] = e + 1; else res[n] = 1;
                return;
            }
            var d = PollardRho(n);
            FactorRec(d, res);
            FactorRec(n / d, res);
        }

        /// <summary>
        /// Probable-prime test: quick small trial division, then Miller–Rabin with a fixed base set.
        /// Passing this fixed base set does not prove primality, even across the full 64-bit range.
        /// </summary>
        public static bool IsProbablePrime(BigInteger n)
        {
            if (n < 2) return false;

            // NOTE: if your target C# version doesn’t support collection expressions,
            // replace with: new int[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 }
            int[] small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
            for (int i = 0; i < small.Length; i++)
            {
                int p = small[i];
                if (n == p) return true;
                if (n % p == 0) return n == p;
            }

            // Legacy probable-prime check only. Certified factorization uses NativeNumberTheory.
            int[] bases = [2, 3, 5, 7, 11, 13, 17];
            BigInteger d = n - 1; int s = 0;
            while ((d & 1) == 0) { d >>= 1; s++; }
            for (int i = 0; i < bases.Length; i++)
            {
                int a = bases[i];
                if (n == a) return true;
                if (MillerRabinCheck(n, d, s, a) == false) return false;
            }
            return true;
        }

        /// <summary>
        /// One Miller–Rabin round for given odd n, with factorization n−1 = d·2^s and base a.
        /// Returns false if n is definitely composite, true if it passes this round.
        /// </summary>
        public static bool MillerRabinCheck(BigInteger n, BigInteger d, int s, int a)
        {
            BigInteger x = BigInteger.ModPow(a, d, n);
            if (x == 1 || x == n - 1) return true;
            for (int r = 1; r < s; r++)
            {
                x = x * x % n;
                if (x == n - 1) return true;
            }
            return false;
        }

        /// <summary>
        /// Pollard–Rho with a simple f(x)=x^2+c map and Brent-like cycle detection.
        /// Returns a nontrivial divisor of odd composite n (heuristic, but very effective in practice).
        /// </summary>
        public static BigInteger PollardRho(BigInteger n)
        {
            if ((n & 1) == 0) return 2;
            var rnd = new Random(1234567);
            while (true)
            {
                BigInteger c = RandomBelow(n, rnd);
                BigInteger x = RandomBelow(n, rnd);
                BigInteger y = x;
                BigInteger d = 1;
                while (d == 1)
                {
                    x = (x * x + c) % n;
                    y = (y * y + c) % n;
                    y = (y * y + c) % n; // 2 steps
                    d = BigInteger.GreatestCommonDivisor(BigInteger.Abs(x - y), n);
                }
                if (d != n) return d; // nontrivial factor found
            }
        }

        /// <summary>
        /// Return a random integer r with 0 &lt; r &lt; n using the provided PRNG.
        /// The distribution is not cryptographically secure (intended for math utilities).
        /// </summary>
        public static BigInteger RandomBelow(BigInteger n, Random rng)
        {
            var bytes = n.ToByteArray();
            BigInteger r;
            do
            {
                rng.NextBytes(bytes);
                bytes[bytes.Length - 1] &= 0x7F; // keep positive
                r = new BigInteger(bytes);
            } while (r <= 0 || r >= n);
            return r;
        }

        /// <summary>
        /// Enumerate all square divisors y2 of |Δ| using its prime factorization:
        /// if |Δ|=∏ p_i^{e_i}, then y2 ranges over ∏ p_i^{2f_i} with 0≤f_i≤⌊e_i/2⌋.
        /// </summary>
        public static IEnumerable<BigInteger> EnumerateSquareDivisors(Dictionary<BigInteger, int> fact)
        {
            var primes = new List<BigInteger>(fact.Keys);
            var exps = new List<int>(primes.Count);
            for (int i = 0; i < primes.Count; i++) exps.Add(fact[primes[i]] / 2);

            // Build all products Π p_i^{2f_i} iteratively.
            var accs = new List<BigInteger> { BigInteger.One };
            for (int i = 0; i < primes.Count; i++)
            {
                var next = new List<BigInteger>();
                BigInteger p = primes[i];
                int e = exps[i];
                for (int j = 0; j < accs.Count; j++)
                {
                    BigInteger baseVal = accs[j];
                    BigInteger pow = BigInteger.One;
                    for (int t = 0; t <= e; t++)
                    {
                        next.Add(baseVal * pow);
                        pow *= p * p;
                    }
                }
                accs = next;
            }
            for (int i = 0; i < accs.Count; i++) yield return accs[i];
        }

        /// <summary>
        /// Enumerate all positive divisors of |n| from its factorization (includes 1 and |n|).
        /// Returns 0 only if input n was 0 (by convention).
        /// </summary>
        public static IEnumerable<BigInteger> EnumerateDivisorsAbs(BigInteger n)
        {
            n = n >= 0 ? n : -n;
            if (n.IsZero) { yield return BigInteger.Zero; yield break; }

            var fact = FactorAbs(n);
            var primes = new List<BigInteger>(fact.Keys);
            var exps = new List<int>(primes.Count);
            for (int i = 0; i < primes.Count; i++) exps.Add(fact[primes[i]]);

            var accs = new List<BigInteger> { BigInteger.One };
            for (int i = 0; i < primes.Count; i++)
            {
                var next = new List<BigInteger>();
                BigInteger p = primes[i];
                int e = exps[i];
                for (int j = 0; j < accs.Count; j++)
                {
                    BigInteger baseVal = accs[j];
                    BigInteger pow = BigInteger.One;
                    for (int t = 0; t <= e; t++)
                    {
                        next.Add(baseVal * pow);
                        pow *= p;
                    }
                }
                accs = next;
            }
            for (int i = 0; i < accs.Count; i++) yield return accs[i];
        }
    }
}

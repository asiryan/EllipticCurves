using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;

namespace EllipticCurves
{
    // Self-initializing quadratic sieve, with CRT polynomial families, Gray-code
    // root updates and paired large cofactors. See docs/native-arithmetic.md.
    // Floating-point logarithms select candidates only; relations and divisors
    // are computed with exact integers. No primality decision is made here.
    internal sealed class NativeQuadraticSieve
    {
        private readonly BigInteger n, sievedNumber;
        private readonly CancellationToken token;
        private readonly List<BasePrime> factorBase = new List<BasePrime>();
        private readonly List<SievePower> powers = new List<SievePower>();
        private readonly Dictionary<BigInteger, Relation> partials = new Dictionary<BigInteger, Relation>();
        private readonly List<Relation> relations = new List<Relation>();
        private readonly NativeQuadraticSieve relationOwner;
        private readonly object relationLock = new object();
        private readonly Row[] basis;
        private readonly int radius, words;
        private readonly long largeBound;
        private readonly float[] scores;
        // Contributions of the smallest prime powers repeat over this wheel.
        // A CRT sign change rotates the table; block copies replace most sieve writes.
        private const int PresieveWheel = 32 * 9 * 5 * 7;
        private readonly float[] presieve = new float[PresieveWheel];
        private int[] presieveShifts;
        private int presieveOffset;
        private readonly BigInteger targetA;
        private BigInteger a, b;
        private BigInteger[] crtParts;
        private int[] aFactors;
        private int polynomialIndex, familySize;
        private uint randomState = 0x9e3779b9;

        private sealed class BasePrime
        {
            internal int P, SquareRoot, Root1, Root2;
            internal float Log;
        }

        private sealed class SievePower
        {
            internal int Modulus, PrimeIndex, InverseA;
            internal int[] SquareRoots, Roots, Shifts;
            internal float Log;
        }

        private sealed class Relation
        {
            internal BigInteger X, Square;
            // Factor-base indices, repeated with multiplicity; zero means -1.
            internal int[] Factors;
        }

        private sealed class Row
        {
            internal ulong[] Parity, Combination;
        }

        internal static BigInteger FindDivisor(BigInteger n, CancellationToken token, int maxWorkers = 0)
        {
            token.ThrowIfCancellationRequested();
            if (maxWorkers < 0) throw new ArgumentOutOfRangeException(nameof(maxWorkers));
            int digits = (int)(BigInteger.Log10(n) + 1);
            int bound = digits < 30 ? 2000 : digits < 40 ? 8000 : digits < 50 ? 25000
                : digits < 60 ? 60000 : digits < 70 ? 150000 : 300000;
            var primes = NativeNumberTheory.SievePrimes(bound);
            foreach (int p in primes)
            {
                token.ThrowIfCancellationRequested();
                if (n % p == 0) return p;
            }
            int multiplier = ChooseMultiplier(n, primes);
            int workers = WorkerCount(digits, maxWorkers, Environment.ProcessorCount);
            if (workers == 1) return new NativeQuadraticSieve(n, multiplier, primes, token).Run();
            using (var stop = CancellationTokenSource.CreateLinkedTokenSource(token))
            {
                var owner = new NativeQuadraticSieve(n, multiplier, primes, stop.Token);
                BigInteger result = 1;
                try
                {
                    Parallel.For(0, workers, new ParallelOptions
                    { MaxDegreeOfParallelism = workers, CancellationToken = stop.Token }, index =>
                    {
                        try
                        {
                            var worker = index == 0 ? owner : new NativeQuadraticSieve(n, multiplier, primes, stop.Token, owner);
                            worker.randomState ^= unchecked((uint)index * 0x85ebca6bU);
                            var divisor = worker.Run();
                            lock (owner.relationLock) if (result.IsOne) result = divisor;
                            stop.Cancel();
                        }
                        catch
                        {
                            // A failed worker must stop its peers before Parallel.For joins them.
                            stop.Cancel();
                            throw;
                        }
                    });
                }
                catch (OperationCanceledException) when (!token.IsCancellationRequested && result > 1) { }
                token.ThrowIfCancellationRequested();
                return result;
            }
        }

        internal static int WorkerCount(int digits, int requested, int available)
        {
            if (requested < 0) throw new ArgumentOutOfRangeException(nameof(requested));
            if (available < 1) throw new ArgumentOutOfRangeException(nameof(available));
            // Small jobs cannot amortize extra sieve instances. Large residuals
            // can use all CPUs, while explicit limits also reach recursive proofs.
            if (digits < 45) return 1;
            return Math.Min(available, requested == 0 ? (digits < 70 ? 4 : available) : requested);
        }

        private NativeQuadraticSieve(BigInteger n, int multiplier, int[] primes, CancellationToken token,
            NativeQuadraticSieve owner = null)
        {
            this.n = n;
            this.token = token;
            relationOwner = owner ?? this;
            sievedNumber = n * multiplier;
            radius = primes[primes.Length - 1] < 25000 ? 16384 : 65536;
            scores = new float[2 * radius + 1];
            if (owner == null) foreach (int p in primes)
            {
                token.ThrowIfCancellationRequested();
                int residue = (int)(sievedNumber % p);
                if (p > 2 && residue != 0 && InternalMath.Legendre(residue, p) != 1) continue;
                factorBase.Add(new BasePrime { P = p, SquareRoot = SquareRootModPrime(residue, p), Log = (float)Math.Log(p) });
            }
            else foreach (var prime in owner.factorBase)
                factorBase.Add(new BasePrime { P = prime.P, SquareRoot = prime.SquareRoot, Log = prime.Log });
            basis = new Row[factorBase.Count + 1];
            words = (basis.Length + 1 + 63) / 64;
            largeBound = (long)primes[primes.Length - 1] * 64;
            targetA = InternalMath.IntegerSqrt(2 * sievedNumber) / radius;
            if (owner == null) PreparePrimePowers();
            else foreach (var power in owner.powers)
                powers.Add(new SievePower { Modulus = power.Modulus, PrimeIndex = power.PrimeIndex,
                    SquareRoots = power.SquareRoots, Log = power.Log });
        }

        private static int ChooseMultiplier(BigInteger n, int[] primes)
        {
            int best = 1;
            double bestScore = double.NegativeInfinity;
            foreach (int k in new[] { 1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73 })
            {
                var kn = n * k;
                int mod8 = (int)(kn % 8);
                double score = -0.5 * Math.Log(k) + Math.Log(2) * (mod8 == 1 ? 2 : mod8 == 5 ? 1 : 0.5);
                foreach (int p in primes)
                {
                    if (p == 2) continue;
                    if (p > 1000) break;
                    int residue = (int)(kn % p);
                    if (residue == 0) score += Math.Log(p) / p;
                    else if (InternalMath.Legendre(residue, p) == 1) score += 2 * Math.Log(p) / (p - 1);
                }
                if (score > bestScore) { bestScore = score; best = k; }
            }
            return best;
        }

        private BigInteger Run()
        {
            while (true)
            {
                token.ThrowIfCancellationRequested();
                // A is a product of distinct factor-base primes. Changing one
                // CRT sign changes B by +/-2*Bi, so all nonsingular sieve roots
                // move by precomputed offsets instead of being recomputed.
                int changedPart = -1, direction = 0;
                if (polynomialIndex == familySize) StartFamily();
                else
                {
                    int gray = polynomialIndex ^ (polynomialIndex >> 1);
                    int previous = (polynomialIndex - 1) ^ ((polynomialIndex - 1) >> 1);
                    int changed = gray ^ previous;
                    changedPart = 0;
                    while ((changed & (1 << changedPart)) == 0) changedPart++;
                    direction = (gray & changed) != 0 ? 1 : -1;
                    b -= direction * 2 * crtParts[changedPart];
                }
                polynomialIndex++;
                BigInteger divisor;
                var c = (b * b - sievedNumber) / a;
                SievePolynomial(c, changedPart, direction);
                // Extra tolerance covers log rounding and prime powers beyond
                // the sieve interval. Missing a candidate affects speed only.
                double threshold = BigInteger.Log(BigInteger.Abs(c) + a * radius * radius
                    + 2 * BigInteger.Abs(b) * radius) - Math.Log(largeBound) - 2;
                for (int index = 0; index < scores.Length; index++)
                {
                    if ((index & 1023) == 0) token.ThrowIfCancellationRequested();
                    if (scores[index] < threshold) continue;
                    int x = index - radius;
                    var value = (a * x + 2 * b) * x + c;
                    var factors = new List<int>();
                    if (value.Sign < 0) { factors.Add(0); value = -value; }
                    if (value.IsZero)
                    {
                        divisor = BigInteger.GreatestCommonDivisor(a * x + b, n);
                        if (divisor > 1 && divisor < n) return divisor;
                        continue;
                    }
                    for (int i = 0; i < factorBase.Count; i++)
                    {
                        var prime = factorBase[i];
                        int residue = index % prime.P;
                        if (residue != prime.Root1 && residue != prime.Root2) continue;
                        while (value % prime.P == 0) { factors.Add(i + 1); value /= prime.P; }
                    }
                    if (value > largeBound) continue;
                    factors.AddRange(aFactors);
                    var relation = new Relation { X = NativeNumberTheory.Mod(a * x + b, n), Square = 1, Factors = factors.ToArray() };
                    lock (relationOwner.relationLock) divisor = relationOwner.CollectRelation(relation, value);
                    if (divisor > 1) return divisor;
                }
            }
        }

        private BigInteger CollectRelation(Relation relation, BigInteger cofactor)
        {
            if (!cofactor.IsOne)
            {
                if (!partials.TryGetValue(cofactor, out var previous))
                {
                    // Keep memory bounded even for difficult inputs.
                    if (partials.Count >= 100000) partials.Clear();
                    partials[cofactor] = relation;
                    return BigInteger.One;
                }
                partials.Remove(cofactor);
                var combined = new int[relation.Factors.Length + previous.Factors.Length];
                Array.Copy(relation.Factors, combined, relation.Factors.Length);
                Array.Copy(previous.Factors, 0, combined, relation.Factors.Length, previous.Factors.Length);
                relation.Factors = combined;
                relation.X = relation.X * previous.X % n;
                relation.Square = relation.Square * previous.Square * cofactor % n;
            }
            return AddRelation(relation);
        }

        private void PreparePrimePowers()
        {
            for (int index = 0; index < factorBase.Count; index++)
            {
                token.ThrowIfCancellationRequested();
                var prime = factorBase[index];
                int p = prime.P;
                var roots = new List<int> { prime.SquareRoot };
                if (p != 2 && prime.SquareRoot != 0) roots.Add(p - prime.SquareRoot);
                int modulus = p;
                while (true)
                {
                    powers.Add(new SievePower { Modulus = modulus, PrimeIndex = index,
                        SquareRoots = roots.ToArray(), Log = prime.Log });
                    if (modulus > scores.Length / p) break;
                    int nextModulus = modulus * p;
                    var lifted = new List<int>();
                    foreach (int r in roots)
                    {
                        if (p == 2)
                        {
                            for (int next = r; next < nextModulus; next += modulus)
                                if (((long)next * next - sievedNumber) % nextModulus == 0) lifted.Add(next);
                        }
                        else
                        {
                            int derivative = 2 * r % p;
                            if (derivative == 0) continue;
                            var quotient = (sievedNumber - (long)r * r) / modulus;
                            int adjustment = (int)NativeNumberTheory.Mod(quotient * Inverse(derivative, p), p);
                            lifted.Add(r + modulus * adjustment);
                        }
                    }
                    if (lifted.Count == 0) break;
                    roots = lifted;
                    modulus = nextModulus;
                }
            }
        }

        private void StartFamily()
        {
            double logTarget = BigInteger.Log(targetA);
            int count = Math.Max(2, Math.Min(12, (int)Math.Round(logTarget / Math.Log(4096))));
            double desired = Math.Exp(Math.Min(logTarget / count, Math.Log(factorBase[factorBase.Count - 1].P)));
            int first = LowerBound((int)(desired / 2)), last = LowerBound((int)(desired * 2));
            int minimum = LowerBound(11);
            first = Math.Max(minimum, first);
            last = Math.Max(first + count + 1, Math.Min(factorBase.Count, last));
            last = Math.Min(factorBase.Count, last);
            aFactors = new int[count];
            a = 1;
            for (int i = 0; i < count; i++)
            {
                token.ThrowIfCancellationRequested();
                int index;
                if (i == count - 1)
                {
                    var wanted = targetA / a;
                    index = LowerBound((int)BigInteger.Min(int.MaxValue, wanted));
                    index = Math.Max(minimum, Math.Min(factorBase.Count - 1, index));
                }
                else index = first + (int)(NextRandom() % (uint)(last - first));
                while (factorBase[index].SquareRoot == 0 || Array.IndexOf(aFactors, index + 1, 0, i) >= 0)
                {
                    index++;
                    if (index == factorBase.Count) index = minimum;
                }
                aFactors[i] = index + 1;
                a *= factorBase[index].P;
            }
            crtParts = new BigInteger[count];
            b = 0;
            for (int i = 0; i < count; i++)
            {
                var prime = factorBase[aFactors[i] - 1];
                var quotient = a / prime.P;
                int gamma = (int)((long)prime.SquareRoot * Inverse((int)(quotient % prime.P), prime.P) % prime.P);
                if (gamma > prime.P / 2) gamma = prime.P - gamma;
                crtParts[i] = quotient * gamma;
                b += crtParts[i];
            }
            Array.Clear(presieve, 0, presieve.Length);
            foreach (var power in powers)
            {
                token.ThrowIfCancellationRequested();
                int modulus = power.Modulus;
                power.InverseA = Inverse((int)(a % modulus), modulus);
                power.Roots = new int[power.InverseA == 0 ? 1 : power.SquareRoots.Length];
                if (power.InverseA == 0) continue;
                int bmod = (int)NativeNumberTheory.Mod(b, modulus);
                for (int i = 0; i < power.Roots.Length; i++)
                    power.Roots[i] = (int)(((long)(power.SquareRoots[i] - bmod + modulus) * power.InverseA + radius) % modulus);
                power.Shifts = new int[count - 1];
                for (int i = 0; i < power.Shifts.Length; i++)
                    power.Shifts[i] = (int)(2 * crtParts[i] % modulus * power.InverseA % modulus);
                if (PresieveWheel % modulus == 0)
                    foreach (int root in power.Roots)
                        for (int i = root; i < presieve.Length; i += modulus) presieve[i] += power.Log;
            }
            int inverseWheel = Inverse((int)(a % PresieveWheel), PresieveWheel);
            presieveShifts = new int[count - 1];
            for (int i = 0; i < presieveShifts.Length; i++)
                presieveShifts[i] = (int)(2 * crtParts[i] % PresieveWheel * inverseWheel % PresieveWheel);
            presieveOffset = 0;
            polynomialIndex = 0;
            familySize = 1 << (count - 1);
        }

        private void SievePolynomial(BigInteger c, int changedPart, int direction)
        {
            if (changedPart >= 0)
            {
                presieveOffset -= direction * presieveShifts[changedPart];
                if (presieveOffset < 0) presieveOffset += PresieveWheel;
                else if (presieveOffset >= PresieveWheel) presieveOffset -= PresieveWheel;
            }
            int offset = presieveOffset;
            for (int index = 0; index < scores.Length;)
            {
                int length = Math.Min(PresieveWheel - offset, scores.Length - index);
                Array.Copy(presieve, offset, scores, index, length);
                index += length; offset = 0;
            }
            foreach (var power in powers)
            {
                token.ThrowIfCancellationRequested();
                int modulus = power.Modulus;
                var prime = factorBase[power.PrimeIndex];
                if (power.InverseA == 0)
                {
                    int p = prime.P;
                    int inverse = Inverse((int)NativeNumberTheory.Mod(2 * b, p), p);
                    int r = (int)NativeNumberTheory.Mod(-c * inverse, p);
                    for (int m = p; m < modulus; m *= p)
                        r += m * (int)NativeNumberTheory.Mod(-((a * r + 2 * b) * r + c) / m * inverse, p);
                    power.Roots[0] = (r + radius) % modulus;
                }
                else if (changedPart >= 0)
                {
                    int shift = direction * power.Shifts[changedPart];
                    for (int i = 0; i < power.Roots.Length; i++)
                    {
                        int root = power.Roots[i] + shift;
                        if (root < 0) root += modulus;
                        else if (root >= modulus) root -= modulus;
                        power.Roots[i] = root;
                    }
                }
                if (modulus == prime.P)
                {
                    prime.Root1 = power.Roots[0];
                    prime.Root2 = power.Roots[power.Roots.Length - 1];
                }
                if (PresieveWheel % modulus != 0)
                    foreach (int root in power.Roots)
                        for (int i = root; i < scores.Length; i += modulus) scores[i] += power.Log;
            }
        }

        private int LowerBound(int value)
        {
            int lo = 0, hi = factorBase.Count;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (factorBase[mid].P < value) lo = mid + 1; else hi = mid;
            }
            return lo;
        }

        private uint NextRandom()
        {
            randomState ^= randomState << 13;
            randomState ^= randomState >> 17;
            randomState ^= randomState << 5;
            return randomState;
        }

        private static int Inverse(int value, int modulus)
        {
            int r = modulus, next = value;
            long t = 0, coefficient = 1;
            while (next != 0)
            {
                int quotient = r / next;
                int remainder = r - quotient * next;
                r = next; next = remainder;
                long previous = t;
                t = coefficient; coefficient = previous - quotient * coefficient;
            }
            return r == 1 ? (int)((t % modulus + modulus) % modulus) : 0;
        }

        private BigInteger AddRelation(Relation relation)
        {
            var row = new Row { Parity = new ulong[words], Combination = new ulong[words] };
            foreach (int index in relation.Factors) row.Parity[index / 64] ^= 1UL << (index % 64);
            int id = relations.Count;
            row.Combination[id / 64] = 1UL << (id % 64);
            for (int pivot = basis.Length - 1; pivot >= 0; pivot--)
            {
                if ((pivot & 255) == 0) token.ThrowIfCancellationRequested();
                if ((row.Parity[pivot / 64] & (1UL << (pivot % 64))) == 0) continue;
                var previous = basis[pivot];
                if (previous == null)
                {
                    basis[pivot] = row;
                    relations.Add(relation);
                    return BigInteger.One;
                }
                for (int i = 0; i < words; i++)
                {
                    row.Parity[i] ^= previous.Parity[i];
                    row.Combination[i] ^= previous.Combination[i];
                }
            }
            var exponents = new int[basis.Length];
            BigInteger x = 1, y = 1;
            for (int i = 0; i <= id; i++)
            {
                token.ThrowIfCancellationRequested();
                if ((row.Combination[i / 64] & (1UL << (i % 64))) == 0) continue;
                var selected = i == id ? relation : relations[i];
                x = x * selected.X % n;
                y = y * selected.Square % n;
                foreach (int factor in selected.Factors) exponents[factor]++;
            }
            for (int i = 1; i < exponents.Length; i++)
                if (exponents[i] != 0) y = y * BigInteger.ModPow(factorBase[i - 1].P, exponents[i] / 2, n) % n;
            // Both square roots are tried; even a trivial dependency is harmless.
            var divisor = BigInteger.GreatestCommonDivisor(x - y, n);
            if (divisor > 1 && divisor < n) return divisor;
            divisor = BigInteger.GreatestCommonDivisor(x + y, n);
            return divisor > 1 && divisor < n ? divisor : BigInteger.One;
        }

        private static int SquareRootModPrime(int value, int p)
        {
            if (p == 2 || value == 0) return value;
            if (p % 4 == 3) return InternalMath.ModPow(value, (p + 1) / 4, p);
            int odd = p - 1, twos = 0;
            while ((odd & 1) == 0) { odd >>= 1; twos++; }
            int nonresidue = 2;
            while (InternalMath.Legendre(nonresidue, p) != -1) nonresidue++;
            long c = InternalMath.ModPow(nonresidue, odd, p);
            long x = InternalMath.ModPow(value, (odd + 1) / 2, p);
            long t = InternalMath.ModPow(value, odd, p);
            while (t != 1)
            {
                int i = 0;
                long power = t;
                while (power != 1) { power = power * power % p; i++; }
                long b = c;
                for (int j = 0; j < twos - i - 1; j++) b = b * b % p;
                x = x * b % p;
                c = b * b % p;
                t = t * c % p;
                twos = i;
            }
            return (int)x;
        }
    }
}

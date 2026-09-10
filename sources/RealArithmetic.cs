using System;
using System.Numerics;
using System.Threading;

namespace EllipticCurves
{
    // Fixed-grid rational interval arithmetic. Every operation rounds outwards.
    // Log uses 2 atanh(z) with an explicit geometric tail; pi uses Machin's formula.
    internal sealed class RealArithmetic
    {
        internal readonly BigInteger Scale;
        internal readonly RealComputationOptions Options;
        internal readonly CancellationToken Token;
        internal readonly BigRational Tolerance;
        private RealEnclosure logTwo, pi;
        internal RealArithmetic(RealComputationOptions options, CancellationToken token)
        { Options = options; Token = token; Scale = BigInteger.One << options.PrecisionBits; Tolerance = new BigRational(1, BigInteger.Pow(10, options.DecimalDigits)); }
        internal static int BitLength(BigInteger n)
        { var b = BigInteger.Abs(n).ToByteArray(); int k = b.Length - 1; while (k > 0 && b[k] == 0) k--; int count = 8 * k; for (int v = b[k]; v > 0; v >>= 1) count++; return count; }
        internal RealEnclosure I(BigRational r) => new RealEnclosure(r, r);
        internal RealEnclosure Round(BigRational l, BigRational u) => new RealEnclosure(new BigRational(DescentPolynomial.Floor(l * new BigRational(Scale)), Scale), new BigRational(DescentPolynomial.Ceiling(u * new BigRational(Scale)), Scale));
        internal RealEnclosure Add(RealEnclosure a, RealEnclosure b) => Round(a.LowerBound + b.LowerBound, a.UpperBound + b.UpperBound);
        internal RealEnclosure Neg(RealEnclosure a) => new RealEnclosure(-a.UpperBound, -a.LowerBound);
        internal RealEnclosure Sub(RealEnclosure a, RealEnclosure b) => Add(a, Neg(b));
        internal RealEnclosure Mul(RealEnclosure a, RealEnclosure b)
        {
            var x = a.LowerBound * b.LowerBound; var y = a.LowerBound * b.UpperBound; var z = a.UpperBound * b.LowerBound; var w = a.UpperBound * b.UpperBound;
            return Round(Min(Min(x, y), Min(z, w)), Max(Max(x, y), Max(z, w)));
        }
        internal RealEnclosure Div(RealEnclosure a, RealEnclosure b)
        {
            if (b.Contains(0)) throw new ArithmeticException("Insufficient precision: denominator enclosure contains zero.");
            return Mul(a, new RealEnclosure(1 / b.UpperBound, 1 / b.LowerBound));
        }
        internal RealEnclosure Abs(RealEnclosure a) => a.LowerBound >= 0 ? a : a.UpperBound <= 0 ? Neg(a) : new RealEnclosure(0, Max(-a.LowerBound, a.UpperBound));
        internal RealEnclosure Sqrt(RealEnclosure a)
        {
            if (a.LowerBound < 0) throw new ArithmeticException("Insufficient precision to enclose a real square root.");
            BigInteger Root(BigRational x) => InternalMath.IntegerSqrt(DescentPolynomial.Floor(x * new BigRational(Scale * Scale)));
            var lo = Root(a.LowerBound); var hi = Root(a.UpperBound);
            if (new BigRational(hi * hi, Scale * Scale) < a.UpperBound) hi++;
            return new RealEnclosure(new BigRational(lo, Scale), new BigRational(hi, Scale));
        }
        internal RealEnclosure Log(RealEnclosure a)
        {
            if (a.LowerBound <= 0) throw new ArithmeticException("Insufficient precision to enclose a logarithm.");
            return new RealEnclosure(LogPoint(a.LowerBound).LowerBound, LogPoint(a.UpperBound).UpperBound);
        }
        private RealEnclosure UnitLog(BigRational q)
        {
            var z = I((q - 1) / (q + 1)); var z2 = Mul(z, z); var power = z; var sum = I(0);
            int count = Options.PrecisionBits / 3 + 12;
            for (int k = 0; k < count; k++)
            { Token.ThrowIfCancellationRequested(); sum = Add(sum, Div(power, I(2 * k + 1))); power = Mul(power, z2); }
            var tail = Div(Mul(I(2), Abs(power)), Mul(I(2 * count + 1), Sub(I(1), z2))).UpperBound;
            return Add(Mul(I(2), sum), new RealEnclosure(0, tail));
        }
        private RealEnclosure LogPoint(BigRational q)
        {
            int exponent = BitLength(q.Num) - BitLength(q.Den);
            var factor = exponent >= 0 ? new BigRational(BigInteger.One << exponent) : new BigRational(1, BigInteger.One << -exponent);
            q /= factor;
            if (q < 1) { q *= 2; exponent--; }
            if (q >= 2) { q /= 2; exponent++; }
            logTwo ??= UnitLog(2);
            return Add(UnitLog(q), Mul(I(exponent), logTwo));
        }
        internal RealEnclosure Pi()
        {
            if (pi != null) return pi;
            RealEnclosure Atan(int denominator)
            {
                var z = I(new BigRational(1, denominator)); var z2 = Mul(z, z); var power = z; var sum = I(0);
                int count = Options.PrecisionBits / 4 + 12;
                for (int k = 0; k < count; k++)
                { Token.ThrowIfCancellationRequested(); var term = Div(power, I(2 * k + 1)); sum = k % 2 == 0 ? Add(sum, term) : Sub(sum, term); power = Mul(power, z2); }
                var tail = Div(Abs(power), I(2 * count + 1)).UpperBound;
                return Add(sum, new RealEnclosure(-tail, tail));
            }
            pi = Sub(Mul(I(16), Atan(5)), Mul(I(4), Atan(239))); return pi;
        }
        internal RealEnclosure Agm(RealEnclosure a, RealEnclosure b)
        {
            if (a.LowerBound <= 0 || b.LowerBound <= 0) throw new ArithmeticException("AGM requires positive enclosures.");
            for (int k = 0; k < Options.MaxIterations; k++)
            {
                Token.ThrowIfCancellationRequested();
                var bounds = new RealEnclosure(Min(a.LowerBound, b.LowerBound), Max(a.UpperBound, b.UpperBound));
                // Resolve the AGM to working precision: an early absolute-width
                // stop can be amplified arbitrarily by the reciprocal period formula.
                if (bounds.Width < Min(Tolerance / 100, new BigRational(1, BigInteger.One << (Options.PrecisionBits / 2)))) return bounds;
                var nextA = Div(Add(a, b), I(2)); b = Sqrt(Mul(a, b)); a = nextA;
            }
            throw new ArithmeticException("AGM precision or iteration limit reached.");
        }
        internal RealEnclosure Finish(RealEnclosure a)
        {
            if (a.Width > Tolerance) throw new ArithmeticException("Requested accuracy was not reached; increase PrecisionBits and, for a large regulator, GuardDigits.");
            return a;
        }
        internal static BigRational Min(BigRational a, BigRational b) => a < b ? a : b;
        internal static BigRational Max(BigRational a, BigRational b) => a > b ? a : b;
    }
}

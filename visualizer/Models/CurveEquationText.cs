#nullable enable
using System.Numerics;
using System.Text;

namespace EllipticCurves.Visualizer.Models;

/// <summary>A bounded polynomial parser for rational Weierstrass equations.</summary>
public static class CurveEquationText
{
    // Text length and nesting alone do not bound expressions such as ((2^3)^3)^3.
    // Limit intermediate rationals as well, since parsing runs on the editor thread.
    private const int MaxCoefficientBits = 32768;

    public static bool TryParse(string text, out EllipticCurveQ? curve, out string error)
    {
        curve = null;
        error = "";
        if (string.IsNullOrWhiteSpace(text)) { error = "Enter an equation, such as y^2 = x^3 - x."; return false; }
        if (text.Length > 4096) { error = "Keep the equation within 4096 characters."; return false; }
        if (text.Contains('²') || text.Contains('³')) { error = "Enter powers using ^, for example y^2 = x^3 - x."; return false; }
        try
        {
            var normalized = text.Replace('−', '-').Replace('–', '-').Replace('·', '*').Replace('×', '*').ToLowerInvariant();
            var sides = normalized.Split('=');
            if (sides.Length != 2) throw new FormatException("Use one equals sign, for example y^2 = x^3 - x.");
            var polynomial = Add(new Parser(sides[0]).Parse(), new Parser(sides[1]).Parse(), -1);
            var allowed = new[] { (0, 2), (3, 0), (1, 1), (2, 0), (0, 1), (1, 0), (0, 0) };
            if (polynomial.Keys.Any(term => !allowed.Contains(term)))
                throw new FormatException("Use Weierstrass terms: y^2, x^3, xy, x^2, y, x and a constant.");
            var leading = Coefficient(polynomial, (0, 2));
            if (leading.IsZero || Coefficient(polynomial, (3, 0)) != -leading)
                throw new FormatException("The equation must reduce to y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6.");
            curve = new EllipticCurveQ(Bounded(Coefficient(polynomial, (1, 1)) / leading),
                Bounded(-Coefficient(polynomial, (2, 0)) / leading), Bounded(Coefficient(polynomial, (0, 1)) / leading),
                Bounded(-Coefficient(polynomial, (1, 0)) / leading), Bounded(-Coefficient(polynomial, (0, 0)) / leading));
            return true;
        }
        catch (FormatException exception) { error = exception.Message; return false; }
    }

    public static string Format(EllipticCurveQ curve)
    {
        var text = new StringBuilder("y^2");
        Append(text, curve.A1, "x*y");
        Append(text, curve.A3, "y");
        text.Append(" = x^3");
        Append(text, curve.A2, "x^2");
        Append(text, curve.A4, "x");
        Append(text, curve.A6, "");
        return text.ToString();
    }

    private static void Append(StringBuilder text, BigRational coefficient, string term)
    {
        if (coefficient.IsZero) return;
        text.Append(coefficient.Sign < 0 ? " - " : " + ");
        var magnitude = coefficient.Abs();
        if (term.Length == 0 || magnitude != BigRational.One)
        {
            text.Append(RationalText.Format(magnitude));
            if (term.Length > 0) text.Append('*');
        }
        text.Append(term);
    }

    private static BigRational Coefficient(Dictionary<(int X, int Y), BigRational> polynomial, (int, int) term)
        => polynomial.TryGetValue(term, out var value) ? value : BigRational.Zero;

    private static BigRational Bounded(BigRational value)
    {
        if (BigInteger.Abs(value.Num).GetBitLength() > MaxCoefficientBits || value.Den.GetBitLength() > MaxCoefficientBits)
            throw new FormatException("The equation produces numbers that are too large. Simplify the coefficients.");
        return value;
    }

    private static Dictionary<(int X, int Y), BigRational> Add(Dictionary<(int X, int Y), BigRational> left,
        Dictionary<(int X, int Y), BigRational> right, int sign)
    {
        var result = new Dictionary<(int X, int Y), BigRational>(left);
        foreach (var (term, value) in right)
        {
            var sum = Bounded(Coefficient(result, term) + sign * value);
            if (sum.IsZero) result.Remove(term);
            else result[term] = sum;
        }
        return result;
    }

    private static Dictionary<(int X, int Y), BigRational> Multiply(Dictionary<(int X, int Y), BigRational> left,
        Dictionary<(int X, int Y), BigRational> right)
    {
        var result = new Dictionary<(int X, int Y), BigRational>();
        foreach (var (a, av) in left)
            foreach (var (b, bv) in right)
            {
                var term = (X: a.X + b.X, Y: a.Y + b.Y);
                if (term.X + term.Y > 3) throw new FormatException("Only polynomials of degree at most 3 are supported.");
                var sum = Bounded(Coefficient(result, term) + Bounded(av * bv));
                if (sum.IsZero) result.Remove(term);
                else result[term] = sum;
            }
        return result;
    }

    private static Dictionary<(int X, int Y), BigRational> Constant(BigRational value)
        => value.IsZero ? new() : new() { [(0, 0)] = Bounded(value) };

    private sealed class Parser(string input)
    {
        private int position, depth;

        public Dictionary<(int X, int Y), BigRational> Parse()
        {
            var value = Sum();
            if (Peek() != '\0') throw new FormatException($"Unexpected character '{Peek()}'. Use x, y, numbers and arithmetic operators.");
            return value;
        }

        private char Peek()
        {
            while (position < input.Length && char.IsWhiteSpace(input[position])) position++;
            return position < input.Length ? input[position] : '\0';
        }

        private Dictionary<(int X, int Y), BigRational> Sum()
        {
            var value = Product();
            while (Peek() is '+' or '-')
            {
                var sign = input[position++] == '+' ? 1 : -1;
                value = Add(value, Product(), sign);
            }
            return value;
        }

        private Dictionary<(int X, int Y), BigRational> Product()
        {
            var value = Unary();
            while (true)
            {
                var next = Peek();
                if (next == '*') { position++; value = Multiply(value, Unary()); }
                else if (next == '/')
                {
                    position++;
                    var divisor = Unary();
                    if (divisor.Keys.Any(term => term != (0, 0))) throw new FormatException("Divide only by a nonzero numeric constant.");
                    var constant = Coefficient(divisor, (0, 0));
                    if (constant.IsZero) throw new FormatException("A denominator cannot be zero.");
                    value = Multiply(value, Constant(constant.Reciprocal()));
                }
                else if (next is 'x' or 'y' or '(') value = Multiply(value, Unary());
                else return value;
            }
        }

        private Dictionary<(int X, int Y), BigRational> Unary()
        {
            if (++depth > 64) throw new FormatException("The equation is nested too deeply.");
            try
            {
                if (Peek() is '+' or '-')
                {
                    var sign = input[position++] == '+' ? 1 : -1;
                    return Multiply(Constant(sign), Unary());
                }
                var value = Atom();
                if (Peek() == '^')
                {
                    position++;
                    var digit = Peek();
                    if (digit is < '0' or > '3') throw new FormatException("Use integer powers from 0 to 3.");
                    position++;
                    var power = digit - '0';
                    var result = Constant(BigRational.One);
                    for (var i = 0; i < power; i++) result = Multiply(result, value);
                    value = result;
                }
                return value;
            }
            finally { depth--; }
        }

        private Dictionary<(int X, int Y), BigRational> Atom()
        {
            var next = Peek();
            if (next == '(')
            {
                position++;
                var value = Sum();
                if (Peek() != ')') throw new FormatException("Close the parenthesis to finish the expression.");
                position++;
                return value;
            }
            if (next is 'x' or 'y')
            {
                position++;
                return new() { [next == 'x' ? (1, 0) : (0, 1)] = BigRational.One };
            }
            var start = position;
            while (position < input.Length && (char.IsAsciiDigit(input[position]) || input[position] is '.' or ',')) position++;
            if (position > start && position < input.Length && input[position] == 'e')
            {
                position++;
                if (position < input.Length && input[position] is '+' or '-') position++;
                while (position < input.Length && char.IsAsciiDigit(input[position])) position++;
            }
            if (position == start || !RationalText.TryParse(input[start..position], out var number))
                throw new FormatException("Finish the term with x, y or a number (for example 8.325, 2/7 or 1e-5).");
            return Constant(number);
        }
    }
}

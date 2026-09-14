using System.Diagnostics;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;

// Exact specialization and certification of independently constructed bisections.
static class BisectionSearch
{
    sealed record Bisection(int Index, BigRational[] H, BigRational[] L, BigRational[] A,
        BigRational[] C, BigRational[] U, BigRational[] S);
    static BigRational Q(string s)
    {
        var p = s.Split('/');
        if (p.Length > 2) throw new FormatException("Expected a rational number");
        return new(BigInteger.Parse(p[0]), p.Length == 1 ? BigInteger.One : BigInteger.Parse(p[1]));
    }
    static BigRational Eval(BigRational[] a, BigRational t)
    {
        BigRational result = 0;
        for (int i = a.Length - 1; i >= 0; i--) result = result*t + a[i];
        return result;
    }
    static BigInteger Sqrt(BigInteger n)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
        if (n < 2) return n;
        var x = BigInteger.One << checked((int)((n.GetBitLength()+1)/2));
        while (true) { var y = (x+n/x)/2; if (y >= x) return x; x = y; }
    }
    static bool Square(BigRational q, out BigRational root)
    {
        root = 0;
        if (q.Sign < 0) return false;
        var n = Sqrt(q.Num); if (n*n != q.Num) return false;
        var d = Sqrt(q.Den); if (d*d != q.Den) return false;
        root = new(n, d); return true;
    }
    static readonly bool[] SquareMod64 = Enumerable.Range(0,64).Select(i=>Enumerable.Range(0,64).Any(j=>j*j%64==i)).ToArray();
    static EllipticCurvePoint? Extra(Bisection b, BigRational t, BigRational lFamily, BigRational bFamily,
        BigInteger n2, BigInteger nd, BigInteger d2)
    {
        BigRational r;
        if (b.H.Length==3 && b.H.All(c=>c.Den.IsOne))
        {
            var homogeneous = b.H[2].Num*n2+b.H[1].Num*nd+b.H[0].Num*d2;
            if (homogeneous.Sign<0 || !SquareMod64[(int)(homogeneous & 63)]) return null;
            var root = Sqrt(homogeneous); if (root*root != homogeneous) return null;
            r = new(root,t.Den);
        }
        else if (!Square(Eval(b.H,t), out r)) return null;
        var l = Eval(b.L,t); if (l.IsZero) return null;
        var x = (-Eval(b.U,t)+Eval(b.S,t)*r)/2;
        var y = (Eval(b.A,t)*x-Eval(b.C,t))/l + (lFamily*x+bFamily)/2;
        return new(x * new BigRational(BigInteger.Pow(t.Den,4)), y * new BigRational(BigInteger.Pow(t.Den,6)));
    }

    public static void Certify(string input, string output, int limit, CancellationToken token)
    {
        using var bis = JsonDocument.Parse(File.ReadAllText(Path.Combine(input,"bisections.json")));
        var all = bis.RootElement.GetProperty("bisections").EnumerateArray().Select(b => {
            BigRational[] Coeff(string name) => b.GetProperty(name).EnumerateArray().Select(c=>Q(c.GetString()!)).ToArray();
            return new Bisection(b.GetProperty("index").GetInt32(),Coeff("h"),Coeff("l"),Coeff("a"),Coeff("c"),Coeff("quadratic_u"),Coeff("sqrt_factor"));
        }).ToArray();
        var unique = all.GroupBy(b=>string.Join(",",b.H)).Select(g=>g.First()).ToArray();
        using var advanced = JsonDocument.Parse(File.ReadAllText(Path.Combine(input,"advanced.json")));
        Directory.CreateDirectory(output);
        var watch = Stopwatch.StartNew(); var rows = new List<object>(); int best = 0, index = 0;
        foreach (var candidate in advanced.RootElement.GetProperty("candidates").EnumerateArray().Take(limit))
        {
            token.ThrowIfCancellationRequested();
            var timer = Stopwatch.StartNew();
            var p = candidate.GetProperty("parameter");
            var t = new BigRational(BigInteger.Parse(p[0].GetString()!),BigInteger.Parse(p[1].GetString()!));
            var (curve, seeds) = Structured302.Create(t.Num,t.Den);
            var f = Family302.Parameters(t.Num,t.Den);
            var l = new BigRational(f[0],BigInteger.Pow(t.Den,2));
            var b = new BigRational(f[5],BigInteger.Pow(t.Den,6));
            var points = seeds.ToList(); var onConics = new List<int>();
            var n2=t.Num*t.Num;var nd=t.Num*t.Den;var d2=t.Den*t.Den;
            foreach (var section in unique)
            {
                var extra = Extra(section,t,l,b,n2,nd,d2);
                if (extra is null) continue;
                if (!curve.IsOnCurve(extra.Value)) throw new InvalidOperationException("Bisection point is off the curve");
                points.Add(extra.Value); onConics.Add(section.Index);
            }
            var cert = curve.GetRankLowerBound(points,1009,token);
            best = Math.Max(best,cert.LowerBound);
            string name = $"curve_{index:D4}.json";
            Hunt.Save(Path.Combine(output,name),new {
                ainvs = new[]{curve.A1,curve.A2,curve.A3,curve.A4,curve.A6}.Select(a=>a.ToString()).ToArray(),
                points = points.Select(q=>new[]{q.X.ToString(),q.Y.ToString()}).ToArray(),
                parameter = new[]{t.Num.ToString(),t.Den.ToString()},
                source = candidate.Clone(), bisection_indices = onConics,
                rank_lower_bound = cert.LowerBound, cert.ImageDimension, cert.NoTwoTorsionPrime,
                hypotheses = Array.Empty<string>() });
            rows.Add(new { index, file = name, lower_bound = cert.LowerBound, point_count = points.Count,
                conics = onConics, parameter_digits = Math.Max(t.Num.ToString().TrimStart('-').Length,t.Den.ToString().Length),
                seconds = timer.Elapsed.TotalSeconds });
            if (index%25==0 || cert.LowerBound>19)
                Console.WriteLine($"{index}: bound {cert.LowerBound}, {points.Count} points, {onConics.Count} conics, {timer.Elapsed.TotalSeconds:F3}s");
            Hunt.Save(Path.Combine(output,"summary.json"),new { input, best_lower_bound = best,
                count = rows.Count, seconds = watch.Elapsed.TotalSeconds, rows });
            index++;
        }
        Console.WriteLine($"Certified {rows.Count} curves in {watch.Elapsed.TotalSeconds:F3}s; best lower bound {best}");
    }
}

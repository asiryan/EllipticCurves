using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using EllipticCurves;

// Run from any directory. Every invocation recomputes the result and proves all
// prime factors; no factor hints or cached factorizations enter the measurement.
CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
Console.WriteLine("case,digits,first_ms,median_ms");
void Measure(string label, int digits, Action<CancellationToken> calculation)
{
    double Run()
    {
        using var timeout = new CancellationTokenSource(TimeSpan.FromMinutes(2));
        var watch = Stopwatch.StartNew();
        calculation(timeout.Token);
        return watch.Elapsed.TotalMilliseconds;
    }
    double first = Run();
    var samples = new[] { Run(), Run(), Run() };
    Array.Sort(samples);
    Console.WriteLine($"{label},{digits},{first:F1},{samples[1]:F1}");
}

var curve = new EllipticCurveQ(0, 1, 0,
    new BigRational(BigInteger.Parse("-221556180740323405132844117936")),
    new BigRational(BigInteger.Parse("35386140191724122461245294467670188433973860")));
var conductor = BigInteger.Parse("1103561624055499058867562340698878392772504928025988266715523317532246643920");
Measure("conductor", 90, token =>
{
    if (curve.GetConductor(token) != conductor) throw new InvalidOperationException("Conductor mismatch.");
});

int index = 0;
foreach (var line in File.ReadLines(Path.Combine(AppContext.BaseDirectory, "factorization.csv")).Skip(1))
{
    var fields = line.Split(',');
    var p = BigInteger.Parse(fields[0]);
    var q = BigInteger.Parse(fields[1]);
    var n = p * q;
    Measure((++index).ToString(), n.ToString().Length, token =>
    {
        var factors = NativeNumberTheory.Factor(n, token);
        if (factors.Count != 2 || factors[p] != 1 || factors[q] != 1)
            throw new InvalidOperationException("Factorization mismatch.");
    });
}

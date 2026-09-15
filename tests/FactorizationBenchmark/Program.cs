using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using EllipticCurves;

// Run from any directory. Every invocation recomputes the result and proves all
// prime factors; no factor hints or cached factorizations enter the measurement.
CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
if (args.Length > 0 && args[0] == "--conductor-worker")
{
    await ConductorWorkerBenchmark.RunAsync(args);
    return;
}
if (args.Length > 0 && args[0] == "--conductor-cpu")
{
    var largeCurve = new EllipticCurveQ(1, 0, 0,
        new BigRational(BigInteger.Parse("-20820207864197471248300179976626")),
        new BigRational(BigInteger.Parse("36732936589138673862895758597955508398047757956")));
    var expected = BigInteger.Parse("21785392458764315483988614758901932764423833726404496768609056369259382575196330");
    var limits = args.Length > 1 ? args.Skip(1).Select(int.Parse).ToArray() : new[] { 4, 8, 12, 0 };
    Console.WriteLine("worker_limit,logical_cpus,elapsed_ms,cpu_ms");
    foreach (int workers in limits)
    {
        using var timeout = new CancellationTokenSource(TimeSpan.FromMinutes(5));
        using var process = Process.GetCurrentProcess();
        var cpuStarted = process.TotalProcessorTime;
        var watch = Stopwatch.StartNew();
        var result = largeCurve.GetConductor(new FactorizationOptions { MaxDegreeOfParallelism = workers }, timeout.Token);
        if (result != expected) throw new InvalidOperationException("Conductor mismatch.");
        Console.WriteLine($"{workers},{Environment.ProcessorCount},{watch.Elapsed.TotalMilliseconds:F1},{(process.TotalProcessorTime - cpuStarted).TotalMilliseconds:F1}");
    }
    return;
}
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

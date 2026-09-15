using System.Diagnostics;
using System.Text.Json;

internal static class ConductorWorkerBenchmark
{
    internal static async Task RunAsync(string[] args)
    {
        if (args.Length < 2) throw new ArgumentException("Supply the Explorer executable path, then optional worker and sample counts.");
        string executable = Path.GetFullPath(args[1]);
        int workers = args.Length > 2 ? int.Parse(args[2]) : 12;
        int samples = args.Length > 3 ? int.Parse(args[3]) : 1;
        const string expected = "21785392458764315483988614758901932764423833726404496768609056369259382575196330";
        var request = new
        {
            OperationId = "EllipticCurveQ.GetConductor(CancellationToken)",
            Equation = "y^2 + xy = x^3 - 20820207864197471248300179976626*x + 36732936589138673862895758597955508398047757956",
            Arguments = new Dictionary<string, string> { ["options.MaxDegreeOfParallelism"] = workers.ToString() },
            TimeoutSeconds = 300, MaxItems = 1000
        };
        Console.WriteLine("sample,worker_limit,startup_ms,calculation_ms,format_and_exit_ms,total_ms,cpu_ms");
        for (int sample = 1; sample <= samples; sample++)
        {
            using var timeout = new CancellationTokenSource(TimeSpan.FromMinutes(5));
            var start = new ProcessStartInfo(executable)
            {
                UseShellExecute = false, CreateNoWindow = true,
                RedirectStandardInput = true, RedirectStandardOutput = true, RedirectStandardError = true
            };
            start.ArgumentList.Add("--compute-worker");
            using var process = new Process { StartInfo = start };
            var watch = Stopwatch.StartNew();
            if (!process.Start()) throw new InvalidOperationException("Cannot start Explorer worker.");
            try
            {
                var errors = process.StandardError.ReadToEndAsync(timeout.Token);
                await process.StandardInput.WriteLineAsync(JsonSerializer.Serialize(request).AsMemory(), timeout.Token);
                await process.StandardInput.FlushAsync(timeout.Token);
                double started = -1, calculated = -1;
                string result = null;
                // Keep stdin open just as CalculationRunner does. Each sample is
                // a fresh Explorer process, including its runtime and JIT startup.
                while (await process.StandardOutput.ReadLineAsync(timeout.Token) is { } line)
                {
                    using var message = JsonDocument.Parse(line);
                    var root = message.RootElement;
                    string kind = root.GetProperty("Kind").GetString();
                    string text = root.GetProperty("Message").GetString();
                    if (kind == "error") throw new InvalidOperationException(text);
                    if (kind == "progress" && text.StartsWith("Computing ·")) started = watch.Elapsed.TotalMilliseconds;
                    if (kind == "progress" && text == "Collecting and formatting results") calculated = watch.Elapsed.TotalMilliseconds;
                    if (kind == "completed") result = root.GetProperty("Result").GetString();
                }
                await process.WaitForExitAsync(timeout.Token);
                double total = watch.Elapsed.TotalMilliseconds;
                string stderr = await errors;
                if (process.ExitCode != 0 || result == null || !result.Contains(expected) || started < 0 || calculated < started)
                    throw new InvalidOperationException("Worker result or timing markers are invalid. " + stderr);
                Console.WriteLine($"{sample},{workers},{started:F1},{calculated - started:F1},{total - calculated:F1},{total:F1},{process.TotalProcessorTime.TotalMilliseconds:F1}");
            }
            finally
            {
                if (!process.HasExited) { process.Kill(entireProcessTree: true); await process.WaitForExitAsync(); }
            }
        }
    }
}

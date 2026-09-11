#nullable enable
using System.IO;
using System.Text.Json;

namespace EllipticCurves.Visualizer.Computations;

public static class CalculationWorker
{
    public static async Task<int> RunAsync(TextReader input, TextWriter output)
    {
        var writeLock = new object();
        void Send(CalculationUpdate update)
        {
            lock (writeLock) { output.WriteLine(JsonSerializer.Serialize(update)); output.Flush(); }
        }
        try
        {
            var line = await input.ReadLineAsync().ConfigureAwait(false);
            if (line == null || line.Length > 2_000_000) throw new FormatException("Invalid calculation request.");
            var request = JsonSerializer.Deserialize<CalculationRequest>(line) ?? throw new FormatException("Missing calculation request.");
            var calculation = Task.Run(async () =>
            {
                try
                {
                    var result = await CalculationEngine.ExecuteAsync(request, Send).ConfigureAwait(false);
                    Send(new("completed", "Calculation complete", 100, result));
                    return 0;
                }
                catch (Exception error) { Send(new("error", error.Message)); return 1; }
            });
            // The host keeps stdin open. A disconnected host must not leave a long
            // non-cooperative library operation running as an orphan process.
            // Console.In is a synchronized reader: its ReadLineAsync may block
            // synchronously. Keep the disconnect monitor off the main thread.
            var disconnected = Task.Run(() => input.ReadLine());
            var finished = await Task.WhenAny(calculation, disconnected).ConfigureAwait(false);
            return finished == calculation ? await calculation.ConfigureAwait(false) : 2;
        }
        catch (Exception error) { Send(new("error", error.Message)); return 1; }
    }
}

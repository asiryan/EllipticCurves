#nullable enable
using System.IO;
using System.ComponentModel;
using System.Diagnostics;
using System.Text;
using System.Text.Json;

namespace EllipticCurves.Explorer.Computations;

public sealed record CalculationOutcome(string Status, string Message, string? Text = null);

/// <summary>Owns only the process it starts. Cancellation also stops methods without tokens.</summary>
public sealed class CalculationRunner(Func<ProcessStartInfo>? createStartInfo = null)
{
    private static ProcessStartInfo DefaultStartInfo()
    {
        var path = Environment.ProcessPath ?? throw new InvalidOperationException("Cannot locate the calculation executable.");
        var start = new ProcessStartInfo(path);
        if (Path.GetFileNameWithoutExtension(path).Equals("dotnet", StringComparison.OrdinalIgnoreCase))
            start.ArgumentList.Add(typeof(CalculationRunner).Assembly.Location);
        start.ArgumentList.Add(CalculationProtocol.WorkerArgument);
        return start;
    }

    public async Task<CalculationOutcome> RunAsync(CalculationRequest request, Action<CalculationUpdate> report, CancellationToken cancellationToken = default)
    {
        if (request.TimeoutSeconds is < 0 or > 86400) return new(CalculationStatus.Failed, "Time limit must be 0–86,400 seconds; 0 means unlimited.");
        using var timeout = new CancellationTokenSource();
        if (request.TimeoutSeconds > 0) timeout.CancelAfter(TimeSpan.FromSeconds(request.TimeoutSeconds));
        using var linked = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken, timeout.Token);
        using var process = new Process();
        var started = false;
        try
        {
            linked.Token.ThrowIfCancellationRequested();
            var start = (createStartInfo ?? DefaultStartInfo)();
            start.UseShellExecute = false;
            start.CreateNoWindow = true;
            start.RedirectStandardInput = start.RedirectStandardOutput = start.RedirectStandardError = true;
            start.StandardInputEncoding = start.StandardOutputEncoding = start.StandardErrorEncoding = new UTF8Encoding(false);
            process.StartInfo = start;
            started = process.Start();
            if (!started) return new(CalculationStatus.Failed, "The calculation process could not start.");
            using var stop = linked.Token.Register(() => Stop(process));
            // Drain stderr concurrently, so a startup diagnostic cannot block the child.
            var errors = process.StandardError.ReadToEndAsync(linked.Token);
            await process.StandardInput.WriteLineAsync(JsonSerializer.Serialize(request).AsMemory(), linked.Token).ConfigureAwait(false);
            await process.StandardInput.FlushAsync(linked.Token).ConfigureAwait(false);
            CalculationUpdate? final = null;
            while (await process.StandardOutput.ReadLineAsync(linked.Token).ConfigureAwait(false) is { } line)
            {
                var update = JsonSerializer.Deserialize<CalculationUpdate>(line) ?? throw new IOException("The calculation returned an empty response.");
                if (update.Kind is CalculationProtocol.Completed or CalculationProtocol.Error) final = update;
                report(update);
            }
            await process.WaitForExitAsync(linked.Token).ConfigureAwait(false);
            var stderr = await errors.ConfigureAwait(false);
            linked.Token.ThrowIfCancellationRequested();
            if (final?.Kind == CalculationProtocol.Completed && process.ExitCode == 0) return new(CalculationStatus.Completed, final.Message, final.Result);
            return new(CalculationStatus.Failed, final?.Message ?? ("The calculation process exited unexpectedly. " + stderr[..Math.Min(stderr.Length, 2000)]).Trim());
        }
        catch (Exception) when (linked.IsCancellationRequested)
        {
            return cancellationToken.IsCancellationRequested ? new(CalculationStatus.Cancelled, "Calculation stopped by you.")
                : new(CalculationStatus.TimedOut, "Stopped after " + request.TimeoutSeconds + " seconds. Increase the time limit to try again.");
        }
        catch (Exception error) when (error is IOException or InvalidOperationException or Win32Exception or JsonException or UnauthorizedAccessException)
        { return new(CalculationStatus.Failed, error.Message); }
        finally
        {
            if (started)
            {
                Stop(process);
                try { await process.WaitForExitAsync().ConfigureAwait(false); }
                catch (InvalidOperationException) { }
            }
        }
    }

    private static void Stop(Process process)
    {
        try { if (!process.HasExited) process.Kill(entireProcessTree: true); }
        catch (Exception error) when (error is InvalidOperationException or Win32Exception) { }
    }
}

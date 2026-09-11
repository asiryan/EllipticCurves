using EllipticCurves.Visualizer.Computations;

if (args.Contains("--unresponsive"))
{
    Console.WriteLine(System.Text.Json.JsonSerializer.Serialize(new CalculationUpdate("progress", "Unresponsive test operation")));
    await Task.Delay(Timeout.InfiniteTimeSpan);
    return 0;
}
return await CalculationWorker.RunAsync(Console.In, Console.Out);

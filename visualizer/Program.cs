using System.Text;
using System.IO;
using EllipticCurves.Visualizer.Computations;

namespace EllipticCurves.Visualizer;

public static class Program
{
    [STAThread]
    public static int Main(string[] args)
    {
        if (args.Length == 1 && args[0] == "--compute-worker")
        {
            // A WinExe has no console. Use the redirected handles directly;
            // changing Console.OutputEncoding tries to access a nonexistent console.
            using var input = new StreamReader(Console.OpenStandardInput(), new UTF8Encoding(false));
            using var output = new StreamWriter(Console.OpenStandardOutput(), new UTF8Encoding(false)) { AutoFlush = true };
            return CalculationWorker.RunAsync(input, output).GetAwaiter().GetResult();
        }
        var app = new App();
        app.InitializeComponent();
        return app.Run();
    }
}

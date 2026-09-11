using EllipticCurves.Explorer.Models;

namespace EllipticCurves.ConsoleApp;

internal static class ConsoleApplication
{
    internal static int Run(string[] args, TextWriter output, TextWriter error)
    {
        CommandLineOptions options;
        try { options = CommandLineOptions.Parse(args); }
        catch (ArgumentException exception) { return InvalidArguments(exception.Message, error); }

        if (options.ShowHelp)
        {
            output.WriteLine(CommandLineOptions.HelpText);
            return 0;
        }

        if (!CurveEquationText.TryParse(options.Equation, out var curve, out var parseError))
            return InvalidArguments(parseError, error);
        if (curve!.Discriminant.IsZero)
            return InvalidArguments("The equation is singular (discriminant = 0). Enter a nonsingular elliptic curve.", error);

        try
        {
            CurveReport.Write(curve, options.UseLmfdb, output);
            return 0;
        }
        catch (Exception exception)
        {
            error.WriteLine($"Error: {exception.Message}");
            return 1;
        }
    }

    private static int InvalidArguments(string message, TextWriter error)
    {
        error.WriteLine($"Error: {message}");
        error.WriteLine("Run EllipticCurves.Console --help for usage.");
        return 2;
    }
}

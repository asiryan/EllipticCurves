#nullable enable
namespace EllipticCurves.ConsoleApp;

internal sealed record CommandLineOptions(string Equation, bool UseLmfdb, bool ShowHelp)
{
    public const string DefaultEquation = "Y^2 = X^3 - 17 X^2 + 72 X";

    public static string HelpText => $"""
        Elliptic Curves Console

        Usage: EllipticCurves.Console [--curve "equation"] [--lmfdb [true|false]]

        Options:
          --curve "equation"   Weierstrass equation over Q.
                               Default: {DefaultEquation}
          --lmfdb [true|false] Fetch and compare LMFDB data. Default: false.
                               --lmfdb alone means true; requires internet access.
          --help, -h          Show this help without running calculations.

        Quote the entire equation. Use ^ for powers; * is optional.
        Fractions, decimals and scientific notation are supported, as in Explorer.
        Both --curve="equation" and --lmfdb=true are also accepted.

        Examples:
          EllipticCurves.Console
          EllipticCurves.Console --curve "y^2 + y = x^3 - x"
          EllipticCurves.Console --curve "y^2 = x^3 - 1.5x + 1/2" --lmfdb

        Calculations run locally unless --lmfdb is enabled. Use Ctrl+C to stop.
        Exit codes: 0 = success/help, 1 = calculation or lookup failed, 2 = invalid input.
        """;

    public static CommandLineOptions Parse(string[] args)
    {
        if (args.Any(argument => argument is "--help" or "-h"))
            return new(DefaultEquation, false, true);

        var equation = DefaultEquation;
        var useLmfdb = false;
        var seen = new HashSet<string>(StringComparer.Ordinal);
        for (var index = 0; index < args.Length; index++)
        {
            var argument = args[index];
            var separator = argument.IndexOf('=');
            var name = separator < 0 ? argument : argument[..separator];
            var value = separator < 0 ? null : argument[(separator + 1)..];
            if (name is not ("--curve" or "--lmfdb"))
                throw new ArgumentException($"Unknown argument '{argument}'. Quote the entire equation after --curve.");
            if (!seen.Add(name)) throw new ArgumentException($"Option '{name}' was supplied more than once.");

            if (name == "--curve")
            {
                value ??= NextValue(args, ref index);
                if (string.IsNullOrWhiteSpace(value))
                    throw new ArgumentException("--curve requires a quoted equation.");
                equation = value;
            }
            else
            {
                value ??= NextValue(args, ref index);
                if (value is null) useLmfdb = true;
                else if (!bool.TryParse(value, out useLmfdb))
                    throw new ArgumentException("--lmfdb accepts only true or false, or no value to enable it.");
            }
        }
        return new(equation, useLmfdb, false);
    }

    private static string? NextValue(string[] args, ref int index)
    {
        if (index + 1 >= args.Length || args[index + 1].StartsWith("--", StringComparison.Ordinal)) return null;
        return args[++index];
    }
}

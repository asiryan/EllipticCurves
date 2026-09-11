using EllipticCurves.ConsoleApp;
using EllipticCurves.Explorer.Models;
using Xunit;

namespace EllipticCurves.Tests;

public sealed class ConsoleCommandLineTests
{
    [Fact]
    public void DefaultsToTheRequestedCurveWithoutLmfdb()
    {
        var options = CommandLineOptions.Parse([]);
        Assert.False(options.ShowHelp);
        Assert.False(options.UseLmfdb);
        Assert.True(CurveEquationText.TryParse(options.Equation, out var curve, out var error), error);
        Assert.Equal(new EllipticCurveQ(0, -17, 0, 72, 0), curve);
    }

    [Theory]
    [InlineData("--lmfdb", true)]
    [InlineData("--lmfdb=true", true)]
    [InlineData("--lmfdb=false", false)]
    [InlineData("--lmfdb=FALSE", false)]
    public void LmfdbIsAnExplicitOptIn(string argument, bool enabled)
        => Assert.Equal(enabled, CommandLineOptions.Parse([argument]).UseLmfdb);

    [Theory]
    [InlineData("true", true)]
    [InlineData("false", false)]
    public void AcceptsSeparatedBooleanValues(string value, bool enabled)
        => Assert.Equal(enabled, CommandLineOptions.Parse(["--lmfdb", value]).UseLmfdb);

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void CurveArgumentPreservesTheEqualsSignAndExactCoefficients(bool inline)
    {
        const string equation = "Y^2 + XY + Y = X^3 - 1.5 X + 1/2";
        string[] args = inline ? ["--curve=" + equation, "--lmfdb", "false"]
            : ["--lmfdb", "false", "--curve", equation];
        var options = CommandLineOptions.Parse(args);
        Assert.Equal(equation, options.Equation);
        Assert.False(options.UseLmfdb);
        Assert.True(CurveEquationText.TryParse(options.Equation, out var curve, out var error), error);
        Assert.Equal(new EllipticCurveQ(1, 0, 1, new BigRational(-3, 2), new BigRational(1, 2)), curve);
    }

    [Fact]
    public void BareLmfdbFlagDoesNotConsumeTheFollowingCurveOption()
    {
        var options = CommandLineOptions.Parse(["--lmfdb", "--curve", "y^2 = x^3 - x"]);
        Assert.True(options.UseLmfdb);
        Assert.Equal("y^2 = x^3 - x", options.Equation);
    }

    public static IEnumerable<object[]> InvalidArguments()
    {
        yield return [new[] { "--curve" }];
        yield return [new[] { "--curve=" }];
        yield return [new[] { "--curve", " " }];
        yield return [new[] { "--curve", "--lmfdb" }];
        yield return [new[] { "--curve", "y^2 = x^3 - x", "--curve", "y^2 = x^3 + 1" }];
        yield return [new[] { "--lmfdb=yes" }];
        yield return [new[] { "--lmfdb=" }];
        yield return [new[] { "--lmfdb", "0" }];
        yield return [new[] { "--lmfdb", "--lmfdb=false" }];
        yield return [new[] { "--unknown" }];
        yield return [new[] { "--curve", "y^2", "=", "x^3", "-", "x" }];
        yield return [new[] { "y^2 = x^3 - x" }];
        yield return [new[] { "--curve", "y^2 = x^3 + 1/0" }];
        yield return [new[] { "--curve", "y^2 = x^3" }];
    }

    [Theory]
    [MemberData(nameof(InvalidArguments))]
    public void InvalidInputReturnsCodeTwoWithoutStartingAReport(string[] args)
    {
        using var output = new StringWriter();
        using var error = new StringWriter();
        Assert.Equal(2, ConsoleApplication.Run(args, output, error));
        Assert.Equal("", output.ToString());
        Assert.Contains("Error:", error.ToString());
        Assert.Contains("--help", error.ToString());
    }

    [Theory]
    [InlineData("--help")]
    [InlineData("-h")]
    public void HelpDoesNotValidateOrComputeTheCurveOrFetchLmfdb(string help)
    {
        using var output = new StringWriter();
        using var error = new StringWriter();
        Assert.Equal(0, ConsoleApplication.Run(["--curve", "not a curve", "--lmfdb", help], output, error));
        Assert.Contains(CommandLineOptions.DefaultEquation, output.ToString());
        Assert.Contains("--lmfdb", output.ToString());
        Assert.DoesNotContain("Computing", output.ToString());
        Assert.Equal("", error.ToString());
    }
}

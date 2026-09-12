using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;

namespace EllipticCurves.Explorer;

public static class ExplorerInfo
{
    public const string WindowTitle = "Elliptic Curves · Explorer";
    public const string RepositoryUrl = "https://github.com/asiryan/EllipticCurves";
    public const string UserGuideUrl = RepositoryUrl + "/blob/main/explorer/README.md";
    public const string ReportIssueUrl = RepositoryUrl + "/issues/new";
    public const string LmfdbUrl = "https://www.lmfdb.org/";

    private static readonly Assembly Library = typeof(EllipticCurveQ).Assembly;
    public static string LibraryVersion => Library.GetName().Version!.ToString(3);
    public static string Author => Library.GetCustomAttribute<AssemblyCompanyAttribute>()!.Company;
    public static string Copyright => Library.GetCustomAttribute<AssemblyCopyrightAttribute>()!.Copyright;
    public static string VersionInfo => $"{WindowTitle}\nEllipticCurves library: {LibraryVersion}\n"
        + $"Explorer build: {typeof(ExplorerInfo).Assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>()!.InformationalVersion}\n"
        + $"{RuntimeInformation.OSDescription}\n{RuntimeInformation.FrameworkDescription} · {RuntimeInformation.ProcessArchitecture}";
    public static string LicenseText
    {
        get
        {
            using var stream = typeof(ExplorerInfo).Assembly.GetManifestResourceStream("EllipticCurves.Explorer.LICENSE")
                ?? throw new InvalidOperationException("The bundled MIT license is missing.");
            using var reader = new StreamReader(stream);
            return reader.ReadToEnd();
        }
    }
}

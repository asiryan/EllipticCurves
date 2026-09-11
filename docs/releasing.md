# Release preparation

Run the commands below from the repository root. They build, test or package the
current checkout; they do not upload a NuGet package or publish a GitHub release.
Use the .NET 8 SDK. A newer SDK can also target these projects, but running the
tests and framework-dependent applications still requires the .NET 8 runtime;
WPF requires the Windows Desktop runtime. Explorer and its presentation check
run on Windows.

## Build script

On Windows, run [build.bat](../build.bat) to build the solution in Release, run
the arithmetic, Console and portable Explorer tests and the WPF presentation
check, pack the library, and publish Explorer and the console application:

```powershell
.\build.bat
# Optional target for both applications:
.\build.bat win-arm64
```

The default target for both applications is `win-x64`. Each invocation creates a fresh
`artifacts/release/TIMESTAMP-RUNTIME` directory containing `nuget/`,
`explorer-RUNTIME/`, `console-RUNTIME/`, `EllipticCurves.Explorer-RUNTIME.zip` and
`EllipticCurves.Console-RUNTIME.zip`. Both application archives include the license
and their .NET runtimes; neither requires a separate .NET installation to run.
Versions come from the project files described below.

The script works from any current directory, stops on the first failed command
and returns a nonzero exit code. A failed run can leave partial output in its
directory. It does not delete previous releases or upload artifacts. Windows
PowerShell is used for timestamps and ZIP creation. Tests run on the build
machine; an ARM64 package still needs validation on an ARM64 Windows machine.
Run `build.bat --help` for usage. The individual commands follow below.

## Components and versions

| Component | Project | Target | Distribution |
| --- | --- | --- | --- |
| Library | [sources/EllipticCurves.csproj](../sources/EllipticCurves.csproj) | .NET Standard 2.0 | NuGet package |
| Console | [console/EllipticCurves.Console.csproj](../console/EllipticCurves.Console.csproj) | .NET 8 | Complete self-contained publish folder in a ZIP |
| Desktop Explorer | [explorer/EllipticCurves.Explorer.csproj](../explorer/EllipticCurves.Explorer.csproj) | .NET 8, Windows | Complete publish folder in a ZIP |

The library sets `Version`, `AssemblyVersion` and `FileVersion` in its project
file. Use that file as the source of truth for the NuGet version.
The Console and Explorer projects currently have no explicit version settings
and use the SDK's default `Version`. A project reference does not
inherit the library's version. If all release artifacts should share a version,
set the applications' version properties explicitly before packaging them.

The NuGet README is [docs/nuget-readme.md](nuget-readme.md), not the root README.
The library project packs it as `README.md`, together with `LICENSE.md`,
`ec_logo.png`, the library and its generated XML API documentation.
Console and Explorer are excluded from NuGet packaging. Pack the library project explicitly
rather than packing the entire solution.

## Validation

Run all arithmetic, Console command-line and portable Explorer tests:

```powershell
dotnet test tests/EllipticCurves.Tests.csproj -c Release
```

On Windows, also run the WPF presentation check:

```powershell
dotnet run --project tests/PresentationHost/PresentationHost.csproj -c Release
```

The second command checks compiled XAML, layout and selected interaction paths;
it is not part of the solution's `dotnet test` run and does not show application
windows. Tests use committed fixtures and simulated HTTP responses. Restoring
SDK/NuGet dependencies can require internet access; the arithmetic tests do not
call LMFDB. Console only performs a live LMFDB lookup when `--lmfdb` is enabled.

## NuGet package

```powershell
dotnet pack sources/EllipticCurves.csproj -c Release -p:GeneratePackageOnBuild=false -o artifacts/nuget
```

This command builds and explicitly packs the library. Disabling automatic packing
for this invocation avoids coupling the pack operation to the project's normal
`GeneratePackageOnBuild` behavior. The package is written to `artifacts/nuget/`,
with its filename derived from the package ID and version in the project file.

Before uploading, inspect the package archive for the `lib/netstandard2.0` DLL
and XML documentation, the README, license and icon. Confirm that the `.nuspec`
version and dependency declarations match the project. Package inspection and a
successful install in a consumer project are separate from the source tests.

## Explorer archive

Use a fresh output directory for each release and architecture so that files
from an older publish are not included in the archive.

```powershell
dotnet publish explorer/EllipticCurves.Explorer.csproj -c Release -r win-x64 --self-contained true -o artifacts/explorer-win-x64
Copy-Item -LiteralPath LICENSE -Destination artifacts/explorer-win-x64/EllipticCurves.LICENSE.txt
Compress-Archive -Path artifacts/explorer-win-x64/* -DestinationPath artifacts/EllipticCurves.Explorer-win-x64.zip
```

The ZIP must contain the entire publish folder, including runtime files and
`EllipticCurves.dll`. This is a self-contained folder deployment, not a standalone
single-file executable. Unpack it and run `EllipticCurves.Explorer.exe` on Windows.
For ARM64, replace `win-x64` with `win-arm64` in the runtime, output directory and
archive name, and validate that build on an appropriate Windows machine.

Check the published application itself before uploading: open Explorer and run
a torsion calculation, verify that its result appears, stop a running calculation,
save a report and export PNGs from both **Real locus** and **Complex torus**.
For the complex view, check period preparation, shared point selection, rotation,
zoom and **Reset view**; a singular cubic should show an explanation instead of a
torus. Inspect both exported images for dark backgrounds, correct bounds and
unclipped visible content. Repeat at a compact window size, where the complex
view can scroll and export captures only its visible viewport.

Also check panel resizing/folding and the Clear and Reset confirmations.
Exercise live LMFDB fetching separately when internet access is available.
Testing the packaged executable verifies the calculation worker's startup and
published dependencies as well as the UI.

## Console archive

```powershell
dotnet publish console/EllipticCurves.Console.csproj -c Release -r win-x64 --self-contained true -o artifacts/console-win-x64
Copy-Item -LiteralPath LICENSE -Destination artifacts/console-win-x64/EllipticCurves.LICENSE.txt
Compress-Archive -Path artifacts/console-win-x64/* -DestinationPath artifacts/EllipticCurves.Console-win-x64.zip
```

Use a fresh output directory and keep the entire publish folder, including the
runtime and library files. From a terminal in the extracted folder, run
`EllipticCurves.Console.exe`; a separate .NET installation is not required.
The executable uses Explorer's icon. For ARM64, replace `win-x64` with
`win-arm64` in the runtime, output directory and archive name.

The default equation is `Y^2 = X^3 - 17 X^2 + 72 X`, and LMFDB lookup defaults to
`false`. Check the packaged application with no arguments, a custom
`--curve "y^2 + y = x^3 - x"`, `--lmfdb false` and `--help`. Test an invalid
equation and confirm exit code `2`, with a message on standard error and no
calculation started. Test `--lmfdb` separately with internet access: the local
report precedes the lookup, and lookup failure returns exit code `1`.
The [Console README](../console/README.md) documents all argument forms and output.

## Publication order

Publish the corresponding source and documentation before uploading the package.
The packaged README currently links to the repository's `main` branch, including
`explorer/README.md`. Those paths must exist publicly when the package is released;
files present only in a local checkout or development branch do not make the links
work. Alternatively, use links to the published release tag in the packaged README.

Publish the package and archives only after their checks pass. An unversioned
`dotnet add package EllipticCurves` installs the latest published stable package,
which can be older than the checkout being documented. Release notes should name
the component versions, distinguish numerical estimates and stored database values
from certified results, and retain the documented limits on rank determination,
point searches and saturation.

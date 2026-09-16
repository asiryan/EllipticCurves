# Elliptic Curves Console

A command-line application built on the [EllipticCurves library](../README.md).
It computes invariants, rational torsion points, a minimal model, the conductor
and proved rank bounds. It also estimates the analytic rank and attempts a
rigorous rank 0/1 certificate. Calculations run locally by default; LMFDB lookup
and comparison are optional.

The default curve is **y² = x³ − 17x² + 72x**, studied in
[this paper](https://arxiv.org/abs/2510.11768).

## Run

Unpack the complete `EllipticCurves.Console.VERSION-win-x64.zip` release archive,
where `VERSION` is the library package version. Choose `win-arm64` for Windows ARM64.
No separate .NET installation is needed.
Open a terminal in the extracted folder and run:

```powershell
.\EllipticCurves.Console.exe
.\EllipticCurves.Console.exe --curve "y^2 + y = x^3 - x"
.\EllipticCurves.Console.exe --curve "y^2 = x^3 - 1.5x + 1/2" --lmfdb
.\EllipticCurves.Console.exe --help
```

Keep the runtime files alongside the executable. The Windows executable uses
the same icon as Explorer. The application does not pause at the end, so run it
from a terminal to keep the output visible. Use `Ctrl+C` to stop a calculation.

## Arguments

| Argument | Default | Meaning |
| --- | --- | --- |
| `--curve "equation"` | `Y^2 = X^3 - 17 X^2 + 72 X` | The nonsingular Weierstrass equation to study. |
| `--lmfdb [true\|false]` | `false` | Fetch LMFDB metadata and compare it with the local results. The bare flag enables lookup. |
| `--help`, `-h` | | Print help without running calculations or accessing the network. |

Both `--curve="equation"` and `--lmfdb=true` are accepted. Use `--lmfdb false`
or `--lmfdb=false` to explicitly disable lookup. Options may appear in either
order; unknown options, repeated options and invalid values are rejected.

Quote the entire equation. Parsing is shared with Explorer: `x` and `y` are
case-insensitive, powers use `^`, multiplication may use `*` or be implicit,
and coefficients may be exact fractions, decimals or scientific notation.
General equations such as `y^2 + xy + y = x^3 - x` are supported. Singular
equations (zero discriminant) are rejected before calculations start.

The program writes results to standard output and errors to standard error.
Exit codes are `0` for success/help, `1` for a calculation or LMFDB failure
(including a contradiction with the native arithmetic results), and `2` for
invalid arguments or equations. An inconclusive rank estimate remains a result,
with its evidence status and explanation; it is not reported as a proved rank.

Only `--lmfdb` requires internet access. The local report is printed first, so it
remains available if the lookup fails or the curve is absent from the database.
Analytic-rank comparisons report `Unknown` if either estimate is missing;
a disagreement with a numerical estimate alone is not treated as a proof failure.

## Run from source

Install the .NET 8 SDK, then run from the repository root:

```sh
dotnet run --project console/EllipticCurves.Console.csproj
dotnet run --project console/EllipticCurves.Console.csproj -- --curve "y^2 + y = x^3 - x" --lmfdb
```

Arguments after `--` are passed to the application. Alternatively, run
`dotnet run -- --help` from `console`. The project references the library source
directly and shares Explorer's equation parser without depending on WPF.
Restoring SDK/NuGet dependencies may require internet access.

## Build a release

Run [build.bat](../build.bat) from the repository root. It builds and tests the
solution and publishes both Console and Explorer with their .NET runtimes for
`win-x64` (or `win-arm64` when supplied as the argument).

To publish only Console:

```powershell
dotnet publish console/EllipticCurves.Console.csproj -c Release -r win-x64 --self-contained true -o artifacts/console-win-x64
```

See [release preparation](../docs/releasing.md#console-archive) for packaging.
The SDK is needed on the build machine, not on the user's machine.

## Example output

The default local report begins as follows (torsion point order may vary).
It continues with the native minimal model, conductor, rank bounds and
analytic-rank estimate, including their proof status and explanations.

```text
E: y^2 = x^3 - 17*x^2 + 72*x
Short Weierstrass: y^2 = x^3 - 73/3*x + 1190/27
b2 = -68
b4 = 144
b6 = 0
b8 = -5184
D  = 82944
c4 = 1168
c6 = -38080
j  = 1556068/81

Computing rational torsion...
Torsion: Z/2Z x Z/4Z
Torsion points:
O
(0, 0)
(8, 0)
(9, 0)
(6, 6)
(6, -6)
(12, 12)
(12, -12)
```

With `--lmfdb`, an additional section follows the local report. For the default
curve it includes:

```text
Fetching LMFDB data...
LMFDB: 48.a3
Url: https://www.lmfdb.org/EllipticCurve/Q/48.a3/
Minimal Weierstrass model: y^2 = x^3 + x^2 - 24*x + 36
Torsion: Z/2Z x Z/4Z
Rank(E) = 0
Analytic rank(E) = 0
Cond(E) = 48
Isomorphic to E: True
Native minimal model matches LMFDB: True
Native conductor matches LMFDB: True
LMFDB rank is within native bounds: True
```

See the [library README](../README.md#rank-bounds-and-conductor) for the native
rank APIs and the meaning of their proof and certification statuses.

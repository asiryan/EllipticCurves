# Elliptic Curves Explorer

A WPF desktop application on .NET 8 for exploring the real geometry and basic
arithmetic of elliptic curves. All calculations use the local EllipticCurves
library. The application makes no LMFDB requests or other network calls.

![Elliptic Curves Explorer](../docs/png/visualizer.png)

## Run

On Windows, install the .NET 8 SDK, then run this command from the repository root:

```powershell
dotnet run --project visualizer/EllipticCurves.Visualizer.csproj -c Release
```

Alternatively, open `EllipticCurves.sln` in an IDE that supports .NET 8 and choose
`EllipticCurves.Visualizer` as the startup project. The application is separate from
the library's NuGet package.

To produce a folder that runs without an installed .NET runtime:

```powershell
dotnet publish visualizer/EllipticCurves.Visualizer.csproj -c Release -r win-x64 --self-contained true -o artifacts/visualizer-win-x64
```

Run `EllipticCurves.Visualizer.exe` from that folder. Use `win-arm64` instead of
`win-x64` for an ARM64 build.

## Explore

- Type the entire equation in one field, for example
  `y^2 = x^3 - 106.16*x - 0.32` or `y^2 + xy + y = x^3 - x`.
  Simple and general Weierstrass forms are recognized automatically.
- Enter powers only with `^`, such as `x^3` and `y^2`. Multiplication can be
  explicit (`2*x*y`) or implicit (`2xy`). Parentheses, rational constant division,
  and rearranged polynomial equations are supported. The result must reduce to
  `y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6`; variables in denominators,
  other curve families, and named functions are rejected with an explanation.
- Enter exact decimals (`8.325` or `8,325`), fractions (`-2/7`) or scientific
  notation (`1e-5`). There is no ±100 coefficient limit or one-decimal-place
  restriction. Input is parsed directly into rational numbers, without rounding.
- Open **Experiment** for optional sliders. Each slider moves up to 50 exact
  steps on either side of its anchor. Set a positive **Slider step** there, or
  choose `1`, `0.1` or `0.01`. Changing a coefficient in the formula or changing
  the step recenters its range. Sliders update the formula using exact arithmetic.
  **Right-click a slider to reset its coefficient to 0**, including when the
  formula field contains invalid input.
- See the discriminant, j-invariant, c₄, c₆, real component count and short model
  update after a 300 ms pause in valid input; **Enter** applies it immediately.
  The graph always uses the entered coordinates; the short
  model is shown separately in transformed coordinates.
- Drag the plot to pan and use the mouse wheel to zoom about the pointer. Both
  axes use the same scale. Edits preserve the viewport. **Fit curve**, a double-click
  or **Ctrl+F** recenters the view around the real branch points and part of the
  unbounded branch. With the plot focused, **Home** fits and **+ / −** zoom.
  Slider arrow keys move by the chosen exact step.
- Hover near a branch to inspect approximate coordinates, or near a gold marker
  to see an exact rational point.
- Choose one of four built-in examples, including a singular cubic. Preset names
  such as `37.a1` are static labels and do not trigger database lookups.
- Copy the equation and exact invariants, or export the current plot to PNG.

## Calculation scope

The real locus uses double-precision drawing coordinates and is a numerical
visualization. Completing the square handles both branches and the linear y
terms; sampling is split at real roots to preserve disconnected components.
The point at infinity is not drawn in the affine plot.
Inputs outside the numerical plot's representable range retain exact invariants
and show a plot precision message. To bound input processing, text is limited to
4096 characters and scientific exponents to ±4096.

Gold markers are exact affine rational points found with `RationalPoints(12, 4)`:
their x-coordinates are `m/n`, with `|m| ≤ 12` and `1 ≤ n ≤ 4`. This is a bounded
sample, not a complete list of rational points or a Mordell–Weil basis. The search
runs in the background after a 120 ms debounce. Superseded searches are cancelled,
and markers are cleared when a new curve is applied. While editing, existing
markers continue to belong to the displayed curve.

For Δ = 0 the cubic is still drawn, but j and the elliptic real-component count
are marked undefined, and rational sampling is disabled. Invalid or incomplete
text input retains the last valid curve with a neutral editing status. Validation
is shown after leaving the field or pressing **Enter**, so partial input such as
`-` or `1/` is not highlighted while typing.

Rank, conductor, torsion enumeration, heights and periods are not automatically
computed by this application. Coefficient updates invoke only the inexpensive native
invariants and the bounded sample search described above.

## Development

The XAML theme, plot control, immutable calculation snapshots and view models are
separate files. No external charting or UI package is required. The portable model
and view-model sources are linked into the existing test project, so their tests
also run without WPF on non-Windows systems:

```powershell
dotnet test tests/EllipticCurves.Tests.csproj -c Release --filter "FullyQualifiedName~VisualizerModelTests|FullyQualifiedName~VisualizerEquationTests"
```

The desktop project uses the [Microsoft .NET Desktop SDK settings](https://learn.microsoft.com/en-us/dotnet/core/project-sdk/msbuild-props-desktop).

The application logo is `ec_logo.png`; it is embedded as a WPF resource for the
header and window icon. After replacing the PNG, run `./visualizer/tools/Update-Icon.ps1`
from the repository root to regenerate the executable's multi-size `ec_logo.ico`,
then rebuild the application.

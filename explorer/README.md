# Elliptic Curves Explorer

A WPF desktop application on .NET 8 for exploring elliptic-curve geometry and
the full computational API of the EllipticCurves library. Native calculations run
locally. The explicit LMFDB search and fetch commands are the only operations that use the
network; plotting and editing do not make network requests.

![Elliptic Curves Explorer](../docs/png/ec_explorer.png)

## Run

On Windows, install the .NET 8 SDK, then run this command from the repository root:

```powershell
dotnet run --project explorer/EllipticCurves.Explorer.csproj -c Release
```

Alternatively, open `EllipticCurves.sln` in an IDE that supports .NET 8 and choose
`EllipticCurves.Explorer` as the startup project. The application is separate from
the library's NuGet package.

To produce a folder that runs without an installed .NET runtime:

```powershell
dotnet publish explorer/EllipticCurves.Explorer.csproj -c Release -r win-x64 --self-contained true -o artifacts/explorer-win-x64
```

Run `EllipticCurves.Explorer.exe` from that folder. Distribute the entire publish
folder, including its libraries and runtime files; the executable alone is not
self-contained. For an ARM64 build, use `-r win-arm64` and a separate output
folder such as `artifacts/explorer-win-arm64`.
See [release preparation](../docs/releasing.md) for version settings, packaging
and the checks to run before publishing a release.

## Sessions

**File**, before **Edit** and **Tools** in the title bar, contains **New**,
**Open**, **Save**, **Save as** and **Exit**. New starts with the classic curve, empty history
and default visualization and panel settings. Open and Save use `.ec` session files;
Exit closes Explorer.
Use **Ctrl+N**, **Ctrl+O**, **Ctrl+S** and **Ctrl+Shift+S** for New, Open, Save and
Save as; the shortcuts are also displayed beside their menu items. Save and Ctrl+S
are available only for an existing file with unsaved changes or a failed save to retry.
They are disabled while the status is **Saved**. Use Save as for the first save of
a new session. Save writes to the current file without a picker. Save as always
opens the file picker and makes the chosen file the current session after a successful save.
For an existing session, the picker opens in its folder with its current file name.
Selecting an existing file asks for overwrite confirmation. Confirming replaces
that exact file; cancelling leaves it untouched. Explorer never adds a numeric suffix
automatically. To keep a separate copy, choose a different name. The title bar and
subsequent Save commands use the exact path selected in the picker. This also applies
to the first save and to a renamed session from the unsaved-changes dialog.
The title bar shows its file name (the full path appears on hover), an asterisk for
unsaved changes, and **New**, **Unsaved**, **Saving…**, **Saved** or
**Save failed**. File writing runs in the background; edits made during a save remain
unsaved. Other session commands are disabled while writing.
The menu uses the same header and popup styling as Tools, with a single vertical
list. All four menus close on Escape, another menu, an outside click (including the
title bar), window movement, resizing or deactivation.

A session file preserves the exact equation and up to 50 calculation reports with their original
inputs, limits and timestamps. Opening a session selects the newest result;
browsing other results does not trigger an unsaved-changes warning. Adding,
deleting or changing results still does. Samples and periods are recomputed
locally as needed. Open calculation parameter windows are not saved and close when
another session is opened.
Slider step and positions, the preset name, visualization mode, cameras, grid and sample visibility, and selected torus point
belong to the current window and its undo history. Panel layout and scroll positions are local to the window
and are not restored by Undo or Redo. None of these settings are written to the file. Opening a session starts in **Real locus** with default display settings and fits the
curve to the current plot size. Panning, zooming,
rotating the torus and resetting the view do not trigger an unsaved-changes warning.
Saving does not move or reset the graph currently on screen.
Coefficients and a matching preset are recovered from the equation. The slider step
starts at **0.01**, with every slider centered at **0** around its coefficient.

New and Open are disabled during a calculation; stop it first. Save can capture a running
calculation, which reopens as **Interrupted**, without automatically starting work
or network requests. Use **Repeat** to reopen its parameters. Finish an incomplete
equation before saving; a local slider-step error does not block saving curve data.

Files contain version 1 JSON with only `Format`, `Version`, `Equation` and `History`.
Data is validated before replacing the workspace; no format migrations are performed.
An invalid or unsupported file leaves the current session intact. Saving writes
a temporary file before replacing the destination. Session files are limited to
256 MB. Before New, Open or closing the application (including the window close
button and Alt+F4), unsaved changes prompt **Save**, **Discard** or **Cancel** in
the application's dark dialog. Open shows this confirmation before the file picker.
The file name in this confirmation is editable. An unchanged name saves directly to
the current file. New or renamed sessions open the file picker with `.ec` added when
needed and the current session's folder preserved. Cancelling does not change the
session's name or location. If new edits arrive during saving, the pending New,
Open or Close is cancelled so those edits remain available.
Enter, Escape and the dialog close button cancel.
Cancelling the file picker or failing to save also cancels the pending action.
An untouched session or an unchanged saved/opened session does not prompt.
Changes to the equation and results are tracked. Changing the slider step or switching
between 2D and 3D, toggling the grid or samples, selecting a torus point, scrolling
or rearranging panels does not count as a data edit or enable Save. Save and Save as
write the document data and leave the current visualization unchanged.
Background sample generation does not count as an edit. Sessions are saved
explicitly; there is no automatic saving on exit.
The Results panel's **Export** exports an individual text report.

## Help

**Help → User Guide** (**F1**) opens this guide in the browser. **Keyboard Shortcuts**
opens a local reference for the main window's commands. The menu also links to
**LMFDB Website**, **Project on GitHub** and **Report an Issue**; the last item opens
GitHub's new-issue page for you to fill in and submit.

**About EllipticCurves** shows the version and author metadata from the loaded
EllipticCurves library. **Copy version info** copies the library version, Explorer
build, operating system, .NET runtime and process architecture for issue reports.
**MIT License** opens the full project license bundled with the application.
About, shortcuts and license windows work offline and are modeless, so the main
window remains usable. Reopening a page activates its existing window. **Esc** closes it.

## Undo and redo

**Edit → Undo** (**Ctrl+Z**) and **Edit → Redo** (**Ctrl+Y**) restore up to 50 previous
workspace edits: equations (including incomplete input), coefficients, exact slider
steps and positions, calculation reports and the selected result, graph mode, grid and point visibility,
the real plot's position and zoom, and the torus camera and selected point.
Typing and repeated wheel/slider changes are grouped after a short pause; a mouse
drag remains one edit until it ends. The shortcuts also work while editing an equation.
Each selection in **Results** is a separate undo step, so stepping back through reports
does not skip straight to an earlier graph edit. Automatic selection when adding,
deleting or clearing reports belongs to that action and does not add another step.

History uses in-memory mementos with shared computed results. Undo and Redo restore
existing curve snapshots, samples, period data and reports without rerunning calculation
operations or network requests. Work that had not finished preparing a graph when you
left it is restored as incomplete. Reopen **Complex torus** to resume its unfinished preparation;
Undo and Redo themselves do not start work. Deleting a report, clearing all reports and adding
a completed, failed, timed-out or stopped calculation are undoable. Undo also restores
the oldest report if a new calculation displaced it from the 50-report session limit.
Undo and Redo are disabled while a calculation or session file operation is active.

A new edit after Undo clears Redo. New and Open start a fresh undo history;
Save preserves it. Undoing data changes back to the saved state clears the unsaved
indicator. Graph navigation stays separate from the session's unsaved-data indicator.
Undo history is kept only for the current run and is not written to `.ec` files.

## Explore

- Type the entire equation in one field, for example
  `y^2 = x^3 - 106.16*x - 0.32` or `y^2 + xy + y = x^3 - x`.
  Simple and general Weierstrass forms are recognized automatically.
- Enter powers only with `^` and integer exponents from 0 to 3, such as `x^3`
  and `y^2`. Intermediate polynomials must also have degree at most 3.
  Multiplication can be explicit (`2*x*y`) or implicit (`2xy`). Parentheses, rational constant division,
  and rearranged polynomial equations are supported. The result must reduce to
  `y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6`; variables in denominators,
  other curve families, and named functions are rejected with an explanation.
- Enter exact decimals (`8.325` or `8,325`), fractions (`-2/7`) or scientific
  notation (`1e-5`). There is no ±100 coefficient limit or one-decimal-place
  restriction. Input is parsed directly into rational numbers, without rounding.
  **EQUATION OVER ℚ** shows the applied equation with exact rational coefficients,
  for example `3/2` for an input coefficient of `1.5`.
- Open **Coefficients** for optional sliders. Each slider moves up to 50 exact
  steps on either side of its anchor. Set a positive **Slider step** there, or
  choose `1`, `0.1` or `0.01`. Changing a coefficient in the formula or changing
  the step recenters its range. Leaving the step unchanged, including an equivalent
  value such as `1/100` for `0.01`, preserves the slider positions and ranges.
  Sliders update the formula using exact arithmetic.
  **Right-click a slider to reset its coefficient to 0**, including when the
  formula field contains invalid input.
- See the discriminant, j-invariant, c₄, c₆, real component count and short model
  update after a 300 ms pause in valid input; **Enter** applies it immediately.
  The real-locus graph uses the entered coordinates; the short
  model is shown separately in transformed coordinates. Confirming unchanged input
  keeps the selected example and current samples. Equivalent equations reuse the
  samples; opening a calculation for the same curve preserves the torus selection.
- Drag the plot to pan and use the mouse wheel to zoom about the pointer. Ordinary
  curves use equal axis units. When the coordinate ranges differ greatly, fitting
  uses separate x/y scales so that large curves remain visible.
  The plot says **Independent axis scales** in that case; labels
  and point coordinates still refer to the entered equation. Small edits preserve
  the viewport; entering a curve with a substantially different coordinate range
  automatically fits it after the normal input delay or pressing Enter. This also
  works when returning from Complex torus to Real locus. **Reset view**, a double-click
  or **Ctrl+F** recenters the view around the real branch points and part of the
  unbounded branch. With the plot focused, **Home** fits and **+ / −** zoom.
  Slider arrow keys move by the chosen exact step.
- Hover near a branch to inspect approximate coordinates, or near a gold marker
  to see an exact rational point.
- Choose one of four built-in examples, including a singular cubic. Preset names
  such as `37.a1` are static labels and do not trigger database lookups.
  Choosing an example fits the real plot and resets the torus camera. If the real
  plot is hidden, it is fitted when you next switch to **Real locus**.
- **Reset** in the Equation header asks for confirmation before restoring
  `y^2 = x^3 - x` and recentering the plot. The custom dark dialog opens centered
  on the main window without highlighting either button. Enter or Escape cancels
  from the initial state; the close button also cancels.
- Copy the equation and exact invariants, or export the current plot to PNG.

## Complex torus

The plot's mode selector switches between **Real locus · E(ℝ)** and
**Complex torus · E(ℂ)**. The complex view has two linked panels:

- **Period lattice** shows a fundamental parallelogram with its opposite edges
  identified. The drawing is scaled by ω₁, so its two basis vectors are 1 and
  τ = ω₂/ω₁. This basis is adapted to the real points; τ is not reduced to the
  modular fundamental domain. The numerical periods above it refer to the global minimal model;
  their tooltip gives that model's equation.
- **Complex torus** shows a fixed topological embedding of ℂ/Λ. The ring's
  radii do not represent the curve's complex structure; that information is in
  the lattice and τ. The turquoise and blue cycles correspond to ω₁ and ω₂.

The legend shows **Selected point** with a white dot, followed by
**Exact rational samples** with a gold dot, using the same styling as the real plot.
Gold markers represent the same bounded rational samples as the real plot;
the selected point is white in both panels. The point at infinity **O** is included
at the lattice origin. Select a point from the dropdown or click a marker in
either panel to highlight it in both.
Repeated boundary markers in the parallelogram represent the same point after
edge identification. Exact x/y coordinates remain in the entered curve's model;
their period coordinates u and v are numerical elliptic logarithms, with
z ≡ uω₁ + vω₂. Up to 32 distinct affine samples are mapped, with the mapped count
and any numerical failures shown below the view. This is not a complete point list.

Drag the torus to rotate it and scroll to zoom. With it focused, arrow keys rotate,
**Home** resets the camera and **+ / −** zoom. **Reset view**, a double-click on
the torus, or **Ctrl+F** restores the camera. Arrow keys in the lattice cycle
through the markers. **Grid** toggles the subdivisions in both panels;
**Rational samples** hides affine markers while retaining O. The real plot keeps
its viewport when switching modes unless an example or **Reset** was selected
while it was hidden. Small coefficient edits preserve that viewport in either mode;
substantially different coordinate ranges are fitted when the real plot is visible.
**Reset view** and **Ctrl+F** affect only the active mode. The **−**, **+** and **Reset view** buttons
remain in the same position in both modes. In a small window the complex view
scrolls vertically to keep its diagrams readable.

Period and point mapping calculations run locally in the background only while
the complex mode is open. They use a 180 ms debounce, a 15-second cancellation
deadline and bounded root-isolation/iteration work. Changing the curve or leaving
the mode cancels outstanding work; late results cannot replace the current curve.
The period lattice remains usable if subsequent point mapping reaches its limit.
Singular cubics (Δ = 0) show an explanation
instead of a smooth torus. No internet connection or extra graphics package is needed.

## PNG export

**Export plot** saves the active visualization with a dark background at twice
its WPF layout dimensions (192 DPI). In real-locus mode, the image contains the
graph, axes and visible rational markers. In complex mode, it contains the visible
period-lattice and torus panels, point selector and notes, preserving the current
camera and selection. Complex export becomes available once the period lattice
is ready, even if point mapping is still running.

The image excludes the surrounding equation, results and invariant panels, as
well as the shared plot toolbar, legend and navigation buttons. For a compact
complex view, the current scroll position determines what is captured; content
outside the scroll viewport is not included. Enlarge the plot area to fit more
content before exporting. Export does not reset the camera or change the live layout.

## Calculation scope

The real locus uses double-precision drawing coordinates and is a numerical
visualization. Completing the square handles both branches and the linear y
terms; sampling is split at real roots to preserve disconnected components.
The point at infinity is not drawn in the affine plot.
Inputs outside the numerical plot's representable range retain exact invariants
and show a plot precision message. To bound input processing, text is limited to
100,000 characters per equation and 20,000 per number, with scientific exponents
limited to ±4096. Each side of an equation accepts up to 4096 terms. Numerators and denominators
of intermediate expressions and normalized coefficients are limited to 32768 bits
(roughly 9800 decimal digits); parser nesting is also limited to 64 levels.
These checks bound expression growth during input processing. Long terminating
decimals may be formatted as fractions so generated equations remain readable by
the editor, calculations and session files. Slider results that exceed the coefficient
limit are reported as invalid input and cannot be saved until corrected.

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

In **Real locus** mode, coefficient updates invoke only the inexpensive native
invariants and the bounded sample search described above. The optional complex
view additionally computes periods and numerical point mappings. Rank, conductor,
torsion enumeration and heights require an explicit calculation from Tools.

## Tools

### Verify rank from supplied points

Choose **Tools → Ranks and arithmetic → Verify rank from supplied points**.
Enter one `x; y` pair per line in the current equation's coordinates. Fractions
are accepted; `O` denotes the point at infinity. Set **reduction Prime Bound**
to control the tested primes, then run the calculation.

The report contains a proved lower bound, the number of supplied points,
**Independence Certified**, the character image dimension, reduction-prime
diagnostics and an explanation of the result. Independence is certified only
when the lower bound equals the number of supplied points. A smaller bound does
not prove dependence, and zero does not prove rank zero. The method neither
searches for points nor computes an upper bound or saturation.

The calculation uses the standard results panel, session history, Repeat,
cancellation and time limits. See [the certificate construction and 31-point
example](../docs/rank31-verification.md).

### Import a curve from LMFDB

Choose **Tools → LMFDB · internet → Import curve from LMFDB**. Enter a conductor
such as `37`, or an inclusive range such as `11-100`, then choose **Search**.
Inputs run from 1 to 500000; LMFDB's complete catalog covers conductors strictly
below 500000. A conductor can have several curves, so select a labelled equation
from the list and choose **Import formula**.

Search requests only the curve labels and five exact integer coefficients. It
loads up to 100 formulas per page; **Previous** and **Next** browse the range
without downloading the whole catalog. Search requires internet access and has
a 30-second timeout; **Stop**, closing the picker or changing the input cancels
the pending request. Merely opening the picker does not make a network request.
LMFDB may occasionally return a CAPTCHA page instead of data; the picker reports
this and leaves the current curve unchanged so the search can be retried later.

Import replaces the current equation and fits the graph. Only the equation is
imported: no rank, points or other database metadata are fetched or added to the
calculation history. **Undo** restores the previous equation and view together;
**Redo** reuses the imported formula without a network request. Saving the session
preserves the equation in its normal `.ec` format. Closing the picker without
importing leaves the workspace unchanged.

### Run calculations

Open **Tools** in the title bar, choose a category or search for an operation.
Each operation opens a movable, modeless parameter window. The curve is captured
when the window opens, so editing the plot later does not silently change a pending
calculation. Finite-extension curves and rational-number tools have independent
inputs. **Run calculation** opens the results panel on the right.

The title bar contains **File**, **Edit**, **Tools** and **Help**. Clicking the logo or **ELLIPTIC CURVES** opens
the project's GitHub repository in the default browser. **Export plot** is in the
plot panel's own toolbar.

| Category | Available calculations |
| --- | --- |
| Curve and models | Exact coefficients and invariants, real components, CM discriminant, short and global minimal models, twists, construction from j |
| Rational points and torsion | Membership, addition, subtraction, negation, doubling, scalar multiplication, bounded rational/integral searches, all torsion points, torsion order and group structure |
| Ranks and arithmetic | Both rank-bound interfaces, analytic rank and certification, conductor, root number, local reduction data, Tamagawa product |
| Heights and periods | Naive, canonical, local and archimedean heights, height pairing and matrix, regulator, Faltings and stable Faltings heights, certified periods and numerical elliptic logarithms |
| Isomorphisms and isogenies | Isomorphism tests and maps, coordinate changes and inverse maps, minimal-model maps, Vélu and 2-isogenies with duals, point mapping, division and prime-by-prime saturation |
| Fourier coefficients and reduction | Individual or ranged Fourier coefficients, Frobenius traces, minimal-model point counts, reduction of curves and points |
| Prime and extension fields | Curve construction and invariants, points, group operations, enumeration and counting; point order over Fp |
| Finite-field arithmetic | Field construction, elements, addition, subtraction, negation, multiplication, division, inverses and powers |
| Rational arithmetic | Exact rational representation, arithmetic, comparison, powers, square testing and exact square roots |
| LMFDB | Fetch all supported database metadata, map database generators onto the captured model, import stored JSON without internet |

Point inputs have separate x/y fields and an infinity checkbox. Point lists use
one `x; y` pair per line (`O` for infinity). Field elements and defining
polynomials use coefficients in ascending powers of t separated by semicolons:
`0; 1` is t and `2; 0; 1` is t^2 + 2. Over Fp, curve tools reduce the **entered**
coefficients; the separate reduction operations use a **global minimal model**.

Both rank-bound actions enable parallel general descent by default, using up to
four workers (fewer when fewer processors are available). In **Precision and work
limits**, set **Max Degree Of Parallelism** to `1` for sequential execution or to
another positive worker limit. All workers share the same work allowances;
increasing the worker count does not increase **Max Descent Work**. The 2-isogeny
method remains sequential, and incomplete descent still returns an unknown upper
bound. The library itself defaults to sequential execution.

Open **Precision and work limits** for each algorithm's options. Every run also has
a wall-clock time limit (120 seconds by default; 0 means unlimited, with a maximum
nonzero setting of 86,400 seconds) and an output item limit (1,000 by default,
configurable from 1 to 100,000). The item limit applies separately to each list
or matrix in the result. Lists exceeding it are explicitly marked as truncated.
The formatted result has an additional fixed budget of approximately 2 million
characters; the report header, input parameters and truncation notice are extra.
Increasing the item limit does not raise that text budget. These limits do not
convert partial searches into completeness claims.

The progress bar shows the current stage and elapsed time. Where the library
does not report completed work, the bar remains indeterminate; item formatting
can show measured progress when the collection size is known. **Stop** and the
time limit terminate the calculation process, including methods without cooperative
cancellation. Only one calculation runs at a time; plotting remains interactive.

Results retain their input curve, parameters and proof/certification status.
Height results include exact enclosure bounds; database decimals are labelled
as approximations. A completed calculation does not imply a proved rank or a
complete Mordell–Weil basis: the library's status and reason are preserved.
Use **Copy** or **Export** to export the displayed text report. **Repeat** reopens
the parameter window with that result's original curve, inputs and limits; it
does not start another calculation until you choose **Run calculation**.
To remove a result, right-click its entry in the history dropdown and choose
**Delete**. This deletes that entry,
even when another result is displayed; stop an active calculation before deleting it. History keeps the last
50 calculations for the current session; use **File → Save** to keep the session
or save individual reports before closing the app.
**Clear** in the Results header removes the entire session history after
confirmation in the same dark dialog. **Clear** is disabled while a calculation
is running or the history is empty. Saved report files are unaffected.
Both side panels start at their minimum widths. Drag either divider to resize its
panel; matching gaps and dividers keep the two sides aligned.
The **Equation** and **Results** panels each have a header chevron that folds the
panel into a narrow, full-height tab on its own side, freeing space for the plot.
The folded and expanded versions share the same top, bottom and outer edge,
including the window margin. Either panel can be folded independently; click its
tab to restore it. Both panels retain their resized width, and Equation keeps its
settings and Coefficients state. Folding uses a short animation when Windows allows
interface animations. A dot on the folded Results tab indicates a running
calculation; the result selection and history remain intact.

## Development

The XAML theme, plot control, immutable calculation snapshots and view models are
separate files. No external charting or UI package is required. The portable model
and view-model sources are linked into the existing test project, so their tests
also run without WPF on non-Windows systems. To run only the Explorer-related tests:

```powershell
dotnet test tests/EllipticCurves.Tests.csproj -c Release --filter "FullyQualifiedName~Explorer"
```

Omit the filter to run all arithmetic and portable Explorer tests. After the
initial dependency restore, these tests use local fixtures and simulated HTTP
responses rather than live LMFDB requests.

The calculation catalog covers public mathematical methods on rational curves,
Fp/Fq curves, finite fields and rational numbers; explicit entries cover computed
properties, returned isomorphism/isogeny maps and LMFDB. Object identity methods,
formatting and duplicate aliases are not separate actions. Coverage tests check
every public method, parameter conversion and representative offline invocation.

The executable's private `--compute-worker` entry point uses redirected UTF-8
streams before creating WPF. `CalculationRunner` owns its child process and
terminates it on cancellation, timeout or app shutdown. The worker also exits if
the host disconnects. A portable test host exercises this protocol, errors,
non-cooperative cancellation and disconnect behavior without opening any UI.

On Windows, also run the compiled-XAML regression check. This standalone host is
not included in `dotnet test EllipticCurves.sln`. It loads the real theme and main
workspace and verifies Tools menu click/focus scrolling, history-menu deletion,
acceptance/rejection of Clear and Reset, the presence of Repeat, both sidebars'
folding, aligned bounds at different window sizes, animation and PNG rendering.
PNG checks cover both embedded views, offsets within the window, fractional layout
sizes, dark backgrounds and preservation of content at the image edges. It also
checks the complex view's shared point selection and compact layout, fitting
presets selected while the real plot is hidden, and viewport preservation across
mode changes and coefficient edits.
It uses an invisible native layout host and shows no application windows:

```powershell
dotnet run --project tests/PresentationHost/PresentationHost.csproj -c Release
```

## Code organization

- `ExplorerInfo` owns application information, Help links, library metadata and the bundled license shared by C# and XAML.
- `Theme.xaml` owns the shared title-bar action and shortcut styles. `TitleBarPopup` handles menu dismissal;
  `MainWindow.CloseTitleBarMenus` closes all menus when a keyboard command runs.
- `SessionFile` owns `.ec` naming, file-picker filters, format identifiers, validation and file I/O.
  `ExplorerSession` contains the saved data and the history limit. `SessionMessages` contains shared session labels.
- `SessionState` defines command availability for both the menu and keyboard shortcuts.
  `SessionStatusViewModel` presents that state; `SessionChanges` compares only data that requires saving.
- `MainWindow.Session` manages session operations, `MainWindow.SessionStatus` observes changes,
  `MainWindow.Sidebars` manages panel layout, and `MainWindow.Plot` handles view navigation and image export.
- `CalculationProtocol` and `CalculationStatus` name the existing worker messages and stored status values.
  `CurvePreset.Classic` and `ClassicEquation` define the initial curve used by the editor, sessions and calculation inputs.
- `ClipboardActions` handles copying and clipboard errors for the equation and calculation reports.
- `BrowserActions` opens external links and reports browser-launch errors. `MainWindow.Help` manages the modeless Help windows.

Shared values belong with their owning feature. Text used in only one place stays
beside that UI or operation. Refactoring these definitions must preserve the session
format, worker messages and displayed text; the tests use literal expectations for those contracts.

The desktop project uses the [Microsoft .NET Desktop SDK settings](https://learn.microsoft.com/en-us/dotnet/core/project-sdk/msbuild-props-desktop).

The shared application logo is [docs/png/ec_logo_v3a.png](../docs/png/ec_logo_v3a.png).
Explorer links it as the embedded WPF resource `ec_logo.png` for the header and
window icon; the library also uses it as its NuGet package icon. After replacing
the PNG, run `./explorer/tools/Update-Icon.ps1` from the repository root to regenerate
[docs/ico/ec_logo.ico](../docs/ico/ec_logo.ico), the multi-size executable icon used
by both Explorer and Console, then rebuild the applications. The files are embedded
during build or packaging; installed applications do not need the `docs` directory.

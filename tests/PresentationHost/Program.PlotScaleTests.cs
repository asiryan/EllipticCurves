using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using EllipticCurves;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;

internal static partial class Program
{
    private static void CheckLargeCurvePlotScales()
    {
        CheckManualCurveAutoFit();
        var search = CurveSearchEngine.Run(new(new()));
        var examples = search.Results.Select(c => (Name: $"elkies-{c.Numerator}-{c.Denominator}", Curve: ElkiesSearchFamily.Create(c.Numerator, c.Denominator).Curve)).ToList();
        examples.Add(("elkies-record-parameter", ElkiesSearchFamily.Create(-9529, 5471).Curve));
        var large = new BigRational(System.Numerics.BigInteger.Pow(10, 80));
        examples.Add(("large-two-components", new EllipticCurveQ(0, 0, 0, -large, 0)));
        examples.Add(("large-one-component", new EllipticCurveQ(0, 0, 0, large, 0)));
        var first = true;
        foreach (var (name, curve) in examples)
        {
            var plot = new CurvePlot { Snapshot = new(curve), ShowGrid = true };
            plot.Measure(new Size(800, 500));
            plot.Arrange(new Rect(0, 0, 800, 500));
            plot.Fit();
            plot.UpdateLayout();
            var saved = plot.CaptureView();
            Require(saved.HorizontalScaleRatio > 1e6, "Large curve coordinates must be fitted with separate axis scales.");
            var data = plot.Snapshot.Plot;
            double left = data.Roots[0], right = data.Roots[^1];
            if (left == right) right += data.CharacteristicScale / 4;
            var visibleWidth = plot.ToScreen(right, 0).X - plot.ToScreen(left, 0).X;
            Require(visibleWidth > 80, "The curve still collapses into a vertical line after Fit.");
            foreach (var x in data.Roots)
                Require(new Rect(plot.RenderSize).Contains(plot.ToScreen(x, data.CenterY(x))), "Fit loses a branch point on a large curve.");
            if (data.Roots.Count == 3)
            {
                double middle = data.Roots[0] / 2 + data.Roots[1] / 2;
                Require(data.TryEvaluate(middle, out var upper, out var lower), "The bounded real component must evaluate.");
                Require(plot.ToScreen(middle, lower).Y - plot.ToScreen(middle, upper).Y > 20, "The oval must have visible height.");
            }

            var anchor = new Point(280, 190);
            var world = plot.ToWorld(anchor);
            plot.ZoomAt(0.5, anchor);
            Require((plot.ToScreen(world.X, world.Y) - anchor).Length < 1e-6, "Unequal-axis zoom moved the world point under the pointer.");
            plot.Pan(61, -37);
            Require((plot.ToScreen(world.X, world.Y) - (anchor + new Vector(61, -37))).Length < 1e-6, "Pan must move by the requested pixels on both axes.");
            plot.RestoreView(saved);
            Require(plot.CaptureView() == saved, "View restore lost the horizontal scale ratio.");
            var roundTrip = plot.ToScreen(plot.ToWorld(anchor).X, plot.ToWorld(anchor).Y);
            Require((roundTrip - anchor).Length < 1e-6, "Screen/world conversion must preserve separate axis scales.");

            var bitmap = PlotImageExporter.Render(plot);
            var pixels = new byte[bitmap.PixelWidth * bitmap.PixelHeight * 4];
            bitmap.CopyPixels(pixels, bitmap.PixelWidth * 4, 0);
            var columns = new HashSet<int>();
            for (int i = 0; i < pixels.Length; i += 4)
                if (pixels[i] == 207 && pixels[i + 1] == 230 && pixels[i + 2] == 99) columns.Add(i / 4 % bitmap.PixelWidth);
            Require(columns.Count > 150, "Exported curve strokes are still confined to a vertical line.");
            var screenshotCurve = curve.A4.Num == -System.Numerics.BigInteger.Parse("4808474278973187287428061457296740047792");
            if ((first || screenshotCurve || name == "large-two-components") && Environment.GetEnvironmentVariable("ELLIPTIC_EXPLORER_PREVIEW_DIRECTORY") is { Length: > 0 } directory)
            {
                System.IO.Directory.CreateDirectory(directory);
                using var output = System.IO.File.Create(System.IO.Path.Combine(directory, name + "-plot.png"));
                var encoder = new PngBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create(bitmap));
                encoder.Save(output);
            }
            first = false;
        }

        var ordinary = new CurvePlot { Snapshot = new(new EllipticCurveQ(0, 0, 0, -1, 0)) };
        ordinary.Measure(new Size(800, 500));
        ordinary.Arrange(new Rect(0, 0, 800, 500));
        ordinary.Fit();
        Require(ordinary.CaptureView().HorizontalScaleRatio == 1, "The classic curve should retain equal axis units.");

        var window = CreateMainWindow();
        try
        {
            SettleSession(window);
            window.ResetHistory();
            var plot = (CurvePlot)window.FindName("Plot");
            var original = plot.CaptureView();
            window.OpenSearchCurve(search.Results[0]);
            SettleSession(window);
            var fitted = plot.CaptureView();
            Require(fitted.HorizontalScaleRatio > 1e6, "Opening an Elkies candidate must fit both coordinate ranges.");
            ApplicationCommands.Undo.Execute(null, window);
            SettleSession(window);
            Require(plot.CaptureView() == original, "Undo must restore the previous axis scales.");
            ApplicationCommands.Redo.Execute(null, window);
            SettleSession(window);
            Require(plot.CaptureView() == fitted, "Redo must restore the fitted large-curve view.");
        }
        finally { window.Close(); }
    }

    private static void CheckManualCurveAutoFit()
    {
        const string equation = "y^2 = x^3 - 1375414938269729933430*x + 20780863582673042643051322404516";
        foreach (var debounce in new[] { false, true })
        {
            var window = CreateMainWindow();
            try
            {
                SettleSession(window);
                var plot = (CurvePlot)window.FindName("Plot");
                plot.Zoom(0.7);
                plot.Pan(35, -15);
                window.ResetHistory();
                var original = plot.CaptureView();
                if (debounce)
                    CompleteSession(async () =>
                    {
                        window.ViewModel.Equation.Text = equation;
                        await window.ViewModel.PendingUpdate;
                        return true;
                    });
                else
                {
                    window.ViewModel.Equation.Text = equation;
                    window.ViewModel.Equation.CommitEdit();
                    window.ViewModel.FlushUpdate();
                }
                // Enter and blur capture history immediately, before queued layout work.
                window.CommitHistory();
                SettleSession(window);
                var fitted = plot.CaptureView();
                Require(fitted == plot.GetResetView() && fitted != original, "A pasted large curve must fit without Reset view, with or without Enter.");
                var branch = plot.Snapshot!.Plot.Roots[0];
                Require(new Rect(plot.RenderSize).Contains(plot.ToScreen(branch, plot.Snapshot.Plot.CenterY(branch))), "The entered curve's branch point is outside the automatic viewport.");
                ApplicationCommands.Undo.Execute(null, window);
                SettleSession(window);
                Require(plot.CaptureView() == original && window.ViewModel.Equation.Text == CurvePreset.ClassicEquation,
                    "One Undo must restore both the equation and its original zoom, without an extra auto-fit step.");
                ApplicationCommands.Redo.Execute(null, window);
                SettleSession(window);
                Require(plot.CaptureView() == fitted, "Redo must preserve the automatic fit.");

                plot.Zoom(0.8);
                plot.Pan(5000, 3000); // A deliberately empty viewport must survive local edits.
                var zoomed = plot.CaptureView();
                window.ViewModel.SimpleCoefficients[0].Text = "-1375414938269729933431";
                window.ViewModel.FlushUpdate();
                SettleSession(window);
                Require(plot.CaptureView() == zoomed, "A small coefficient adjustment must preserve the user's zoom and pan.");
                window.ViewModel.Equation.Text = "y^2 =";
                window.ViewModel.FlushUpdate();
                Require(plot.CaptureView() == zoomed, "Incomplete input must not move the viewport.");
                window.ViewModel.Equation.Text = "y^2 = x^3 - x";
                window.ViewModel.FlushUpdate();
                SettleSession(window);
                Require(plot.CaptureView() == plot.GetResetView() && plot.CaptureView().HorizontalScaleRatio == 1,
                    "Returning from a large curve to a small curve must fit the new scale.");

                window.ViewModel.Equation.Text = "y^2 + 2000000*y = x^3 - x - 1000000000000";
                window.ViewModel.FlushUpdate();
                SettleSession(window);
                Require(plot.CaptureView() == plot.GetResetView() && Math.Abs(plot.CaptureView().CenterY + 1000000) < 1,
                    "A curve translated far outside the former viewport must also be fitted.");
            }
            finally { window.Close(); }
        }

        var hidden = CreateMainWindow();
        try
        {
            SettleSession(hidden);
            var mode = (ComboBox)hidden.FindName("ViewMode");
            var plot = (CurvePlot)hidden.FindName("Plot");
            var original = plot.CaptureView();
            mode.SelectedIndex = 1;
            hidden.ViewModel.Equation.Text = equation;
            hidden.ViewModel.FlushUpdate();
            SettleSession(hidden);
            Require(plot.CaptureView() == original, "A collapsed plot must defer its fit until its layout is available.");
            mode.SelectedIndex = 0;
            SettleSession(hidden);
            Require(plot.CaptureView() == plot.GetResetView() && plot.CaptureView() != original,
                "Returning to Real locus must complete the pending automatic fit.");
        }
        finally { hidden.Close(); }
    }
}

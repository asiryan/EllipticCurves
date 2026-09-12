using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media.Imaging;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

internal static partial class Program
{
    private static void CheckLmfdbImport()
    {
        var formulas = new[]
        {
            new LmfdbCurveFormula("37.a1", 37, "y^2 + y = x^3 - x"),
            new LmfdbCurveFormula("37.b1", 37, "y^2 + y = x^3 + x^2 - 1873*x - 31833"),
            new LmfdbCurveFormula("37.b2", 37, "y^2 + y = x^3 + x^2 - 23*x - 50"),
            new LmfdbCurveFormula("37.b3", 37, "y^2 + y = x^3 + x^2 - 3*x + 1")
        };
        using var model = new LmfdbImportViewModel((_, _, _) => Task.FromResult<IReadOnlyList<LmfdbCurveFormula>>(formulas));
        var dialog = new LmfdbImportWindow(model);
        try
        {
            var button = (Button)dialog.FindName("ImportButton");
            var root = (FrameworkElement)dialog.Content;
            dialog.Content = null;
            root.DataContext = model;
            root.Measure(new Size(758, 678));
            root.Arrange(new Rect(0, 0, 758, 678));
            root.UpdateLayout();
            Require(!button.IsEnabled, "Import must be disabled before selecting a loaded formula.");
            CompleteSession(async () => { await model.SearchAsync(); return true; });
            foreach (var size in new[] { new Size(758, 678), new Size(638, 558) })
            {
                using var source = new HwndSource(new HwndSourceParameters("LMFDB import headless layout")
                    { WindowStyle = 0, Width = (int)size.Width, Height = (int)size.Height }) { RootVisual = root };
                root.Measure(size);
                root.Arrange(new Rect(new Point(), size));
                root.UpdateLayout();
                var list = (ListBox)dialog.FindName("CurveList");
                list.SelectedIndex = 1;
                root.UpdateLayout();
                Require(model.Selected == formulas[1] && button.IsEnabled, "The chosen row must enable formula import.");
                Require(list.ActualHeight >= 60, "The formula list must remain usable at the minimum window size.");
                foreach (var control in new FrameworkElement[] { button, list, (TextBox)dialog.FindName("ConductorInput") })
                    Require(new Rect(root.RenderSize).Contains(control.TransformToAncestor(root).TransformBounds(new Rect(control.RenderSize))),
                        "An import control extends outside the dialog.");
                var bitmap = PlotImageExporter.Render(root);
                Require(bitmap.PixelWidth <= size.Width * 2 && bitmap.PixelHeight <= size.Height * 2,
                    "The import dialog grew beyond the requested window size while rendering.");
                // Optional local visual QA uses the same compiled XAML exercised above.
                if (Environment.GetEnvironmentVariable("ELLIPTIC_EXPLORER_PREVIEW_DIRECTORY") is { Length: > 0 } directory)
                {
                    System.IO.Directory.CreateDirectory(directory);
                    using var output = System.IO.File.Create(System.IO.Path.Combine(directory, $"lmfdb-import-{size.Width}.png"));
                    var encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(bitmap));
                    encoder.Save(output);
                }
            }
            model.Query = "0";
            root.UpdateLayout();
            Require(!button.IsEnabled && !((Button)dialog.FindName("SearchButton")).IsEnabled,
                "Editing to an invalid conductor must disable import and search.");
        }
        finally { dialog.Close(); }

        var menu = new ExplorerMenu();
        ((ExplorerMenuViewModel)menu.DataContext).Group = "LMFDB · internet";
        var importRequests = 0;
        menu.ImportCurveRequested += () => importRequests++;
        var menuRoot = (FrameworkElement)((System.Windows.Controls.Primitives.Popup)menu.FindName("MenuPopup")).Child;
        menuRoot.DataContext = menu.DataContext;
        menuRoot.Measure(new Size(680, 560));
        menuRoot.Arrange(new Rect(new Point(), menuRoot.DesiredSize));
        menuRoot.UpdateLayout();
        Descendants(menuRoot).OfType<Button>().Single(button => ReferenceEquals(button.Tag, ExplorerMenuViewModel.ImportCurve))
            .RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
        Require(importRequests == 1, "The LMFDB menu must open the formula picker exactly once.");

        var window = CreateMainWindow();
        try
        {
            SettleSession(window);
            window.ResetHistory();
            var original = window.CaptureSession();
            var plot = (CurvePlot)window.FindName("Plot");
            var originalView = plot.CaptureView();
            window.ImportLmfdbFormula(formulas[1]);
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == formulas[1].Equation && window.HasUnsavedChanges,
                "Import must replace the equation and mark the session as changed.");
            Require(window.Workbench.Jobs.Count == original.History.Count, "Formula import must not add a database report.");
            var path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), Guid.NewGuid() + ".ec");
            try
            {
                SessionFile.Save(path, window.CaptureSession());
                Require(SessionFile.Load(path).Equation == formulas[1].Equation, "The imported equation must round-trip in .ec files.");
            }
            finally { System.IO.File.Delete(path); }
            var importedView = plot.CaptureView();
            ApplicationCommands.Undo.Execute(null, window);
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == original.Equation && plot.CaptureView() == originalView
                && !window.HasUnsavedChanges && !window.EditHistory.CanUndo,
                "One Undo must restore the original equation, viewport and clean state.");
            ApplicationCommands.Redo.Execute(null, window);
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == formulas[1].Equation && plot.CaptureView() == importedView,
                "Redo must restore the imported formula and fitted plot without a network request.");
        }
        finally { window.Close(); }
    }
}

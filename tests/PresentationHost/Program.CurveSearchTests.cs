using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media.Imaging;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

internal static partial class Program
{
    private static void CheckCurveSearch()
    {
        using var workbench = new WorkbenchViewModel();
        using var model = new CurveSearchViewModel(workbench);
        var dialog = new CurveSearchWindow(model);
        var state = CurveSearchEngine.Run(new(new(-10, 10, 3, 31, 1009, 8, 0, 2)));
        try
        {
            var root = (FrameworkElement)dialog.Content;
            dialog.Content = null;
            root.DataContext = model;
            foreach (var size in new[] { new Size(1078, 798), new Size(878, 678) })
            {
                using var source = new HwndSource(new HwndSourceParameters("Elkies search headless layout")
                    { WindowStyle = 0, Width = (int)size.Width, Height = (int)size.Height }) { RootVisual = root };
                root.Measure(size);
                root.Arrange(new Rect(new Point(), size));
                root.UpdateLayout();
                var start = (Button)dialog.FindName("StartButton");
                Require(start.IsEnabled, "The search must start with built-in data, or recheck a loaded report.");
                model.Load(state);
                root.UpdateLayout();
                var list = (ListBox)dialog.FindName("CandidateList");
                var details = (TextBox)dialog.FindName("CandidateDetails");
                list.SelectedIndex = 1;
                root.UpdateLayout();
                Require(details.Text.Contains(state.Results[1].Equation) && details.Text.Contains("SAVED REPORT"), "The selected report and saved-data qualification must be visible.");
                Require(details.ActualHeight >= 100, "The exact-coordinate report is too short at the minimum size.");
                foreach (var control in new FrameworkElement[] { start, list, details, (Button)dialog.FindName("SaveButton"), (Button)dialog.FindName("OpenCurveButton") })
                    Require(new Rect(root.RenderSize).Contains(control.TransformToAncestor(root).TransformBounds(new Rect(control.RenderSize))), "A search control extends beyond its window.");
                ((Expander)dialog.FindName("AdvancedSettings")).IsExpanded = size.Width < 900;
                root.UpdateLayout();
                if (Environment.GetEnvironmentVariable("ELLIPTIC_EXPLORER_PREVIEW_DIRECTORY") is { Length: > 0 } directory)
                {
                    System.IO.Directory.CreateDirectory(directory);
                    using var output = System.IO.File.Create(System.IO.Path.Combine(directory, $"elkies-search-{size.Width}.png"));
                    var encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(PlotImageExporter.Render(root)));
                    encoder.Save(output);
                }
            }
        }
        finally { dialog.Close(); }

        var menu = new ExplorerMenu();
        ((ExplorerMenuViewModel)menu.DataContext).Search = "Elkies";
        int requested = 0;
        menu.SearchCurvesRequested += () => requested++;
        var menuRoot = (FrameworkElement)((System.Windows.Controls.Primitives.Popup)menu.FindName("MenuPopup")).Child;
        menuRoot.DataContext = menu.DataContext;
        menuRoot.Measure(new Size(680, 560));
        menuRoot.Arrange(new Rect(new Point(), menuRoot.DesiredSize));
        menuRoot.UpdateLayout();
        Descendants(menuRoot).OfType<Button>().Single(b => ReferenceEquals(b.Tag, ExplorerMenuViewModel.SearchCurves)).RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
        Require(requested == 1, "The Elkies action must open its search window.");

        var window = CreateMainWindow();
        try
        {
            SettleSession(window);
            window.ResetHistory();
            var original = window.ViewModel.Equation.Text;
            window.OpenSearchCurve(state.Results[0] with { Equation = "corrupt saved display data" });
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == state.Results[0].Equation, "Opening a candidate must reconstruct the exact curve from the family parameter.");
            ApplicationCommands.Undo.Execute(null, window);
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == original, "Opening a search candidate must be undoable in one step.");
        }
        finally { window.Close(); }
    }
}

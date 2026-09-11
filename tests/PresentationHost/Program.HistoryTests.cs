using System.Windows;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.Computations;

internal static partial class Program
{
    private static void CheckEditHistory()
    {
        CheckResultSelectionHistory();
        var window = CreateMainWindow();
        try
        {
            SettleSession(window);
            window.ResetHistory();
            var menu = (EditMenu)window.FindName("Edit");
            var undo = (Button)menu.FindName("UndoButton");
            var redo = (Button)menu.FindName("RedoButton");
            var plot = (CurvePlot)window.FindName("Plot");
            var torus = (ComplexTorusView)window.FindName("TorusView");
            var mode = (ComboBox)window.FindName("ViewMode");
            var equation = Descendants((FrameworkElement)window.Content).OfType<TextBox>().Single(text =>
                AutomationProperties.GetName(text) == "Curve equation");
            Require(!undo.IsEnabled && !redo.IsEnabled, "A new session must start with disabled Undo and Redo.");
            foreach (var (command, key) in new[] { (ApplicationCommands.Undo, Key.Z), (ApplicationCommands.Redo, Key.Y) })
                Require(window.InputBindings.OfType<KeyBinding>().Any(binding => binding.Command == command
                    && binding.Key == key && binding.Modifiers == ModifierKeys.Control), "An Edit shortcut is missing.");

            var initial = window.ViewModel.Snapshot;
            window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
            window.CommitHistory();
            var edited = window.ViewModel.Snapshot;
            Require(undo.IsEnabled && !redo.IsEnabled && window.HasUnsavedChanges, "Equation edits must enable Undo.");
            CheckEditMenuLayout(window, menu);
            // Tunnelled command bindings must beat TextBox's built-in Undo stack.
            Require(ApplicationCommands.Undo.CanExecute(null, equation), "Undo cannot route from the equation editor.");
            ApplicationCommands.Undo.Execute(null, equation);
            Require(ReferenceEquals(initial, window.ViewModel.Snapshot) && !window.HasUnsavedChanges,
                $"Undo must restore the clean equation and reuse its snapshot: equation={window.ViewModel.Equation.Text}, same={ReferenceEquals(initial, window.ViewModel.Snapshot)}, dirty={window.HasUnsavedChanges}, redo={window.EditHistory.CanRedo}.");
            Require(redo.IsEnabled, "Undo must enable Redo.");
            redo.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(ReferenceEquals(edited, window.ViewModel.Snapshot), "Edit → Redo must restore the computed snapshot.");

            var beforePlot = plot.CaptureView();
            plot.RestoreView(new(12, -4, 20));
            window.CommitHistory();
            ApplicationCommands.Undo.Execute(null, window);
            Require(plot.CaptureView() == beforePlot, "Undo lost the 2D view.");
            ApplicationCommands.Redo.Execute(null, undo);
            Require(plot.CaptureView() == new PlotViewState(12, -4, 20), "Redo cannot route out of the Edit popup.");

            var oldCamera = torus.CaptureCamera();
            mode.SelectedIndex = 1;
            torus.RestoreCamera(new(70, 15, 9));
            window.ViewModel.ShowGrid = false;
            window.ViewModel.ShowPoints = false;
            window.CommitHistory();
            ApplicationCommands.Undo.Execute(null, window);
            Require(mode.SelectedIndex == 0 && torus.CaptureCamera() == oldCamera
                && window.ViewModel.ShowGrid && window.ViewModel.ShowPoints, "Undo must restore graph mode, camera, grid and points.");
            ApplicationCommands.Redo.Execute(null, window);
            Require(mode.SelectedIndex == 1 && torus.CaptureCamera() == new TorusCameraState(70, 15, 9)
                && !window.ViewModel.ShowGrid && !window.ViewModel.ShowPoints, "Redo lost graph state.");

            ApplicationCommands.Undo.Execute(null, window);
            window.ViewModel.Equation.Text = "y^2 =";
            Require(!ApplicationCommands.Redo.CanExecute(null, window), "Pending input must invalidate Redo immediately.");
            window.CommitHistory();
            ApplicationCommands.Undo.Execute(null, window);
            ApplicationCommands.Redo.Execute(null, window);
            Require(window.ViewModel.Equation.Text == "y^2 =" && window.ViewModel.HasIncompleteInput,
                "Redo must restore incomplete text with its last valid graph.");

            var session = ExplorerSession.New() with { History = new() { new(new("Q.TorsionStructure", "y^2 = x^3 - x", new()),
                DateTime.Now, CalculationStatus.Completed, "Done", TimeSpan.FromSeconds(2), 100, "Cached report") } };
            window.RestoreSession(session);
            SettleSession(window);
            Require(!undo.IsEnabled && !redo.IsEnabled, "Open/New must clear both history branches.");
            var report = window.Workbench.Selected!.Report;
            window.Workbench.ClearHistory();
            Require(window.Workbench.Jobs.Count == 0 && undo.IsEnabled, "Clearing reports must be undoable.");
            undo.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(window.Workbench.Selected?.Report == report && !window.Workbench.IsBusy,
                "Undo must restore the report without starting a calculation.");
            Require(!window.HasUnsavedChanges, "Undoing deletion must return an opened session to Saved.");
            redo.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(window.Workbench.Jobs.Count == 0, "Redo must reapply report deletion.");
            window.RestoreSession(ExplorerSession.New());
            Require(!undo.IsEnabled && !redo.IsEnabled, "Replacing a session must discard old history.");

            SettleSession(window);
            window.ResetHistory();
            var originalCoefficient = window.ViewModel.SimpleCoefficients[0].ExactValue;
            foreach (var offset in new[] { 1, 4, 8, 12 }) window.ViewModel.SimpleCoefficients[0].SliderOffset = offset;
            CompleteSession(async () => { await Task.Delay(600); return true; });
            ApplicationCommands.Undo.Execute(null, window);
            Require(window.ViewModel.SimpleCoefficients[0].ExactValue == originalCoefficient && !undo.IsEnabled,
                "Rapid slider edits must coalesce into one step.");
            CompleteSession(async () => { await Task.Delay(600); return true; });
            Require(redo.IsEnabled && !undo.IsEnabled, "Late graph updates must not create edits or clear Redo.");

            CompleteSession(async () =>
            {
                var running = window.Workbench.RunAsync(new("Q.TorsionStructure", "y^2 = x^3 - x", new()));
                Require(window.Workbench.IsBusy && !ApplicationCommands.Undo.CanExecute(null, equation)
                    && !ApplicationCommands.Redo.CanExecute(null, equation), "History must be disabled while a worker is active.");
                window.Workbench.Cancel();
                await running;
                return true;
            });
            var stopped = window.Workbench.Selected!.Report;
            ApplicationCommands.Undo.Execute(null, window);
            Require(window.Workbench.Jobs.Count == 0, "Undo must remove a stopped calculation.");
            ApplicationCommands.Redo.Execute(null, window);
            Require(window.Workbench.Selected?.Report == stopped && !window.Workbench.IsBusy,
                "Redo must restore a stopped calculation without restarting it.");
        }
        finally { window.Close(); }
    }

    private static void CheckResultSelectionHistory()
    {
        var window = CreateMainWindow();
        try
        {
            window.RestoreSession(ExplorerSession.New() with
            {
                History = Enumerable.Range(0, 3).Select(index => new CalculationSession(
                    new("Q.TorsionStructure", "y^2 = x^3 - x", new()), DateTime.Now.AddMinutes(-index),
                    CalculationStatus.Completed, "Done", TimeSpan.FromSeconds(2), 100, "Cached report " + index)).ToList()
            });
            SettleSession(window);
            window.ResetHistory();
            var mode = (ComboBox)window.FindName("ViewMode");
            var torus = (ComplexTorusView)window.FindName("TorusView");
            var picker = Descendants((ResultsPanel)window.FindName("Results")).OfType<ComboBox>().Single();
            var report = Descendants((ResultsPanel)window.FindName("Results")).OfType<TextBox>().Single();
            var snapshot = window.ViewModel.Snapshot;
            var camera = new TorusCameraState(70, 15, 9);
            mode.SelectedIndex = 1;
            torus.RestoreCamera(camera);
            // Use the real two-way binding, before the graph edit debounce fires.
            picker.SelectedIndex = 1;
            picker.SelectedIndex = 2;
            ApplicationCommands.Undo.Execute(null, picker);
            RequireResult(1);
            ApplicationCommands.Undo.Execute(null, picker);
            RequireResult(0);
            ApplicationCommands.Undo.Execute(null, picker);
            Require(mode.SelectedIndex == 0 && !window.EditHistory.CanUndo,
                "Graph edits must remain undoable after undoing each result selection.");
            ApplicationCommands.Redo.Execute(null, picker);
            RequireResult(0);
            ApplicationCommands.Redo.Execute(null, picker);
            RequireResult(1);
            ApplicationCommands.Redo.Execute(null, picker);
            RequireResult(2);
            ApplicationCommands.Undo.Execute(null, picker);
            picker.SelectedIndex = 1;
            Require(window.EditHistory.CanRedo, "Selecting the current result must preserve Redo.");
            picker.SelectedIndex = 0;
            Require(!window.EditHistory.CanRedo, "Selecting a different result must create a new history branch.");
            ApplicationCommands.Undo.Execute(null, picker);
            RequireResult(1);
            CompleteSession(async () => { await Task.Delay(600); return true; });
            Require(window.EditHistory.CanRedo, "Deferred bindings must not create phantom selection edits after Undo.");

            foreach (var action in new Action[]
            {
                () => window.Workbench.Delete(window.Workbench.Selected),
                () => window.Workbench.Delete(window.Workbench.Jobs[0]),
                window.Workbench.ClearHistory
            })
            {
                window.ResetHistory();
                action();
                ApplicationCommands.Undo.Execute(null, picker);
                RequireResult(1);
                Require(window.Workbench.Jobs.Count == 3 && !window.EditHistory.CanUndo,
                    "Deleting or clearing results and automatically selecting a replacement must form one undo step.");
            }

            void RequireResult(int index)
            {
                SettleSession(window);
                Require(window.Workbench.Selected == window.Workbench.Jobs[index] && picker.SelectedIndex == index
                    && report.Text == window.Workbench.Selected.Report,
                    "Undo/Redo must restore the selected result in both the picker and the displayed report.");
                Require(mode.SelectedIndex == 1 && torus.CaptureCamera() == camera,
                    "Undo/Redo of result selection must not jump to Real locus or change the torus camera.");
                Require(ReferenceEquals(snapshot, window.ViewModel.Snapshot) && !window.Workbench.IsBusy
                    && !window.HasUnsavedChanges, "Browsing results must reuse cached data and leave the session saved.");
            }
        }
        finally { window.Close(); }
    }

    private static void CheckEditMenuLayout(MainWindow window, EditMenu menu)
    {
        var popup = (Popup)menu.FindName("MenuPopup");
        var root = (FrameworkElement)popup.Child;
        root.Measure(new Size(220, double.PositiveInfinity));
        root.Arrange(new Rect(root.DesiredSize));
        root.UpdateLayout();
        var buttons = Descendants(root).OfType<Button>().ToArray();
        Require(buttons.Select(button => button.Content).SequenceEqual(new[] { "Undo", "Redo" }),
            "Edit must contain Undo and Redo in order.");
        var rights = new List<double>();
        foreach (var (button, key) in buttons.Zip(new[] { "Z", "Y" }))
        {
            var texts = Descendants(button).OfType<TextBlock>().ToArray();
            var label = texts.Single(text => text.Text == (string)button.Content);
            var hint = texts.Single(text => text.Text == "Ctrl + " + key);
            Rect Bounds(FrameworkElement element) => element.TransformToAncestor(root).TransformBounds(new Rect(element.RenderSize));
            Require(Bounds(hint).Left > Bounds(label).Right && new Rect(root.RenderSize).Contains(Bounds(hint)),
                "Edit shortcut hints overlap their labels or are clipped.");
            Require(AutomationProperties.GetAcceleratorKey(button) == "Ctrl+" + key, "Edit accelerator accessibility is missing.");
            rights.Add(Bounds(hint).Right);
        }
        Require(rights.Max() - rights.Min() < 1, "Edit shortcuts must align at the right edge.");
        if (Environment.GetEnvironmentVariable("EC_EDIT_PREVIEW") is not { Length: > 0 } preview) return;
        var content = (FrameworkElement)window.Content;
        System.Windows.Data.BindingOperations.ClearBinding(popup, Popup.IsOpenProperty);
        var toggle = (ToggleButton)menu.FindName("Toggle");
        toggle.IsChecked = true;
        content.UpdateLayout();
        var drawing = new DrawingVisual();
        using (var context = drawing.RenderOpen())
        {
            context.DrawRectangle(new VisualBrush(content), null, new Rect(content.RenderSize));
            var location = menu.TranslatePoint(new Point(0, menu.ActualHeight), content);
            context.DrawRectangle(new VisualBrush(root), null, new Rect(location, root.RenderSize));
        }
        var bitmap = new RenderTargetBitmap((int)content.ActualWidth, (int)content.ActualHeight, 96, 96, PixelFormats.Pbgra32);
        bitmap.Render(drawing);
        var encoder = new PngBitmapEncoder();
        encoder.Frames.Add(BitmapFrame.Create(bitmap));
        using var stream = System.IO.File.Create(preview);
        encoder.Save(stream);
        toggle.IsChecked = false;
    }
}

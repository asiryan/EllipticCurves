using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Threading;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;

internal static partial class Program
{
    private static Task<bool> StartSession(Func<Task<bool>> action)
    {
        var previous = SynchronizationContext.Current;
        SynchronizationContext.SetSynchronizationContext(new DispatcherSynchronizationContext(Dispatcher.CurrentDispatcher));
        try { return action(); }
        finally { SynchronizationContext.SetSynchronizationContext(previous); }
    }

    private static bool CompleteSession(Func<Task<bool>> action) => WaitForSession(StartSession(action));

    private static bool WaitForSession(Task<bool> task)
    {
        if (!task.IsCompleted)
        {
            var dispatcher = Dispatcher.CurrentDispatcher;
            var frame = new DispatcherFrame();
            _ = task.ContinueWith(_ => dispatcher.BeginInvoke(new Action(() => frame.Continue = false)), TaskScheduler.Default);
            Dispatcher.PushFrame(frame);
        }
        return task.GetAwaiter().GetResult();
    }

    private static void SettleSession(MainWindow window, double width = 1438, double height = 918)
    {
        var root = (FrameworkElement)window.Content;
        root.Measure(new Size(width, height));
        root.Arrange(new Rect(0, 0, width, height));
        root.UpdateLayout();
        Dispatcher.CurrentDispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(root.UpdateLayout));
    }

    private static void CheckSessionLifecycle()
    {
        const string edited = "y^2 = x^3 - 5*x + 3";
        var input = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        var output = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        SessionFile.Save(input, ExplorerSession.New() with { Equation = "y^2 = x^3 + 7", Preset = null });
        try
        {
            foreach (var action in new[] { "New", "Open", "Close" })
            foreach (var option in new[] { "Cancel", "Discard", "Save", "Cancel save", "Failed save" })
            {
                File.Delete(output);
                var choice = option == "Cancel" ? SaveChangesChoice.Cancel
                    : option == "Discard" ? SaveChangesChoice.Discard : SaveChangesChoice.Save;
                var prompts = 0;
                var errors = 0;
                var openPickers = 0;
                var window = new MainWindow(null, new SessionDialogs(name =>
                    {
                        Require(name == "untitled.ec", "An untitled session should suggest a default file name.");
                        prompts++;
                        return new(choice, "renamed-study");
                    },
                    suggested =>
                    {
                        Require(suggested == "renamed-study.ec", "The edited session name did not reach the save picker with its extension.");
                        return option == "Cancel save" ? null : option == "Failed save" ? output + ".missing/session.ec" : output;
                    },
                    () =>
                    {
                        Require(prompts == 1, "Open must confirm unsaved changes before showing the file picker.");
                        openPickers++;
                        return input;
                    }, (_, _) => errors++));
                var closed = false;
                window.Closed += (_, _) => closed = true;
                try
                {
                    Require(!window.HasUnsavedChanges, "A fresh untouched session was marked as modified.");
                    window.ViewModel.Equation.Text = edited;
                    window.ViewModel.Step.Text = "0.2";
                    window.ViewModel.ShowGrid = false;
                    window.ViewModel.ShowPoints = false;
                    window.Workbench.Jobs.Add(CalculationJobViewModel.FromSession(new CalculationSession(
                        new CalculationRequest("Q.TorsionStructure", edited, new()), DateTime.Now,
                        "Completed", "Done", TimeSpan.FromSeconds(1), 100, "A result that must be saved")));
                    window.Workbench.Selected = window.Workbench.Jobs[0];
                    ((CurvePlot)window.FindName("Plot")).Zoom(0.7);
                    Require(window.HasUnsavedChanges, "Editing a session did not require a save warning.");
                    bool CloseWindow() { CompleteSession(() => { window.Close(); return window.PendingSessionOperation; }); return closed; }
                    var proceeded = action == "New" ? CompleteSession(window.NewSessionAsync) : action == "Open" ? CompleteSession(window.OpenSessionAsync) : CloseWindow();
                    var expected = option is "Discard" or "Save";
                    Require(proceeded == expected && prompts == 1,
                        $"{action} / {option} did not obey the unsaved-changes choice.");
                    Require(errors == (option == "Failed save" ? 1 : 0), "A save error was hidden or unexpectedly raised.");
                    Require(openPickers == (action == "Open" && expected ? 1 : 0),
                        "Cancelling or failing to save must not open the file picker.");
                    if (!proceeded)
                        Require(!closed && window.ViewModel.Equation.Text == edited && window.HasUnsavedChanges,
                            "Cancelling or failing to save discarded the current session.");
                    else if (!closed)
                    {
                        SettleSession(window);
                        Require(!window.HasUnsavedChanges, "A new or opened session was immediately marked as modified.");
                        Require(window.ViewModel.Equation.Text == (action == "New" ? "y^2 = x^3 - x" : "y^2 = x^3 + 7"),
                            "New/Open restored the wrong curve.");
                        if (action == "New")
                            Require(!window.Workbench.HasResults && window.ViewModel.ShowGrid && window.ViewModel.ShowPoints
                                && window.ViewModel.Step.Text == "0.01", "New did not reset the session defaults.");
                    }
                    if (option == "Save") Require(SessionFile.Load(output).Equation == edited, "Save did not preserve the old session before leaving.");
                }
                finally { choice = SaveChangesChoice.Discard; if (!closed) window.Close(); }
            }

            var answer = SaveChangesChoice.Save;
            var confirmations = 0;
            var sameFile = new MainWindow(null, new SessionDialogs(_ => { confirmations++; return new(answer); },
                _ => input, () => input, (title, message) => throw new Exception(title + ": " + message)));
            try
            {
                Require(CompleteSession(sameFile.OpenSessionAsync) && confirmations == 0, "Opening from a clean session unnecessarily prompted.");
                SettleSession(sameFile);
                Require(!sameFile.HasUnsavedChanges, "The loaded baseline did not survive layout.");
                sameFile.ViewModel.Equation.Text = edited;
                Require(CompleteSession(sameFile.OpenSessionAsync), "Save followed by opening the same file failed.");
                Require(sameFile.ViewModel.Equation.Text == edited && !sameFile.HasUnsavedChanges,
                    "Opening the same file restored the stale copy from before Save.");
                ((CurvePlot)sameFile.FindName("Plot")).Zoom(0.8);
                Require(!sameFile.HasUnsavedChanges && !CompleteSession(() => sameFile.TrySaveSessionAsync()) && !sameFile.HasUnsavedChanges,
                    "Graph navigation must leave the session saved and Save unavailable.");
                sameFile.ViewModel.Equation.Text = "y^2 =";
                Require(sameFile.HasUnsavedChanges && sameFile.ViewModel.Equation.Error == "",
                    "Checking unsaved edits must not commit or silently ignore incomplete input.");
            }
            finally { answer = SaveChangesChoice.Discard; sameFile.Close(); }
        }
        finally { File.Delete(input); File.Delete(output); }
    }

    private static void CheckSessionNavigation()
    {
        var input = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        var output = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        try
        {
            foreach (var complex in new[] { false, true })
            foreach (var action in new[] { "New", "Open", "Close" })
            {
                File.Delete(output);
                // Legacy files containing camera positions must also reopen fitted.
                var legacy = ExplorerSession.New() with
                {
                    Equation = "y^2 + x*y + y = x^3 - 5*x + 3", Preset = null, SliderOffsets = Array.Empty<int>(),
                    ComplexView = complex, Plot = new(1e18, -1e18, 1e8), FitRealViewWhenShown = false,
                    TorusCamera = new(-71, 23, 8.9),
                    History = new()
                    {
                        new(new("Q.TorsionStructure", "y^2 = x^3 - x", new()), new DateTime(2026, 9, 11, 21, 42, 0),
                            "Completed", "Done", TimeSpan.FromSeconds(2), 100, "Newest result"),
                        new(new("Q.TorsionStructure", "y^2 = x^3 + x", new()), new DateTime(2026, 9, 11, 21, 40, 0),
                            "Completed", "Done", TimeSpan.FromSeconds(1), 100, "Older result")
                    }
                };
                SessionFile.Save(input, legacy);
                var prompts = 0;
                var window = new MainWindow(null, new SessionDialogs(_ => { prompts++; return new(SaveChangesChoice.Cancel); },
                    _ => output, () => input, (title, message) => throw new Exception(title + ": " + message)));
                var closed = false;
                window.Closed += (_, _) => closed = true;
                try
                {
                    Require(CompleteSession(window.OpenSessionAsync), "Could not open the legacy graph session.");
                    SettleSession(window);
                    var plot = (CurvePlot)window.FindName("Plot");
                    var torus = (ComplexTorusView)window.FindName("TorusView");
                    if (!complex)
                    {
                        var opened = plot.CaptureView();
                        plot.Fit();
                        Require(opened == plot.CaptureView() && opened != legacy.Plot,
                            "An old file must open at Reset view, not at its stored pan and zoom.");
                    }
                    Require(torus.CaptureCamera() == TorusCameraState.Default && !window.HasUnsavedChanges,
                        "Opening/resetting a graph must leave a clean session and a default torus camera.");
                    var mode = (ComboBox)window.FindName("ViewMode");
                    foreach (var index in new[] { 1, 0, 1 })
                    {
                        mode.SelectedIndex = index;
                        SettleSession(window);
                        Require(!window.HasUnsavedChanges && window.SessionStatus.Status == "Saved"
                            && !window.SessionStatus.DisplayName.EndsWith(" *"),
                            "Switching between 2D and 3D must not mark a saved session as modified.");
                    }
                    window.ViewModel.ShowGrid = false;
                    window.ViewModel.ShowPoints = false;
                    SettleSession(window);
                    Require(!window.HasUnsavedChanges, "Toggling grid and samples must not mark the session dirty.");
                    window.ViewModel.ShowPoints = true;
                    CompleteSession(async () =>
                    {
                        await window.ViewModel.PendingSamples;
                        return true;
                    });
                    SettleSession(window);
                    // Drain sample bindings before starting the hidden view's background model explicitly.
                    CompleteSession(async () =>
                    {
                        torus.Model.Update(window.ViewModel.Snapshot.Curve, window.ViewModel.Samples, true);
                        await torus.Model.PendingUpdate;
                        return true;
                    });
                    SettleSession(window);
                    Require(!window.HasUnsavedChanges && window.SessionStatus.Status == "Saved",
                        "Background torus preparation must leave the session saved.");
                    torus.Model.SelectedPoint = torus.Model.Points.First(point => !point.Point.IsInfinity);
                    var selectedPoint = torus.Model.SessionSelection;
                    var coefficients = (Expander)window.FindName("CoefficientsExpander");
                    coefficients.IsExpanded = true;
                    ((ColumnDefinition)window.FindName("EquationColumn")).Width = new GridLength(320);
                    ((ColumnDefinition)window.FindName("ResultsColumn")).Width = new GridLength(410);
                    var equationScroll = (ScrollViewer)window.FindName("EquationScroll");
                    equationScroll.Height = 280;
                    SettleSession(window);
                    equationScroll.ScrollToBottom();
                    SettleSession(window);
                    Require(equationScroll.VerticalOffset > 0 && !window.HasUnsavedChanges
                        && window.SessionStatus.Status == "Saved",
                        "Display settings, point selection and panel layout must not count as session edits.");
                    Require(window.Workbench.Selected == window.Workbench.Jobs[0],
                        "Opening a session must select the newest result.");
                    var results = (ResultsPanel)window.FindName("Results");
                    var historyPicker = Descendants(results).OfType<ComboBox>().Single();
                    foreach (var index in new[] { 1, 0, 1 })
                    {
                        historyPicker.SelectedIndex = index;
                        Require(window.Workbench.Selected == window.Workbench.Jobs[index] && !window.HasUnsavedChanges,
                            "Browsing the Results list must update the report without marking the session dirty.");
                    }
                    var fitted = plot.GetResetView();
                    plot.RestoreView(new(1e18, -1e18, 1e8));
                    plot.Zoom(0.7);
                    torus.RestoreCamera(legacy.TorusCamera);
                    torus.Zoom(0.8);
                    var navigatedPlot = plot.CaptureView();
                    var navigatedTorus = torus.CaptureCamera();
                    Require(!window.HasUnsavedChanges, "Panning, zooming or rotating a graph triggered unsaved changes.");
                    Require(CompleteSession(() => window.TrySaveSessionAsync(saveAs: true)), "An explicit save after graph navigation failed.");
                    var saved = SessionFile.Load(output);
                    Require(System.Text.Json.JsonSerializer.Serialize(saved.History) == System.Text.Json.JsonSerializer.Serialize(legacy.History)
                        && window.Workbench.Selected == window.Workbench.Jobs[1],
                        "Saving must preserve every result, default to the newest and leave the displayed report unchanged.");
                    Require(saved.Plot == fitted && saved.FitRealViewWhenShown && saved.TorusCamera == TorusCameraState.Default && saved.ComplexView,
                        "The file must contain Reset view regardless of the user's current graph navigation.");
                    Require(!saved.ShowGrid && saved.ShowPoints && saved.SelectedTorusPoint == selectedPoint
                        && saved.CoefficientsExpanded && saved.EquationPanel.Width == 320 && saved.ResultsPanel.Width == 410
                        && saved.EquationScrollOffset == equationScroll.VerticalOffset,
                        "Save as must preserve visual settings and layout even when they did not mark the session dirty.");
                    Require(plot.CaptureView() == navigatedPlot && torus.CaptureCamera() == navigatedTorus && !window.HasUnsavedChanges,
                        "Saving must not move the visible graph or leave the session dirty.");

                    var restored = CreateMainWindow();
                    try
                    {
                        restored.RestoreSession(saved);
                        SettleSession(restored, 1120, 760);
                        var restoredVisuals = restored.CaptureSession();
                        Require(restoredVisuals.ComplexView && !restoredVisuals.ShowGrid && restoredVisuals.ShowPoints
                            && restoredVisuals.SelectedTorusPoint == selectedPoint && restoredVisuals.CoefficientsExpanded
                            && !restored.HasUnsavedChanges,
                            "Opening must restore the saved visual settings without marking the session dirty.");
                        Require(restored.Workbench.Jobs.Count == 2 && restored.Workbench.Selected == restored.Workbench.Jobs[0]
                            && restored.Workbench.Selected.Result == "Newest result",
                            "Reopening must restore the whole history and display its newest result.");
                        var restoredPlot = (CurvePlot)restored.FindName("Plot");
                        // A file saved with the real plot hidden must fit it on first display.
                        ((ComboBox)restored.FindName("ViewMode")).SelectedIndex = 0;
                        SettleSession(restored, 1120, 760);
                        var opened = restoredPlot.CaptureView();
                        restoredPlot.Fit();
                        Require(opened == restoredPlot.CaptureView(), "Opening must fit the curve to the new window size, including a deferred real view.");
                        Require(((ComplexTorusView)restored.FindName("TorusView")).CaptureCamera() == TorusCameraState.Default,
                            "The saved session restored a navigated torus camera.");
                    }
                    finally { restored.Close(); }

                    mode.SelectedIndex = 0;
                    window.ViewModel.ShowGrid = true;
                    torus.Model.SelectedPoint = torus.Model.Points[0];
                    window.ViewModel.ShowPoints = false;
                    coefficients.IsExpanded = false;
                    ((Button)window.FindName("EquationCollapseButton")).RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                    SettleSession(window);
                    bool CloseWindow() { CompleteSession(() => { window.Close(); return window.PendingSessionOperation; }); return closed; }
                    Require(action == "New" ? CompleteSession(window.NewSessionAsync) : action == "Open" ? CompleteSession(window.OpenSessionAsync) : CloseWindow(),
                        $"{action} was blocked after graph or result navigation only.");
                    Require(prompts == 0, $"{action} asked to save temporary graph or result navigation.");
                }
                finally { if (!closed) window.Close(); }
            }
        }
        finally { File.Delete(input); File.Delete(output); }
    }

    private static void CheckSessionSaveName()
    {
        var input = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        var output = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        SessionFile.Save(input, ExplorerSession.New());
        try
        {
            foreach (var name in new[] { "curve study", "curve study.ec", "curve study.EC" })
            {
                File.Delete(output);
                var choice = SaveChangesChoice.Save;
                string? openPath = input, savePath = null, suggestedPath = null;
                var window = new MainWindow(null, new SessionDialogs(currentName =>
                    {
                        Require(currentName == Path.GetFileName(input), "The confirmation must show the current file name.");
                        return new(choice, name);
                    }, suggested => { suggestedPath = suggested; return savePath; }, () => openPath,
                    (title, message) => throw new Exception(title + ": " + message)));
                try
                {
                    Require(CompleteSession(window.OpenSessionAsync), "Could not load the named session fixture.");
                    SettleSession(window);
                    window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
                    openPath = null;
                    Require(!CompleteSession(window.OpenSessionAsync) && window.HasUnsavedChanges, "Cancelling the renamed save must preserve unsaved edits.");
                    var expectedName = name.EndsWith(".ec", StringComparison.OrdinalIgnoreCase) ? name : name + ".ec";
                    Require(suggestedPath == Path.Combine(Path.GetDirectoryName(input)!, expectedName),
                        "Renaming must retain the current folder and include exactly one .ec extension.");
                    Require(!CompleteSession(() => window.TrySaveSessionAsync(saveAs: true)) && suggestedPath == input,
                        "Cancelling the save picker must not rename the current session.");
                    savePath = output;
                    Require(!CompleteSession(window.OpenSessionAsync) && !window.HasUnsavedChanges && SessionFile.Load(output).Equation == "y^2 = x^3 + 7",
                        "The final path chosen in the save picker must receive the session even if Open is then cancelled.");
                    savePath = null;
                    Require(!CompleteSession(() => window.TrySaveSessionAsync(saveAs: true)) && suggestedPath == output,
                        "A subsequent Save must use the actual saved path, not the proposed name.");
                }
                finally { choice = SaveChangesChoice.Discard; window.Close(); }
            }
        }
        finally { File.Delete(input); File.Delete(output); }
    }

    private static void CheckSessionShortcuts()
    {
        var prompts = 0;
        var saves = 0;
        var choice = SaveChangesChoice.Cancel;
        var window = new MainWindow(null, new SessionDialogs(_ => { prompts++; return new(choice); },
            _ => { saves++; return null; }, () => throw new Exception("Cancelled Open must not show a file picker."),
            (title, message) => throw new Exception(title + ": " + message)));
        try
        {
            SettleSession(window);
            var session = (SessionMenu)window.FindName("Session");
            var explorer = (ExplorerMenu)window.FindName("Explorer");
            var sessionToggle = (ToggleButton)session.FindName("Toggle");
            var explorerToggle = (ToggleButton)explorer.FindName("Toggle");
            foreach (var menu in new UserControl[] { session, explorer })
                System.Windows.Data.BindingOperations.ClearBinding((Popup)menu.FindName("MenuPopup"), Popup.IsOpenProperty);
            var equation = Descendants((FrameworkElement)window.Content).OfType<TextBox>().Single(text =>
                System.Windows.Automation.AutomationProperties.GetName(text) == "Curve equation");
            var origins = new IInputElement[] { window, equation, (Button)session.FindName("NewButton"), (TextBox)explorer.FindName("SearchBox") };
            foreach (var (command, key) in new[] { (ApplicationCommands.New, Key.N), (ApplicationCommands.Open, Key.O), (ApplicationCommands.Save, Key.S), (ApplicationCommands.SaveAs, Key.S) })
            {
                var modifiers = command == ApplicationCommands.SaveAs ? ModifierKeys.Control | ModifierKeys.Shift : ModifierKeys.Control;
                var saving = command == ApplicationCommands.Save || command == ApplicationCommands.SaveAs;
                Require(window.InputBindings.OfType<KeyBinding>().Any(binding =>
                    binding.Command == command && binding.Key == key && binding.Modifiers == modifiers),
                    $"The Ctrl+{key} session shortcut is missing.");
                foreach (var origin in origins)
                {
                    window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
                    var oldPrompts = prompts;
                    var oldSaves = saves;
                    if (command == ApplicationCommands.Save)
                    {
                        Require(!command.CanExecute(null, origin), $"Save must be disabled for new sessions from {origin.GetType().Name}.");
                        command.Execute(null, origin);
                        Require(prompts == oldPrompts && saves == oldSaves, "Disabled Save opened a dialog.");
                        continue;
                    }
                    sessionToggle.IsChecked = explorerToggle.IsChecked = true;
                    Require(command.CanExecute(null, origin), $"{command.Name} cannot route from {origin.GetType().Name}.");
                    CompleteSession(() => { command.Execute(null, origin); return window.PendingSessionOperation; });
                    Require(sessionToggle.IsChecked == false && explorerToggle.IsChecked == false,
                        "A session shortcut must dismiss either open title-bar menu.");
                    Require(prompts - oldPrompts == (saving ? 0 : 1)
                        && saves - oldSaves == (saving ? 1 : 0),
                        "A session shortcut bypassed confirmation or executed more than once.");
                    Require(window.HasUnsavedChanges && window.ViewModel.Equation.Text == "y^2 = x^3 + 7",
                        "Cancelling a shortcut's dialog discarded the current session.");
                }
            }
            window.Workbench.Dispose();
            Require(!ApplicationCommands.New.CanExecute(null, window) && !ApplicationCommands.Open.CanExecute(null, window)
                && !ApplicationCommands.Save.CanExecute(null, window) && ApplicationCommands.SaveAs.CanExecute(null, window),
                "Session shortcuts must respect the same CanRun restriction as the menu.");
        }
        finally { choice = SaveChangesChoice.Discard; window.Close(); }
    }

    private static void CheckSaveChangesDialog()
    {
        foreach (var action in new[] { "ConfirmButton", "DiscardButton", "CancelButton", "CloseButton" })
        {
            var dialog = ConfirmationWindow.CreateSaveChangesDialog("untitled.ec");
            try
            {
                var fileName = (TextBox)dialog.FindName("SessionFileName");
                var save = (Button)dialog.FindName("ConfirmButton");
                var nameError = (TextBlock)dialog.FindName("SessionNameError");
                Require(fileName.Text == "untitled.ec" && !fileName.IsReadOnly && save.IsEnabled,
                    "The session file name must be prefilled and editable.");
                foreach (var invalid in new[] { "", "   ", "../session", "bad:name", "session.", "a\\b" })
                {
                    fileName.Text = invalid;
                    Require(!save.IsEnabled && nameError.Visibility == Visibility.Visible,
                        "An empty or invalid file name must not enable Save.");
                }
                fileName.Text = "  curve-study.ec  ";
                Require(save.IsEnabled && nameError.Visibility == Visibility.Collapsed, "Correcting the file name did not re-enable Save.");
                var root = (FrameworkElement)dialog.Content;
                root.Measure(new Size(dialog.Width, double.PositiveInfinity));
                root.Arrange(new Rect(root.DesiredSize));
                root.UpdateLayout();
                var discard = (Button)dialog.FindName("DiscardButton");
                var cancel = (Button)dialog.FindName("CancelButton");
                Require(Equals(save.Content, "Save") && Equals(discard.Content, "Discard") && Equals(cancel.Content, "Cancel"),
                    "The unsaved-changes actions must be in English.");
                Require(discard.Visibility == Visibility.Visible && cancel.IsDefault && cancel.IsCancel
                    && FocusManager.GetFocusedElement(dialog) == root, "Unsaved changes must default to cancellation without highlighting an action.");
                Rect Bounds(FrameworkElement element) => element.TransformToAncestor(root).TransformBounds(new Rect(element.RenderSize));
                Require(Bounds(save).Right < Bounds(discard).Left && Bounds(discard).Right < Bounds(cancel).Left
                    && new Rect(root.RenderSize).Contains(Bounds(cancel)), "Save/Discard/Cancel buttons overlap or extend outside the dialog.");
                if (action == "CancelButton")
                {
                    var bitmap = new RenderTargetBitmap((int)root.ActualWidth, (int)root.ActualHeight, 96, 96, PixelFormats.Pbgra32);
                    bitmap.Render(root);
                    var encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(bitmap));
                    using var stream = File.Create(Path.Combine(AppContext.BaseDirectory, "confirmation-session.png"));
                    encoder.Save(stream);
                }
                ((Button)dialog.FindName(action)).RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                Require(dialog.SaveChoice == (action == "ConfirmButton" ? SaveChangesChoice.Save
                    : action == "DiscardButton" ? SaveChangesChoice.Discard : SaveChangesChoice.Cancel), "The save confirmation returned the wrong choice.");
                Require(dialog.SaveResult.FileName == (action == "ConfirmButton" ? "curve-study.ec" : null),
                    "Only Save should return the edited, trimmed file name.");
            }
            finally { dialog.Close(); }
        }
    }
}

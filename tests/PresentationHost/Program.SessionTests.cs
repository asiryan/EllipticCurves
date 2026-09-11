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
    private static void SettleSession(MainWindow window)
    {
        var root = (FrameworkElement)window.Content;
        root.Measure(new Size(1438, 918));
        root.Arrange(new Rect(0, 0, 1438, 918));
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
                var choice = option == "Cancel" ? SaveChangesChoice.Cancel
                    : option == "Discard" ? SaveChangesChoice.Discard : SaveChangesChoice.Save;
                var prompts = 0;
                var errors = 0;
                var openPickers = 0;
                var window = new MainWindow(null, new SessionDialogs(name =>
                    {
                        Require(name == "session.ec", "An untitled session should suggest a default file name.");
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
                    bool CloseWindow() { window.Close(); return closed; }
                    var proceeded = action == "New" ? window.NewSession() : action == "Open" ? window.OpenSession() : CloseWindow();
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
                Require(sameFile.OpenSession() && confirmations == 0, "Opening from a clean session unnecessarily prompted.");
                SettleSession(sameFile);
                Require(!sameFile.HasUnsavedChanges, "The loaded baseline did not survive layout.");
                sameFile.ViewModel.Equation.Text = edited;
                Require(sameFile.OpenSession(), "Save followed by opening the same file failed.");
                Require(sameFile.ViewModel.Equation.Text == edited && !sameFile.HasUnsavedChanges,
                    "Opening the same file restored the stale copy from before Save.");
                ((CurvePlot)sameFile.FindName("Plot")).Zoom(0.8);
                Require(sameFile.HasUnsavedChanges && sameFile.TrySaveSession() && !sameFile.HasUnsavedChanges,
                    "Viewport changes or successful saves did not update the dirty state.");
                sameFile.ViewModel.Equation.Text = "y^2 =";
                Require(sameFile.HasUnsavedChanges && sameFile.ViewModel.Equation.Error == "",
                    "Checking unsaved edits must not commit or silently ignore incomplete input.");
            }
            finally { answer = SaveChangesChoice.Discard; sameFile.Close(); }
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
                    Require(window.OpenSession(), "Could not load the named session fixture.");
                    SettleSession(window);
                    window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
                    openPath = null;
                    Require(!window.OpenSession() && window.HasUnsavedChanges, "Cancelling the renamed save must preserve unsaved edits.");
                    var expectedName = name.EndsWith(".ec", StringComparison.OrdinalIgnoreCase) ? name : name + ".ec";
                    Require(suggestedPath == Path.Combine(Path.GetDirectoryName(input)!, expectedName),
                        "Renaming must retain the current folder and include exactly one .ec extension.");
                    Require(!window.TrySaveSession() && suggestedPath == input,
                        "Cancelling the save picker must not rename the current session.");
                    savePath = output;
                    Require(!window.OpenSession() && !window.HasUnsavedChanges && SessionFile.Load(output).Equation == "y^2 = x^3 + 7",
                        "The final path chosen in the save picker must receive the session even if Open is then cancelled.");
                    savePath = null;
                    Require(!window.TrySaveSession() && suggestedPath == output,
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
            foreach (var (command, key) in new[] { (ApplicationCommands.New, Key.N), (ApplicationCommands.Open, Key.O), (ApplicationCommands.Save, Key.S) })
            {
                Require(window.InputBindings.OfType<KeyBinding>().Any(binding =>
                    binding.Command == command && binding.Key == key && binding.Modifiers == ModifierKeys.Control),
                    $"The Ctrl+{key} session shortcut is missing.");
                foreach (var origin in origins)
                {
                    window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
                    var oldPrompts = prompts;
                    var oldSaves = saves;
                    sessionToggle.IsChecked = explorerToggle.IsChecked = true;
                    Require(command.CanExecute(null, origin), $"{command.Name} cannot route from {origin.GetType().Name}.");
                    command.Execute(null, origin);
                    Require(sessionToggle.IsChecked == false && explorerToggle.IsChecked == false,
                        "A session shortcut must dismiss either open title-bar menu.");
                    Require(prompts - oldPrompts == (command == ApplicationCommands.Save ? 0 : 1)
                        && saves - oldSaves == (command == ApplicationCommands.Save ? 1 : 0),
                        "A session shortcut bypassed confirmation or executed more than once.");
                    Require(window.HasUnsavedChanges && window.ViewModel.Equation.Text == "y^2 = x^3 + 7",
                        "Cancelling a shortcut's dialog discarded the current session.");
                }
            }
            window.Workbench.Dispose();
            Require(!ApplicationCommands.New.CanExecute(null, window) && !ApplicationCommands.Open.CanExecute(null, window)
                && ApplicationCommands.Save.CanExecute(null, window), "Session shortcuts must respect the same CanRun restriction as the menu.");
        }
        finally { choice = SaveChangesChoice.Discard; window.Close(); }
    }

    private static void CheckSaveChangesDialog()
    {
        foreach (var action in new[] { "ConfirmButton", "DiscardButton", "CancelButton", "CloseButton" })
        {
            var dialog = ConfirmationWindow.CreateSaveChangesDialog("session.ec");
            try
            {
                var fileName = (TextBox)dialog.FindName("SessionFileName");
                var save = (Button)dialog.FindName("ConfirmButton");
                var nameError = (TextBlock)dialog.FindName("SessionNameError");
                Require(fileName.Text == "session.ec" && !fileName.IsReadOnly && save.IsEnabled,
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

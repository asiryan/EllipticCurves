using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;

internal static partial class Program
{
    private static bool ExecuteSessionCommand(MainWindow window, RoutedCommand command) => CompleteSession(() =>
    {
        command.Execute(null, window);
        return window.PendingSessionOperation;
    });

    private static void CheckSessionSaveStatus()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-session-status-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var first = Path.Combine(directory, "curve-study.ec");
        var second = Path.Combine(directory, "curve-study-copy.ec");
        string? choosePath = first;
        var pickers = 0;
        var errors = 0;
        var writes = 0;
        var window = new MainWindow(null, new SessionDialogs(_ => new(SaveChangesChoice.Discard),
            _ => { pickers++; return choosePath; }, () => first, (_, _) => errors++, (path, snapshot, createCopy) =>
            {
                writes++;
                return Task.Run(() => SessionFile.Save(path, snapshot, createCopy));
            }));
        try
        {
            void CheckStatus(string text, string name)
            {
                SettleSession(window);
                Require(window.SessionStatus.Status == text && window.SessionStatus.DisplayName == name,
                    $"Expected session status '{name} / {text}', got '{window.SessionStatus.DisplayName} / {window.SessionStatus.Status}'.");
                Require(((TextBlock)window.FindName("SessionFileName")).Text == name
                    && ((TextBlock)window.FindName("SessionSaveStatus")).Text == text
                    && window.Title == name + " — Elliptic Curves · Explorer",
                    "The visible file name, save status and native title must update together.");
            }

            void CheckSaveAvailability(bool canSave)
            {
                var menu = (SessionMenu)window.FindName("Session");
                Require(ApplicationCommands.Save.CanExecute(null, window) == canSave
                    && window.SessionStatus.CanSave == canSave && ((Button)menu.FindName("SaveButton")).IsEnabled == canSave,
                    "Save and Ctrl+S must require an existing file and unsaved changes or a failed save.");
                Require(ApplicationCommands.SaveAs.CanExecute(null, window) && ((Button)menu.FindName("SaveAsButton")).IsEnabled
                    && ((Button)menu.FindName("ExitButton")).IsEnabled, "Save as and Exit must remain available for new sessions.");
                if (!canSave)
                {
                    var previousWrites = writes;
                    var previousPickers = pickers;
                    ApplicationCommands.Save.Execute(null, window);
                    Require(!CompleteSession(() => window.TrySaveSessionAsync())
                        && writes == previousWrites && pickers == previousPickers,
                        "Disabled Save must not write a file or open a picker, including through the menu handler.");
                }
            }

            CheckStatus("New session", "untitled.ec");
            CheckSaveAvailability(false);
            ApplicationCommands.Save.Execute(null, window);
            Require(!CompleteSession(() => window.TrySaveSessionAsync()) && pickers == 0 && !File.Exists(first),
                "Save without an existing file must not open a picker or write a file.");
            window.ViewModel.ShowGrid = false;
            CheckStatus("New session", "untitled.ec");
            window.ViewModel.ShowGrid = true;
            CheckStatus("New session", "untitled.ec");
            window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
            CheckStatus("Unsaved changes", "untitled.ec *");
            ((ComboBox)window.FindName("ViewMode")).SelectedIndex = 1;
            window.ViewModel.ShowPoints = false;
            CheckStatus("Unsaved changes", "untitled.ec *");
            CheckSaveAvailability(false);
            Require(ExecuteSessionCommand(window, ApplicationCommands.SaveAs) && pickers == 1,
                "The first save must choose a file through Ctrl+Shift+S.");
            CheckStatus("Saved", "curve-study.ec");
            CheckSaveAvailability(false);
            Require(window.SessionStatus.FileLocation == first && SessionFile.Load(first).Equation == "y^2 = x^3 + 7",
                "Saving must display the real file path and write the edited session.");
            window.ViewModel.Step.Text = "0.2";
            CheckStatus("Unsaved changes", "curve-study.ec *");
            CheckSaveAvailability(true);
            window.ViewModel.Step.Text = "0.01";
            CheckStatus("Saved", "curve-study.ec");
            CheckSaveAvailability(false);
            window.ViewModel.Step.Text = "0.2";
            CheckStatus("Unsaved changes", "curve-study.ec *");
            CheckSaveAvailability(true);
            Require(ExecuteSessionCommand(window, ApplicationCommands.Save) && pickers == 1
                && SessionFile.Load(first).SliderStep == "0.2", "Ctrl+S must overwrite the current file without a picker.");
            CheckStatus("Saved", "curve-study.ec");
            CheckSaveAvailability(false);
            window.ViewModel.ShowGrid = false;
            window.ViewModel.ShowPoints = false;
            ((ComboBox)window.FindName("ViewMode")).SelectedIndex = 1;
            CheckStatus("Saved", "curve-study.ec");
            CheckSaveAvailability(false);
            Require(SessionFile.Load(first).ShowGrid, "Disabled Save must leave the stored visual settings unchanged.");
            var original = File.ReadAllBytes(first);

            choosePath = null;
            Require(!ExecuteSessionCommand(window, ApplicationCommands.SaveAs) && pickers == 2,
                "Ctrl+Shift+S must show a cancellable Save as dialog.");
            CheckStatus("Saved", "curve-study.ec");
            Require(window.SessionStatus.FileLocation == first, "Cancelling Save as changed the current path.");
            choosePath = second;
            Require(ExecuteSessionCommand(window, ApplicationCommands.SaveAs) && pickers == 3,
                "Save as must create and select another session file.");
            CheckStatus("Saved", "curve-study-copy.ec");
            CheckSaveAvailability(false);
            var visual = SessionFile.Load(second);
            Require(visual.ComplexView && !visual.ShowGrid && !visual.ShowPoints,
                "Save as must still write visual changes while the session shows Saved.");
            window.ViewModel.Equation.Text = "y^2 = x^3 - 5*x + 3";
            Require(ExecuteSessionCommand(window, ApplicationCommands.Save) && pickers == 3
                && SessionFile.Load(second).Equation == "y^2 = x^3 - 5*x + 3" && File.ReadAllBytes(first).SequenceEqual(original),
                "After Save as, Ctrl+S must update only the newly selected file.");

            choosePath = Path.Combine(directory, "missing", "failed.ec");
            Require(!ExecuteSessionCommand(window, ApplicationCommands.SaveAs) && errors == 1, "A failed save must report its error.");
            CheckStatus("Save failed", "curve-study-copy.ec");
            CheckSaveAvailability(true);
            RenderSessionHeader(window, "session-save-failed.png");
            Require(window.SessionStatus.FileLocation == second && !File.Exists(choosePath),
                "A failed Save as must retain the actual saved file's identity.");
            Require(ExecuteSessionCommand(window, ApplicationCommands.Save), "Retrying Ctrl+S to the existing file failed.");
            CheckStatus("Saved", "curve-study-copy.ec");
            CheckSaveAvailability(false);
            RenderSessionHeader(window, "session-saved.png");

            Require(CompleteSession(window.NewSessionAsync), "New failed after saving.");
            CheckStatus("New session", "untitled.ec");
            CheckSaveAvailability(false);
            choosePath = null;
            Require(!ExecuteSessionCommand(window, ApplicationCommands.SaveAs), "Cancelling the initial Save as must not create a session file.");
            CheckStatus("New session", "untitled.ec");
            CheckSaveAvailability(false);
            Require(CompleteSession(window.OpenSessionAsync), "Opening an existing file failed.");
            CheckStatus("Saved", "curve-study.ec");
            CheckSaveAvailability(false);
            Require(window.SessionStatus.FileLocation == first, "Open did not update the file path tooltip.");

            var root = (FrameworkElement)window.Content;
            var repository = Descendants(root).OfType<Button>().Single(button =>
                System.Windows.Automation.AutomationProperties.GetName(button) == "Open Elliptic Curves on GitHub");
            Require(Equals(repository.ToolTip, "https://github.com/asiryan/EllipticCurves"),
                "The shared repository address must resolve in the compiled title bar.");
            var indicator = (FrameworkElement)window.FindName("SessionFileIndicator");
            var badge = (FrameworkElement)window.FindName("SessionStatusBadge");
            var fileName = (TextBlock)window.FindName("SessionFileName");
            Rect Bounds(FrameworkElement element) => element.TransformToAncestor(root).TransformBounds(new Rect(element.RenderSize));
            var explorer = (FrameworkElement)window.FindName("Explorer");
            var minimize = Descendants(root).OfType<Button>().Single(button => System.Windows.Automation.AutomationProperties.GetName(button) == "Minimize");
            foreach (var width in new[] { 1120, 1440, 1920 })
            foreach (var name in new[] { "untitled.ec", new string('x', 180) + ".ec" })
            {
                window.SessionStatus.Update(new(Path.Combine(directory, name), true, false, false, false, true));
                root.Measure(new Size(width, 760));
                root.Arrange(new Rect(0, 0, width, 760));
                root.UpdateLayout();
                var bounds = Bounds(indicator);
                Require(Math.Abs(bounds.Left + bounds.Width / 2 - root.ActualWidth / 2) <= 1,
                    "The session indicator must be centered on the whole window.");
                Require(bounds.Left > Bounds(explorer).Right && bounds.Right < Bounds(minimize).Left,
                    "The session indicator overlaps the menus or caption buttons.");
                Require(Bounds(fileName).Left - Bounds(badge).Right >= 11 && Bounds(fileName).Right <= bounds.Right + 1,
                    "The status must precede the file name with a fixed gap and no overflow.");
                Require(fileName.TextTrimming == TextTrimming.CharacterEllipsis && fileName.ActualWidth > 0
                    && Equals(indicator.ToolTip, window.SessionStatus.FileLocation), "Long file names must truncate while exposing the full path.");
                if (width == 1120) RenderSessionHeader(window, name.Length > 20 ? "session-header-long-name.png" : "session-header-centered.png");
            }
        }
        finally
        {
            window.Close();
            File.Delete(first);
            File.Delete(second);
            Directory.Delete(directory);
        }
    }

    private static void CheckSessionSaveAsCopies()
    {
        var directory = Path.Combine(Path.GetTempPath(), "ec-session-copy-" + Guid.NewGuid());
        Directory.CreateDirectory(directory);
        var original = Path.Combine(directory, "session.ec");
        var occupied = Path.Combine(directory, "session(1).ec");
        SessionFile.Save(original, ExplorerSession.New());
        File.WriteAllText(occupied, "Existing file");
        var originalBytes = File.ReadAllBytes(original);
        string? chosen = original;
        var window = new MainWindow(null, new SessionDialogs(_ => new(SaveChangesChoice.Discard),
            _ => chosen, () => original, (title, message) => throw new Exception(title + ": " + message)));
        try
        {
            var dialog = MainWindow.CreateSessionSaveDialog(original);
            Require(dialog.FileName == "session(2).ec" && dialog.InitialDirectory == directory,
                "Save as must immediately display an available file name, with its folder configured separately.");
            Require(!dialog.OverwritePrompt && dialog.AddExtension && dialog.DefaultExt == "ec",
                "Save as must retain the session extension and avoid overwrite prompts when creating numbered copies.");
            Require(MainWindow.CreateSessionSaveDialog(occupied).FileName == "session(2).ec",
                "The picker suggestion must continue an existing numeric suffix.");
            var unused = MainWindow.CreateSessionSaveDialog(Path.Combine(directory, "another study.ec"));
            Require(unused.FileName == "another study.ec" && unused.InitialDirectory == directory,
                "An unused file name must not receive a suffix.");
            var untitled = MainWindow.CreateSessionSaveDialog(null);
            var renamed = MainWindow.CreateSessionSaveDialog("renamed study.ec");
            Require(untitled.FileName == "untitled.ec" && untitled.InitialDirectory == ""
                && renamed.FileName == "renamed study.ec" && renamed.InitialDirectory == "",
                "New sessions must preserve the suggested name and the native picker's default folder.");
            Require(Directory.GetFiles(directory).Length == 2,
                "Preparing a Save as dialog must not reserve a file that would remain after cancellation.");
            Require(CompleteSession(window.OpenSessionAsync), "Could not open the original session.");
            window.ViewModel.ShowGrid = false;
            Require(ExecuteSessionCommand(window, ApplicationCommands.SaveAs), "Save as with an occupied name failed.");
            var copy = Path.Combine(directory, "session(2).ec");
            SettleSession(window);
            Require(window.SessionStatus.FileLocation == copy && window.SessionStatus.FileName == "session(2).ec"
                && window.SessionStatus.Status == "Saved" && !SessionFile.Load(copy).ShowGrid,
                "Save as must display the actual numbered file and save its visual settings.");
            Require(File.ReadAllBytes(original).SequenceEqual(originalBytes) && File.ReadAllText(occupied) == "Existing file",
                "Save as overwrote an existing file.");

            window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
            Require(ExecuteSessionCommand(window, ApplicationCommands.Save) && SessionFile.Load(copy).Equation == "y^2 = x^3 + 7"
                && Directory.GetFiles(directory).Length == 3, "Ctrl+S must update the current numbered file without making another copy.");
            chosen = copy;
            Require(ExecuteSessionCommand(window, ApplicationCommands.SaveAs), "Saving a copy of a numbered file failed.");
            var nextCopy = Path.Combine(directory, "session(3).ec");
            Require(window.SessionStatus.FileLocation == nextCopy && SessionFile.Load(nextCopy).Equation == "y^2 = x^3 + 7",
                "Save as must increment the existing suffix and select that new file.");
            chosen = null;
            Require(!ExecuteSessionCommand(window, ApplicationCommands.SaveAs) && window.SessionStatus.FileLocation == nextCopy
                && Directory.GetFiles(directory).Length == 4, "Cancelling Save as must not create a copy or change the current path.");
            Require(CompleteSession(window.NewSessionAsync), "Could not start a new session after saving a copy.");
            chosen = original;
            Require(ExecuteSessionCommand(window, ApplicationCommands.SaveAs)
                && window.SessionStatus.FileName == "session(4).ec" && File.ReadAllBytes(original).SequenceEqual(originalBytes),
                "The first save must also avoid existing names.");
        }
        finally
        {
            window.Close();
            foreach (var file in Directory.GetFiles(directory)) File.Delete(file);
            Directory.Delete(directory);
        }
    }

    private static void CheckSessionSavingConcurrency()
    {
        foreach (var action in new[] { "Save", "New", "Open", "Close" })
        {
            var path = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
            var gate = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
            var choice = SaveChangesChoice.Save;
            var writes = 0;
            var window = new MainWindow(null, new SessionDialogs(_ => new(choice), _ => path,
                () => throw new Exception("New edits must cancel the pending Open before its picker."),
                (title, message) => throw new Exception(title + ": " + message), async (destination, snapshot, createCopy) =>
                {
                    writes++;
                    await gate.Task;
                    return SessionFile.Save(destination, snapshot, createCopy);
                }));
            var closed = false;
            window.Closed += (_, _) => closed = true;
            try
            {
                SettleSession(window);
                window.ViewModel.Equation.Text = "y^2 = x^3 + 7";
                var saving = StartSession(() => action switch
                {
                    "New" => window.NewSessionAsync(),
                    "Open" => window.OpenSessionAsync(),
                    "Close" => CloseWindow(),
                    _ => window.TrySaveSessionAsync(saveAs: true)
                });
                Task<bool> CloseWindow() { window.Close(); return window.PendingSessionOperation; }
                SettleSession(window);
                Require(!saving.IsCompleted && window.SessionStatus.IsSaving && window.SessionStatus.Status == "Saving…"
                    && ((TextBlock)window.FindName("SessionSaveStatus")).Text == "Saving…", "Saving must be visible while writing asynchronously.");
                if (action == "Save") RenderSessionHeader(window, "session-saving.png");
                Require(!ApplicationCommands.Save.CanExecute(null, window) && !ApplicationCommands.SaveAs.CanExecute(null, window)
                    && !ApplicationCommands.New.CanExecute(null, window) && !window.SessionStatus.CanSave,
                    "Session commands must not overlap an active write.");
                Require(!CompleteSession(() => window.TrySaveSessionAsync(saveAs: true)) && writes == 1, "Repeated Save as started a concurrent write.");
                window.Close();
                Require(!closed, "The window closed before its save completed.");
                window.ViewModel.Equation.Text = "y^2 = x^3 - 5*x + 3";
                gate.SetResult();
                Require(WaitForSession(saving) == (action == "Save"), "Edits during Save must cancel a pending New/Open/Close.");
                SettleSession(window);
                Require(!closed && window.ViewModel.Equation.Text == "y^2 = x^3 - 5*x + 3" && window.HasUnsavedChanges
                    && window.SessionStatus.Status == "Unsaved changes" && window.SessionStatus.DisplayName.EndsWith(" *")
                    && SessionFile.Load(path).Equation == "y^2 = x^3 + 7", "Saving marked later edits clean or discarded them.");
                if (action == "Save") RenderSessionHeader(window, "session-unsaved.png");
            }
            finally
            {
                gate.TrySetResult();
                if (!window.PendingSessionOperation.IsCompleted) WaitForSession(window.PendingSessionOperation);
                choice = SaveChangesChoice.Discard;
                if (!closed) window.Close();
                File.Delete(path);
            }
        }
    }

    private static void RenderSessionHeader(MainWindow window, string fileName)
    {
        var root = (FrameworkElement)window.Content;
        root.UpdateLayout();
        var bitmap = new RenderTargetBitmap((int)root.ActualWidth, 44, 96, 96, PixelFormats.Pbgra32);
        bitmap.Render(root);
        var encoder = new PngBitmapEncoder();
        encoder.Frames.Add(BitmapFrame.Create(bitmap));
        using var stream = File.Create(Path.Combine(AppContext.BaseDirectory, fileName));
        encoder.Save(stream);
    }
}

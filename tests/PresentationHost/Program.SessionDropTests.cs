using System.IO;
using System.Reflection;
using System.Windows;
using System.Windows.Controls;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.Computations;

internal static partial class Program
{
    private static DragEventArgs RaiseSessionDrag(UIElement target, IDataObject data, RoutedEvent routedEvent,
        DragDropEffects allowed = DragDropEffects.Copy)
    {
        // WPF creates these arguments internally for OLE drops; raise the actual routed event in tests.
        var args = (DragEventArgs)Activator.CreateInstance(typeof(DragEventArgs),
            BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic, null,
            new object[] { data, DragDropKeyStates.None, allowed, target, new Point() }, null)!;
        args.RoutedEvent = routedEvent;
        target.RaiseEvent(args);
        return args;
    }

    private static DragEventArgs RaiseSessionFileDrag(MainWindow window, string path, RoutedEvent routedEvent)
    {
        SettleSession(window);
        var equation = Descendants((DependencyObject)window.Content).OfType<TextBox>().Single(box => ReferenceEquals(box.DataContext, window.ViewModel.Equation));
        var data = new DataObject(DataFormats.FileDrop, new[] { path });
        data.SetText(path);
        return RaiseSessionDrag(equation, data, routedEvent);
    }

    private static Task<bool> DropSessionFile(MainWindow window, string path)
    {
        RaiseSessionFileDrag(window, path, DragDrop.PreviewDropEvent);
        return window.PendingSessionOperation;
    }

    private static void CheckSessionFileDrop()
    {
        var input = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".EC");
        var invalid = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        var directory = Path.Combine(Path.GetTempPath(), Guid.NewGuid() + ".ec");
        const string equation = "y^2 = x^3 + 7";
        var errors = 0;
        var prompts = 0;
        SessionFile.Save(input, ExplorerSession.New() with
        {
            Equation = equation,
            History = new() { new(new CalculationRequest("Q.TorsionStructure", equation, new()),
                DateTime.Now, "Completed", "Done", TimeSpan.Zero, 100, "Saved report") }
        });
        File.WriteAllText(invalid, "invalid session");
        Directory.CreateDirectory(directory);
        var window = new MainWindow(null, new SessionDialogs(_ => { prompts++; return new(SaveChangesChoice.Discard); },
            _ => throw new Exception("Dropping into a clean session must not prompt to save."),
            () => throw new Exception("A file drop must not open the file picker."), (_, _) => errors++));
        try
        {
            SettleSession(window);
            var editor = Descendants((DependencyObject)window.Content).OfType<TextBox>().Single(box => ReferenceEquals(box.DataContext, window.ViewModel.Equation));
            Require(window.AllowDrop && editor.AllowDrop, "The window and equation editor must allow file drops.");
            foreach (var routedEvent in new[] { DragDrop.PreviewDragEnterEvent, DragDrop.PreviewDragOverEvent })
            {
                var args = RaiseSessionFileDrag(window, input, routedEvent);
                Require(args.Handled && args.Effects == DragDropEffects.Copy, "A single .ec file must show the copy cursor.");
                Require(window.ViewModel.Equation.Text == CurvePreset.ClassicEquation && prompts == 0,
                    "Hovering a file must not open it or prompt to save.");
            }
            Require(CompleteSession(() => DropSessionFile(window, input)), "Dropping a session over the equation editor failed.");
            SettleSession(window);
            Require(window.ViewModel.Equation.Text == equation && editor.Text == equation && !window.HasUnsavedChanges
                && window.SessionStatus.FileName == Path.GetFileName(input) && window.Workbench.Selected?.Result == "Saved report"
                && errors == 0 && prompts == 0, "A dropped session did not restore its equation, reports and file name cleanly.");

            var plot = (UIElement)window.FindName("Plot");
            foreach (var paths in new[] { Array.Empty<string>(), new[] { input, invalid }, new[] { input + ".txt" },
                new[] { directory }, new[] { input + ".missing.ec" } })
            foreach (var routedEvent in new[] { DragDrop.PreviewDragOverEvent, DragDrop.PreviewDropEvent })
            {
                var pending = window.PendingSessionOperation;
                var args = RaiseSessionDrag(plot, new DataObject(DataFormats.FileDrop, paths), routedEvent);
                Require(args.Handled && args.Effects == DragDropEffects.None && ReferenceEquals(pending, window.PendingSessionOperation),
                    "Unsupported file drops must be rejected without starting a session operation.");
            }
            var move = RaiseSessionDrag(editor, new DataObject(DataFormats.FileDrop, new[] { input }),
                DragDrop.PreviewDropEvent, DragDropEffects.Move);
            Require(move.Handled && move.Effects == DragDropEffects.None && File.Exists(input), "Opening a session must not move its file.");
            var text = RaiseSessionDrag(window, new DataObject(DataFormats.UnicodeText, equation), DragDrop.PreviewDragOverEvent);
            Require(!text.Handled, "Session file handling must leave normal text dragging alone.");
            Require(!CompleteSession(() => DropSessionFile(window, invalid)) && errors == 1
                && window.ViewModel.Equation.Text == equation && !window.HasUnsavedChanges,
                "An invalid session must report an error and preserve the current document.");

            CompleteSession(async () =>
            {
                var running = window.Workbench.RunAsync(new("Q.TorsionStructure", equation, new()));
                try
                {
                    Require(window.Workbench.IsBusy, "The calculation must be active for the drop guard check.");
                    var pending = window.PendingSessionOperation;
                    Require(RaiseSessionFileDrag(window, input, DragDrop.PreviewDragOverEvent).Effects == DragDropEffects.None,
                        "File drops must be disabled during calculations.");
                    RaiseSessionFileDrag(window, input, DragDrop.PreviewDropEvent);
                    Require(ReferenceEquals(pending, window.PendingSessionOperation) && prompts == 0,
                        "A file drop replaced a session with a running calculation.");
                }
                finally { window.Workbench.Cancel(); await running; }
                return true;
            });
        }
        finally
        {
            window.Close();
            File.Delete(input);
            File.Delete(invalid);
            Directory.Delete(directory);
        }
    }
}

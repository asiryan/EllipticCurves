using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

internal sealed record SessionDialogs(Func<string, SaveChangesResult> ConfirmUnsaved,
    Func<string?, string?> ChooseSavePath, Func<string?> ChooseOpenPath, Action<string, string> ShowError,
    Func<string, ExplorerSession, Task>? WriteSession = null);

public partial class MainWindow
{
    private string? sessionPath;
    private int sessionRestoreVersion;
    private int cleanSessionVersion;
    private ExplorerSession? cleanSession;
    private SessionDialogs sessionDialogs = null!;
    private bool sessionActionInProgress, isSavingSession, sessionSaveFailed, allowSessionClose;
    internal Task<bool> PendingSessionOperation { get; private set; } = Task.FromResult(true);

    private void SessionCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
    {
        e.CanExecute = !sessionActionInProgress && (e.Command == ApplicationCommands.Save
            ? sessionPath != null : e.Command == ApplicationCommands.SaveAs || Workbench.CanRun);
        e.Handled = true;
    }

    private async void SessionCommandExecuted(object sender, ExecutedRoutedEventArgs e)
    {
        Session.Close();
        Explorer.Close();
        e.Handled = true;
        if (e.Command == ApplicationCommands.New) await NewSessionAsync();
        else if (e.Command == ApplicationCommands.Open) await OpenSessionAsync();
        else if (e.Command == ApplicationCommands.Save) await TrySaveSessionAsync();
        else if (e.Command == ApplicationCommands.SaveAs) await TrySaveSessionAsync(saveAs: true);
    }

    private void InitializeSession(SessionDialogs? dialogs)
    {
        sessionDialogs = dialogs ?? new SessionDialogs(
            name => ConfirmationWindow.AskToSaveChanges(this, name),
            path =>
            {
                var dialog = new SaveFileDialog
                {
                    Filter = "Explorer session (*.ec)|*.ec", DefaultExt = ".ec", AddExtension = true,
                    FileName = path ?? "session.ec", Title = "Save session as", OverwritePrompt = true
                };
                return dialog.ShowDialog(this) == true ? dialog.FileName : null;
            },
            () =>
            {
                var dialog = new OpenFileDialog
                {
                    Filter = "Explorer session (*.ec)|*.ec", DefaultExt = ".ec",
                    Title = "Open session", CheckFileExists = true, Multiselect = false
                };
                return dialog.ShowDialog(this) == true ? dialog.FileName : null;
            },
            (title, message) => ConfirmationWindow.ShowMessage(this, title, message));
        MarkSessionClean();
        InitializeSessionStatus();
    }

    internal bool HasUnsavedChanges => Workbench.IsBusy || HasSessionEdits;
    private bool HasSessionEdits => cleanSession != null &&
        (ViewModel.HasIncompleteInput || !ViewModel.Step.IsValid
         || !SessionChanges.Equal(cleanSession, ReadSession()));

    private void MarkSessionClean()
    {
        cleanSession = ReadSession();
        cleanSessionVersion++;
        QueueSessionStatusRefresh();
    }

    private async Task<bool> ConfirmSessionChangeAsync()
    {
        if (!HasUnsavedChanges) return true;
        var result = sessionDialogs.ConfirmUnsaved(sessionPath == null ? "session.ec" : Path.GetFileName(sessionPath));
        return result.Choice switch
        {
            SaveChangesChoice.Save => await SaveSessionCoreAsync(false, result.FileName) && !HasSessionEdits,
            SaveChangesChoice.Discard => true,
            _ => false
        };
    }

    private Task<bool> RunSessionOperation(Func<Task<bool>> action)
    {
        if (sessionActionInProgress || sessionClosed) return Task.FromResult(false);
        return PendingSessionOperation = RunAsync();

        async Task<bool> RunAsync()
        {
            sessionActionInProgress = true;
            RefreshSessionStatus();
            try { return await action(); }
            finally
            {
                sessionActionInProgress = false;
                RefreshSessionStatus();
            }
        }
    }

    internal Task<bool> NewSessionAsync() => RunSessionOperation(async () =>
    {
        if (!Workbench.CanRun || !await ConfirmSessionChangeAsync()) return false;
        RestoreSession(ExplorerSession.New());
        sessionPath = null;
        sessionSaveFailed = false;
        return true;
    });

    internal Task<bool> OpenSessionAsync() => RunSessionOperation(async () =>
    {
        if (!Workbench.CanRun || !await ConfirmSessionChangeAsync()) return false;
        var path = sessionDialogs.ChooseOpenPath();
        if (path == null) return false;
        try
        {
            var saved = SessionFile.Load(path);
            RestoreSession(saved);
            sessionPath = Path.GetFullPath(path);
            sessionSaveFailed = false;
            return true;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            sessionDialogs.ShowError("Open session", "The session could not be opened.\n\n" + error.Message);
            return false;
        }
    });

    internal Task<bool> TrySaveSessionAsync(bool saveAs = false, string? fileName = null) =>
        RunSessionOperation(() => saveAs || sessionPath != null
            ? SaveSessionCoreAsync(saveAs, fileName) : Task.FromResult(false));

    private async Task<bool> SaveSessionCoreAsync(bool saveAs, string? fileName)
    {
        ExplorerSession saved;
        try { saved = CaptureSession(); }
        catch (InvalidOperationException error)
        {
            sessionSaveFailed = true;
            RefreshSessionStatus();
            sessionDialogs.ShowError("Save session", error.Message);
            return false;
        }
        var suggestedPath = sessionPath;
        if (!string.IsNullOrWhiteSpace(fileName))
        {
            fileName = fileName.Trim();
            if (!fileName.EndsWith(".ec", StringComparison.OrdinalIgnoreCase)) fileName += ".ec";
            suggestedPath = sessionPath == null ? fileName : Path.Combine(Path.GetDirectoryName(sessionPath)!, fileName);
        }
        var choosePath = saveAs || sessionPath == null || !string.Equals(suggestedPath, sessionPath, StringComparison.OrdinalIgnoreCase);
        var path = choosePath ? sessionDialogs.ChooseSavePath(suggestedPath) : sessionPath;
        if (path == null) return false;
        try
        {
            // Capture on the UI thread, then write an immutable snapshot off it.
            // Later edits remain dirty and must not be lost to a pending New/Open/Close.
            saved = CaptureSession();
            isSavingSession = true;
            sessionSaveFailed = false;
            RefreshSessionStatus();
            if (sessionDialogs.WriteSession is { } write) await write(path, saved);
            else await Task.Run(() => SessionFile.Save(path, saved));
            sessionPath = Path.GetFullPath(path);
            cleanSession = saved;
            cleanSessionVersion++;
            return true;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            sessionSaveFailed = true;
            isSavingSession = false;
            RefreshSessionStatus();
            sessionDialogs.ShowError("Save session", "The session could not be saved.\n\n" + error.Message);
            return false;
        }
        finally
        {
            isSavingSession = false;
            RefreshSessionStatus();
        }
    }

    private async Task<bool> FinishSessionCloseAsync(Task<bool> confirmation)
    {
        if (!await confirmation) return false;
        allowSessionClose = true;
        try { Close(); }
        finally { allowSessionClose = false; }
        return sessionClosed;
    }

    internal ExplorerSession CaptureSession()
    {
        ViewModel.Equation.CommitEdit();
        ViewModel.Step.CommitEdit();
        ViewModel.FlushUpdate();
        if (ViewModel.HasIncompleteInput || !ViewModel.Step.IsValid)
            throw new InvalidOperationException("Finish the curve equation and enter a valid slider step before saving the session.");
        // Save the reset view without changing the graph the user is exploring.
        return ReadSession() with
        {
            Plot = Plot.GetResetView(), TorusCamera = TorusCameraState.Default, FitRealViewWhenShown = true
        };
    }

    private ExplorerSession ReadSession()
    {
        return new ExplorerSession
        {
            Equation = ViewModel.Equation.Text,
            SliderStep = ViewModel.Step.Text,
            Preset = ViewModel.SelectedPreset?.Name,
            SliderOffsets = ViewModel.ActiveCoefficients.Select(coefficient => (int)coefficient.SliderOffset).ToArray(),
            ShowGrid = ViewModel.ShowGrid,
            ShowPoints = ViewModel.ShowPoints,
            ComplexView = IsComplexView,
            FitRealViewWhenShown = realViewResetPending,
            CoefficientsExpanded = CoefficientsExpander.IsExpanded,
            EquationScrollOffset = EquationScroll.VerticalOffset,
            TorusScrollOffset = TorusView.ScrollOffset,
            SelectedTorusPoint = TorusView.Model.SessionSelection,
            Plot = Plot.CaptureView(),
            TorusCamera = TorusView.CaptureCamera(),
            EquationPanel = CaptureSidebar(equationSidebar, EquationColumn),
            ResultsPanel = CaptureSidebar(resultsSidebar, ResultsColumn),
            History = Workbench.Jobs.Select(job => job.CaptureSession()).ToList()
        };
    }

    private static SidebarSession CaptureSidebar(SidebarState state, ColumnDefinition column) => new(state.IsVisible,
        Math.Max(state.MinimumWidth, state.IsVisible && !state.IsAnimating ? column.Width.Value : state.ExpandedWidth));

    internal void RestoreSession(ExplorerSession saved)
    {
        SessionFile.Validate(saved);
        if (!Workbench.CanRun) throw new InvalidOperationException("Stop the active calculation before opening a session.");
        foreach (var calculation in OwnedWindows.OfType<CalculationWindow>().ToArray()) calculation.Close();
        ViewModel.RestoreSession(saved);
        Workbench.RestoreHistory(saved.History);
        // Older files may contain a panned/zoomed view. Refit at the current layout
        // when the real plot becomes visible, including after opening in torus mode.
        realViewResetPending = true;
        RestoreSidebar(equationSidebar, saved.EquationPanel, EquationColumn, EquationPanel, EquationTab, EquationSplitter);
        RestoreSidebar(resultsSidebar, saved.ResultsPanel, ResultsColumn, Results, ResultsTab, ResultsSplitter);
        CoefficientsExpander.IsExpanded = saved.CoefficientsExpanded;
        Plot.RestoreView(Plot.GetResetView());
        TorusView.Fit();
        TorusView.RestoreSelection(saved.SelectedTorusPoint);
        ViewMode.SelectedIndex = saved.ComplexView ? 1 : 0;
        if (realViewResetPending && !IsComplexView) QueueRealViewReset();
        UpdateSidebarBounds();
        var version = ++sessionRestoreVersion;
        MarkSessionClean();
        var baselineVersion = cleanSessionVersion;
        Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() =>
        {
            if (version != sessionRestoreVersion) return;
            var equationOffset = EquationScroll.VerticalOffset;
            var torusOffset = TorusView.ScrollOffset;
            EquationScroll.ScrollToVerticalOffset(saved.EquationScrollOffset);
            TorusView.RestoreScroll(saved.TorusScrollOffset);
            // Account for layout clamping, without marking unrelated edits clean.
            EquationScroll.UpdateLayout();
            TorusView.UpdateLayout();
            if (cleanSession != null && baselineVersion == cleanSessionVersion)
                cleanSession = cleanSession with
                {
                    EquationScrollOffset = cleanSession.EquationScrollOffset == equationOffset
                        ? EquationScroll.VerticalOffset : cleanSession.EquationScrollOffset,
                    TorusScrollOffset = cleanSession.TorusScrollOffset == torusOffset
                        ? TorusView.ScrollOffset : cleanSession.TorusScrollOffset
                };
        }));
    }

    private static void RestoreSidebar(SidebarState state, SidebarSession saved, ColumnDefinition column,
        FrameworkElement panel, Button tab, GridSplitter splitter)
    {
        state.AnimationVersion++;
        state.IsAnimating = false;
        state.IsVisible = saved.Visible;
        state.ExpandedWidth = saved.Width;
        column.BeginAnimation(ColumnDefinition.WidthProperty, null);
        column.MinWidth = saved.Visible ? state.MinimumWidth : SidebarTabWidth;
        column.MaxWidth = double.PositiveInfinity;
        column.Width = new GridLength(saved.Visible ? saved.Width : SidebarTabWidth);
        panel.Width = double.NaN;
        panel.HorizontalAlignment = HorizontalAlignment.Stretch;
        panel.Visibility = splitter.Visibility = saved.Visible ? Visibility.Visible : Visibility.Collapsed;
        tab.Visibility = saved.Visible ? Visibility.Collapsed : Visibility.Visible;
    }
}

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
    Func<string, ExplorerSession, Task<string>>? WriteSession = null);

public partial class MainWindow
{
    private string? sessionPath;
    private int sessionRestoreVersion;
    private ExplorerSession? cleanSession;
    private SessionDialogs sessionDialogs = null!;
    private bool sessionActionInProgress, isSavingSession, sessionSaveFailed, allowSessionClose;
    internal Task<bool> PendingSessionOperation { get; private set; } = Task.FromResult(true);

    private void SessionCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
    {
        var state = CurrentSessionState;
        e.CanExecute = e.Command == ApplicationCommands.Save ? state.CanSave
            : e.Command == ApplicationCommands.SaveAs ? state.CanSaveAs : state.CanReplace;
        e.Handled = true;
    }

    private async void SessionCommandExecuted(object sender, ExecutedRoutedEventArgs e)
    {
        Session.Close();
        Explorer.Close();
        Edit.Close();
        Help.Close();
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
                var dialog = CreateSessionSaveDialog(path);
                return dialog.ShowDialog(this) == true ? dialog.FileName : null;
            },
            () =>
            {
                var dialog = new OpenFileDialog
                {
                    Filter = SessionFile.DialogFilter, DefaultExt = SessionFile.Extension,
                    Title = SessionMessages.OpenTitle, CheckFileExists = true, Multiselect = false
                };
                return dialog.ShowDialog(this) == true ? dialog.FileName : null;
            },
            (title, message) => ConfirmationWindow.ShowMessage(this, title, message));
        MarkSessionClean();
        InitializeSessionStatus();
    }

    internal static SaveFileDialog CreateSessionSaveDialog(string? path)
    {
        var dialog = new SaveFileDialog
        {
            Filter = SessionFile.DialogFilter, DefaultExt = SessionFile.Extension, AddExtension = true,
            FileName = SessionFile.GetFileName(path), Title = SessionMessages.SaveAsTitle, OverwritePrompt = true
        };
        if (path != null && Path.IsPathFullyQualified(path))
        {
            dialog.InitialDirectory = Path.GetDirectoryName(path)!;
        }
        return dialog;
    }

    internal bool HasUnsavedChanges => Workbench.IsBusy || HasSessionEdits;
    private bool HasSessionEdits => cleanSession != null &&
        (ViewModel.HasIncompleteInput
         || !SessionChanges.Equal(cleanSession, ReadSession()));

    private void MarkSessionClean()
    {
        cleanSession = ReadSession();
        QueueSessionStatusRefresh();
    }

    private async Task<bool> ConfirmSessionChangeAsync()
    {
        if (!HasUnsavedChanges) return true;
        var result = sessionDialogs.ConfirmUnsaved(SessionFile.GetFileName(sessionPath));
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
            sessionDialogs.ShowError(SessionMessages.OpenTitle, "The session could not be opened.\n\n" + error.Message);
            return false;
        }
    });

    internal Task<bool> TrySaveSessionAsync(bool saveAs = false, string? fileName = null)
    {
        var state = CurrentSessionState;
        if (!(saveAs ? state.CanSaveAs : state.CanSave)) return Task.FromResult(false);
        return RunSessionOperation(() => SaveSessionCoreAsync(saveAs, fileName));
    }

    private async Task<bool> SaveSessionCoreAsync(bool saveAs, string? fileName)
    {
        ExplorerSession saved;
        try { saved = CaptureSession(); }
        catch (InvalidOperationException error)
        {
            sessionSaveFailed = true;
            RefreshSessionStatus();
            sessionDialogs.ShowError(SessionMessages.SaveTitle, error.Message);
            return false;
        }
        var suggestedPath = SessionFile.SuggestSavePath(sessionPath, fileName);
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
            path = sessionDialogs.WriteSession is { } write
                ? await write(path, saved)
                : await Task.Run(() => SessionFile.Save(path, saved));
            sessionPath = Path.GetFullPath(path);
            cleanSession = saved;
            return true;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            sessionSaveFailed = true;
            isSavingSession = false;
            RefreshSessionStatus();
            sessionDialogs.ShowError(SessionMessages.SaveTitle, "The session could not be saved.\n\n" + error.Message);
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
        ViewModel.FlushUpdate();
        if (ViewModel.HasIncompleteInput)
            throw new InvalidOperationException("Finish the curve equation before saving the session.");
        return ReadSession().DataOnly();
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

    internal void RestoreSession(ExplorerSession saved)
    {
        saved = saved.DataOnly();
        SessionFile.Validate(saved);
        if (!Workbench.CanRun) throw new InvalidOperationException(SessionMessages.StopCalculationBeforeOpen);
        restoringHistory = true;
        historyRestoreVersion++;
        try
        {
            foreach (var dialog in OwnedWindows.Cast<Window>()
                .Where(window => window is CalculationWindow or LmfdbImportWindow).ToArray()) dialog.Close();
            // Every document opens in Real locus, fitted to the current layout.
            // Hide the torus before changing its inputs to avoid starting period work.
            ViewMode.SelectedIndex = 0;
            ViewModel.RestoreSession(saved);
            Workbench.RestoreHistory(saved.History);
            realViewResetPending = true;
            RestoreSidebar(equationSidebar, saved.EquationPanel, EquationColumn, EquationPanel, EquationTab, EquationSplitter);
            RestoreSidebar(resultsSidebar, saved.ResultsPanel, ResultsColumn, Results, ResultsTab, ResultsSplitter);
            CoefficientsExpander.IsExpanded = saved.CoefficientsExpanded;
            Plot.RestoreView(Plot.GetResetView());
            TorusView.Fit();
            TorusView.RestoreSelection(saved.SelectedTorusPoint);
            QueueRealViewReset();
            UpdateSidebarBounds();
            var version = ++sessionRestoreVersion;
            MarkSessionClean();
            Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() =>
            {
                if (version != sessionRestoreVersion) return;
                EquationScroll.ScrollToVerticalOffset(saved.EquationScrollOffset);
                TorusView.RestoreScroll(saved.TorusScrollOffset);
                ResetInitialHistory();
            }));
        }
        finally { restoringHistory = false; }
        ResetHistory();
    }

}

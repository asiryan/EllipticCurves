using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

internal sealed record SessionDialogs(Func<string, SaveChangesResult> ConfirmUnsaved,
    Func<string?, string?> ChooseSavePath, Func<string?> ChooseOpenPath, Action<string, string> ShowError);

public partial class MainWindow
{
    private string? sessionPath;
    private int sessionRestoreVersion;
    private int cleanSessionVersion;
    private bool initialSessionRendered;
    private ExplorerSession? cleanSession;
    private SessionDialogs sessionDialogs = null!;

    private void SessionCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
    {
        e.CanExecute = e.Command == ApplicationCommands.Save || Workbench.CanRun;
        e.Handled = true;
    }

    private void SessionCommandExecuted(object sender, ExecutedRoutedEventArgs e)
    {
        Session.Close();
        Explorer.Close();
        if (e.Command == ApplicationCommands.New) NewSession();
        else if (e.Command == ApplicationCommands.Open) OpenSession();
        else if (e.Command == ApplicationCommands.Save) TrySaveSession();
        e.Handled = true;
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
                    FileName = path ?? "session.ec", Title = "Save session", OverwritePrompt = true
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
    }

    internal bool HasUnsavedChanges => cleanSession != null &&
        (Workbench.IsBusy || ViewModel.HasIncompleteInput || !ViewModel.Step.IsValid
         || !SessionChanges.Equal(cleanSession, ReadSession()));

    private void MarkSessionClean()
    {
        cleanSession = ReadSession();
        cleanSessionVersion++;
    }

    private bool ConfirmSessionChange()
    {
        if (!HasUnsavedChanges) return true;
        var result = sessionDialogs.ConfirmUnsaved(sessionPath == null ? "session.ec" : Path.GetFileName(sessionPath));
        return result.Choice switch
        {
            SaveChangesChoice.Save => TrySaveSession(result.FileName),
            SaveChangesChoice.Discard => true,
            _ => false
        };
    }

    internal bool NewSession()
    {
        if (!Workbench.CanRun || !ConfirmSessionChange()) return false;
        RestoreSession(ExplorerSession.New());
        sessionPath = null;
        return true;
    }

    internal bool OpenSession()
    {
        if (!Workbench.CanRun || !ConfirmSessionChange()) return false;
        var path = sessionDialogs.ChooseOpenPath();
        if (path == null) return false;
        try
        {
            var saved = SessionFile.Load(path);
            RestoreSession(saved);
            sessionPath = Path.GetFullPath(path);
            return true;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            sessionDialogs.ShowError("Open session", "The session could not be opened.\n\n" + error.Message);
            return false;
        }
    }

    internal bool TrySaveSession(string? fileName = null)
    {
        ExplorerSession saved;
        try { saved = CaptureSession(); }
        catch (InvalidOperationException error)
        {
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
        var path = sessionDialogs.ChooseSavePath(suggestedPath);
        if (path == null) return false;
        try
        {
            // Capture again after the file dialog: an active calculation may have
            // completed while it was open.
            saved = CaptureSession();
            SessionFile.Save(path, saved);
            sessionPath = Path.GetFullPath(path);
            cleanSession = saved;
            cleanSessionVersion++;
            return true;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            sessionDialogs.ShowError("Save session", "The session could not be saved.\n\n" + error.Message);
            return false;
        }
    }

    internal ExplorerSession CaptureSession()
    {
        ViewModel.Equation.CommitEdit();
        ViewModel.Step.CommitEdit();
        ViewModel.FlushUpdate();
        if (ViewModel.HasIncompleteInput || !ViewModel.Step.IsValid)
            throw new InvalidOperationException("Finish the curve equation and enter a valid slider step before saving the session.");
        return ReadSession();
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
            History = Workbench.Jobs.Select(job => job.CaptureSession()).ToList(),
            SelectedResult = Workbench.Selected == null ? -1 : Workbench.Jobs.IndexOf(Workbench.Selected)
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
        Workbench.RestoreHistory(saved.History, saved.SelectedResult);
        // Supersede pending fits and animations; they must not overwrite the saved view.
        realViewResetPending = saved.FitRealViewWhenShown;
        RestoreSidebar(equationSidebar, saved.EquationPanel, EquationColumn, EquationPanel, EquationTab, EquationSplitter);
        RestoreSidebar(resultsSidebar, saved.ResultsPanel, ResultsColumn, Results, ResultsTab, ResultsSplitter);
        CoefficientsExpander.IsExpanded = saved.CoefficientsExpanded;
        Plot.RestoreView(saved.Plot);
        TorusView.RestoreCamera(saved.TorusCamera);
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

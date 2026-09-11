using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private string? sessionPath;
    private int sessionRestoreVersion;

    private void OpenSession()
    {
        if (!Workbench.CanRun) return;
        var dialog = new OpenFileDialog
        {
            Filter = "Explorer session (*.ec)|*.ec", DefaultExt = ".ec",
            Title = "Open session", CheckFileExists = true, Multiselect = false
        };
        if (dialog.ShowDialog(this) != true) return;
        try
        {
            var saved = SessionFile.Load(dialog.FileName);
            RestoreSession(saved);
            sessionPath = dialog.FileName;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            ConfirmationWindow.ShowMessage(this, "Open session", "The session could not be opened.\n\n" + error.Message);
        }
    }

    private void SaveSession()
    {
        ExplorerSession saved;
        try { saved = CaptureSession(); }
        catch (InvalidOperationException error)
        {
            ConfirmationWindow.ShowMessage(this, "Save session", error.Message);
            return;
        }
        var dialog = new SaveFileDialog
        {
            Filter = "Explorer session (*.ec)|*.ec", DefaultExt = ".ec", AddExtension = true,
            FileName = sessionPath ?? "session.ec", Title = "Save session", OverwritePrompt = true
        };
        if (dialog.ShowDialog(this) != true) return;
        try
        {
            // Capture again after the file dialog: an active calculation may have
            // completed while it was open.
            saved = CaptureSession();
            SessionFile.Save(dialog.FileName, saved);
            sessionPath = dialog.FileName;
        }
        catch (Exception error) when (error is IOException or InvalidDataException or UnauthorizedAccessException or InvalidOperationException)
        {
            ConfirmationWindow.ShowMessage(this, "Save session", "The session could not be saved.\n\n" + error.Message);
        }
    }

    internal ExplorerSession CaptureSession()
    {
        ViewModel.Equation.CommitEdit();
        ViewModel.Step.CommitEdit();
        ViewModel.FlushUpdate();
        if (ViewModel.HasIncompleteInput || !ViewModel.Step.IsValid)
            throw new InvalidOperationException("Finish the curve equation and enter a valid slider step before saving the session.");
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
        Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() =>
        {
            if (version != sessionRestoreVersion) return;
            EquationScroll.ScrollToVerticalOffset(saved.EquationScrollOffset);
            TorusView.RestoreScroll(saved.TorusScrollOffset);
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

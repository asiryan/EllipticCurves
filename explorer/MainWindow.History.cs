using System.ComponentModel;
using System.Windows;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    public EditHistoryViewModel EditHistory { get; } = new();
    private readonly MementoHistory<WorkspaceMemento> editHistory = new();
    private readonly DispatcherTimer editTimer = new() { Interval = TimeSpan.FromMilliseconds(450) };
    private WorkspaceMemento? editCheckpoint;
    private bool restoringHistory, committingHistory;
    private int torusSelectionRevision, historyRestoreVersion;
    private readonly List<INotifyPropertyChanged> editObservers = new();

    private sealed record WorkspaceMemento(MainViewModel.Memento Editor, WorkbenchViewModel.Memento Results,
        PlotViewState Plot, TorusCameraState Camera, bool Complex, bool FitPending, ComplexTorusViewModel.Memento Torus, int SelectionRevision)
    {
        public bool SameEdit(WorkspaceMemento other) => Editor.SameEdit(other.Editor) && Results.SameEdit(other.Results)
            && Plot == other.Plot && Camera == other.Camera && Complex == other.Complex
            && SelectionRevision == other.SelectionRevision;
    }

    private void InitializeHistory()
    {
        Edit.DataContext = EditHistory;
        Edit.UndoRequested += () => ApplicationCommands.Undo.Execute(null, this);
        Edit.RedoRequested += () => ApplicationCommands.Redo.Execute(null, this);
        editTimer.Tick += HistoryTimerTick;
        foreach (var model in new INotifyPropertyChanged[] { ViewModel, ViewModel.Equation, ViewModel.Step }
            .Concat(ViewModel.Coefficients).Concat(ViewModel.SimpleCoefficients))
        {
            editObservers.Add(model);
            model.PropertyChanged += EditorHistoryChanged;
        }
        Workbench.HistoryChanging += CommitHistory;
        Workbench.HistoryChanged += CommitHistory;
        Workbench.HistoryReplaced += ResetHistory;
        Workbench.PropertyChanged += WorkbenchHistoryStateChanged;
        Plot.ViewChanged += QueueHistory;
        TorusView.ViewChanged += QueueHistory;
        TorusView.Model.SelectionEdited += TorusSelectionEdited;
        ResetHistory();
    }

    private WorkspaceMemento CaptureHistory() => new(ViewModel.CaptureMemento(), Workbench.CaptureMemento(),
        Plot.CaptureView(), TorusView.CaptureCamera(), IsComplexView, realViewResetPending, TorusView.CaptureMemento(), torusSelectionRevision);

    private void TorusSelectionEdited() { if (!restoringHistory) torusSelectionRevision++; QueueHistory(); }

    private void ResetInitialHistory()
    {
        if (editCheckpoint != null && !editHistory.CanUndo && !editHistory.CanRedo && ViewModel.CaptureMemento().SameEdit(editCheckpoint.Editor)
            && Workbench.CaptureMemento().SameEdit(editCheckpoint.Results)) ResetHistory();
    }

    private void EditorHistoryChanged(object? sender, PropertyChangedEventArgs e)
    {
        // Worker progress and rendering results enrich existing mementos; they
        // are not user edits and must not add steps or invalidate Redo.
        if (e.PropertyName is "Text" or "SliderOffset" or "SelectedPreset" or "ShowGrid" or "ShowPoints") QueueHistory();
    }

    private void WorkbenchHistoryStateChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName == nameof(WorkbenchViewModel.CanRun)) RefreshHistoryCommands();
    }

    private void QueueHistory()
    {
        if (editCheckpoint == null || restoringHistory || committingHistory || sessionClosed) return;
        editTimer.Stop();
        editTimer.Start();
        RefreshHistoryCommands();
    }

    private void HistoryTimerTick(object? sender, EventArgs e)
    {
        // A pause while dragging is still part of the same gesture.
        if (Mouse.Captured != null) return;
        CommitHistory();
    }

    internal void CommitHistory()
    {
        if (editCheckpoint == null || restoringHistory || committingHistory || sessionClosed) return;
        editTimer.Stop();
        committingHistory = true;
        try
        {
            ViewModel.FlushUpdate();
            var state = CaptureHistory();
            if (!state.SameEdit(editCheckpoint)) editHistory.Record(editCheckpoint);
            editCheckpoint = state;
        }
        finally { committingHistory = false; }
        RefreshHistoryCommands();
    }

    internal void ResetHistory()
    {
        if (restoringHistory || sessionClosed) return;
        editTimer.Stop();
        editHistory.Clear();
        committingHistory = true;
        try { editCheckpoint = CaptureHistory(); }
        finally { committingHistory = false; }
        RefreshHistoryCommands();
    }

    private void RefreshHistoryCommands()
    {
        if (editCheckpoint == null || restoringHistory || committingHistory || sessionClosed) return;
        // Comparing editor values and camera positions is cheap. Do not capture
        // reports or touch the torus worker while WPF queries a command.
        var pending = !ViewModel.CaptureMemento().SameEdit(editCheckpoint.Editor)
            || Plot.CaptureView() != editCheckpoint.Plot || TorusView.CaptureCamera() != editCheckpoint.Camera
            || IsComplexView != editCheckpoint.Complex || torusSelectionRevision != editCheckpoint.SelectionRevision;
        var available = Workbench.CanRun && !sessionActionInProgress;
        var canUndo = available && (editHistory.CanUndo || pending);
        var canRedo = available && editHistory.CanRedo && !pending;
        if (canUndo == EditHistory.CanUndo && canRedo == EditHistory.CanRedo) return;
        EditHistory.Update(canUndo, canRedo);
        CommandManager.InvalidateRequerySuggested();
    }

    private void EditCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
    {
        RefreshHistoryCommands();
        e.CanExecute = e.Command == ApplicationCommands.Undo ? EditHistory.CanUndo : EditHistory.CanRedo;
        e.Handled = true;
    }

    private void EditCommandExecuted(object sender, ExecutedRoutedEventArgs e)
    {
        e.Handled = true;
        Session.Close();
        Explorer.Close();
        Edit.Close();
        if (!Workbench.CanRun || sessionActionInProgress) return;
        CommitHistory();
        var undo = e.Command == ApplicationCommands.Undo;
        if (undo ? !editHistory.CanUndo : !editHistory.CanRedo) return;
        var current = CaptureHistory();
        var state = undo ? editHistory.Undo(current) : editHistory.Redo(current);
        restoringHistory = true;
        historyRestoreVersion++;
        try
        {
            // Restoration bypasses normal setters that schedule computation.
            // Results/cameras are restored directly, including incomplete input.
            ViewModel.RestoreMemento(state.Editor);
            Workbench.RestoreMemento(state.Results);
            realViewResetPending = false;
            ViewMode.SelectedIndex = state.Complex ? 1 : 0;
            Plot.RestoreView(state.Plot);
            TorusView.RestoreCamera(state.Camera);
            TorusView.RestoreMemento(state.Torus, ViewModel.Snapshot, ViewModel.Samples, state.Complex);
            torusSelectionRevision = state.SelectionRevision;
            realViewResetPending = state.FitPending;
            editCheckpoint = state;
        }
        finally { restoringHistory = false; }
        RefreshHistoryCommands();
        RefreshSessionStatus();
    }

    private void DisposeHistory()
    {
        editTimer.Stop();
        editTimer.Tick -= HistoryTimerTick;
        foreach (var model in editObservers) model.PropertyChanged -= EditorHistoryChanged;
        editObservers.Clear();
        Workbench.HistoryChanging -= CommitHistory;
        Workbench.HistoryChanged -= CommitHistory;
        Workbench.HistoryReplaced -= ResetHistory;
        Workbench.PropertyChanged -= WorkbenchHistoryStateChanged;
        Plot.ViewChanged -= QueueHistory;
        TorusView.ViewChanged -= QueueHistory;
        TorusView.Model.SelectionEdited -= TorusSelectionEdited;
        editHistory.Clear();
        editCheckpoint = null;
    }
}

using System.Collections.Specialized;
using System.ComponentModel;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.ViewModels;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    public SessionStatusViewModel SessionStatus { get; } = new();
    private readonly HashSet<INotifyPropertyChanged> sessionObservers = new();
    private readonly HashSet<CalculationJobViewModel> observedJobs = new();
    private DispatcherOperation? sessionStatusRefresh;
    private bool sessionClosed;
    private SessionState CurrentSessionState => new(sessionPath, HasUnsavedChanges, isSavingSession,
        sessionSaveFailed, sessionActionInProgress, Workbench.CanRun);

    private void InitializeSessionStatus()
    {
        foreach (var model in new INotifyPropertyChanged[] { ViewModel, ViewModel.Equation, ViewModel.Step, Workbench, TorusView.Model }
            .Concat(ViewModel.Coefficients).Concat(ViewModel.SimpleCoefficients))
            WatchSessionModel(model);
        Workbench.Jobs.CollectionChanged += SessionHistoryChanged;
        RefreshSessionJobs();
        EquationScroll.ScrollChanged += SessionScrollChanged;
        TorusView.AddHandler(ScrollViewer.ScrollChangedEvent, new ScrollChangedEventHandler(SessionScrollChanged));
        CoefficientsExpander.Expanded += SessionLayoutChanged;
        CoefficientsExpander.Collapsed += SessionLayoutChanged;
        RefreshSessionStatus();
    }

    private void WatchSessionModel(INotifyPropertyChanged model)
    {
        if (sessionObservers.Add(model)) model.PropertyChanged += SessionModelChanged;
    }

    private void RefreshSessionJobs()
    {
        foreach (var job in observedJobs.Where(job => !Workbench.Jobs.Contains(job)).ToArray())
        {
            job.PropertyChanged -= SessionModelChanged;
            observedJobs.Remove(job);
            sessionObservers.Remove(job);
        }
        foreach (var job in Workbench.Jobs)
            if (observedJobs.Add(job)) WatchSessionModel(job);
    }

    private void SessionHistoryChanged(object? sender, NotifyCollectionChangedEventArgs e)
    {
        RefreshSessionJobs();
        QueueSessionStatusRefresh();
    }

    private void SessionModelChanged(object? sender, PropertyChangedEventArgs e) => QueueSessionStatusRefresh();
    private void SessionScrollChanged(object sender, ScrollChangedEventArgs e) => QueueSessionStatusRefresh();
    private void SessionLayoutChanged(object sender, RoutedEventArgs e) => QueueSessionStatusRefresh();

    private void QueueSessionStatusRefresh()
    {
        if (cleanSession == null || sessionClosed || sessionStatusRefresh?.Status == DispatcherOperationStatus.Pending) return;
        sessionStatusRefresh = Dispatcher.BeginInvoke(DispatcherPriority.Background, new Action(RefreshSessionStatus));
    }

    internal void RefreshSessionStatus()
    {
        if (sessionClosed) return;
        SessionStatus.Update(CurrentSessionState);
        RefreshHistoryCommands();
        var title = SessionStatus.DisplayName + " — " + ExplorerInfo.WindowTitle;
        if (Title != title) Title = title;
        CommandManager.InvalidateRequerySuggested();
    }

    private void DisposeSessionStatus()
    {
        sessionClosed = true;
        sessionStatusRefresh?.Abort();
        foreach (var model in sessionObservers) model.PropertyChanged -= SessionModelChanged;
        sessionObservers.Clear();
        observedJobs.Clear();
        Workbench.Jobs.CollectionChanged -= SessionHistoryChanged;
        EquationScroll.ScrollChanged -= SessionScrollChanged;
        TorusView.RemoveHandler(ScrollViewer.ScrollChangedEvent, new ScrollChangedEventHandler(SessionScrollChanged));
        CoefficientsExpander.Expanded -= SessionLayoutChanged;
        CoefficientsExpander.Collapsed -= SessionLayoutChanged;
    }
}

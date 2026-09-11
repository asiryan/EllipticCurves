using System.ComponentModel;
using System.Diagnostics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.ViewModels;
using EllipticCurves.Explorer.Windowing;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer;

public partial class MainWindow : Window
{
    public MainViewModel ViewModel { get; } = new();
    public WorkbenchViewModel Workbench { get; } = new();
    private readonly Func<bool> confirmEquationReset;

    public MainWindow() : this(null) { }
    internal MainWindow(Func<bool>? confirmation, SessionDialogs? sessionDialogs = null)
    {
        confirmEquationReset = confirmation ?? (() => ConfirmationWindow.Confirm(this,
            "Reset equation?", "Restore the classic curve and recenter the plot. Your current equation will be replaced.",
            "Reset equation", CurvePreset.ClassicEquation));
        InitializeComponent();
        DataContext = ViewModel;
        Results.DataContext = Workbench;
        ResultsTab.DataContext = Workbench;
        Session.DataContext = SessionStatus;
        Session.NewRequested += () => ApplicationCommands.New.Execute(null, this);
        Session.OpenRequested += () => ApplicationCommands.Open.Execute(null, this);
        Session.SaveRequested += () => ApplicationCommands.Save.Execute(null, this);
        Session.SaveAsRequested += () => ApplicationCommands.SaveAs.Execute(null, this);
        Session.ExitRequested += Close;
        Explorer.OperationRequested += OpenCalculation;
        Results.HideRequested += () => SetResultsVisible(false);
        Results.RepeatRequested += request => OpenCalculation(CalculationCatalog.Get(request.OperationId), request);
        ViewModel.ViewResetRequested += ResetView;
        ViewModel.CurveResetRequested += ResetCurveViews;
        TorusView.Model.PropertyChanged += TorusStateChanged;
        SourceInitialized += UpdateWindowInsets;
        StateChanged += UpdateWindowInsets;
        LocationChanged += UpdateWindowInsets;
        SizeChanged += UpdateWindowInsets;
        Workspace.SizeChanged += (_, _) => UpdateSidebarBounds();
        EquationHost.SizeChanged += (_, _) => UpdateSidebarBounds();
        ResultsHost.SizeChanged += (_, _) => UpdateSidebarBounds();
        DpiChanged += (_, _) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded,
            new Action(() => AppRoot.Margin = WindowWorkArea.GetContentMargin(this)));
        InitializeSession(sessionDialogs);
    }

    private void UpdateWindowInsets(object? sender, EventArgs e)
    {
        AppRoot.Margin = WindowWorkArea.GetContentMargin(this);
        UpdateSidebarBounds();
    }
    private void WindowLoaded(object sender, RoutedEventArgs e)
    {
        if (sessionRestoreVersion == 0) ResetCurveViews(sender, e);
    }

    protected override void OnClosing(CancelEventArgs e)
    {
        base.OnClosing(e);
        if (e.Cancel || allowSessionClose) return;
        if (sessionActionInProgress) { e.Cancel = true; return; }
        if (!HasUnsavedChanges) return;
        var confirmation = RunSessionOperation(ConfirmSessionChangeAsync);
        if (confirmation.IsCompleted) e.Cancel = !confirmation.GetAwaiter().GetResult();
        else
        {
            e.Cancel = true;
            PendingSessionOperation = FinishSessionCloseAsync(confirmation);
        }
    }

    private void WindowClosed(object? sender, EventArgs e)
    {
        DisposeSessionStatus();
        Workbench.Dispose();
        TorusView.Model.PropertyChanged -= TorusStateChanged;
        TorusView.Dispose();
        ViewModel.ViewResetRequested -= ResetView;
        ViewModel.CurveResetRequested -= ResetCurveViews;
        ViewModel.Dispose();
    }

    private void RepositoryClick(object sender, RoutedEventArgs e)
    {
        try
        {
            using var browser = Process.Start(new ProcessStartInfo(ExplorerInfo.RepositoryUrl) { UseShellExecute = true });
        }
        catch (Exception error) when (error is Win32Exception or InvalidOperationException)
        {
            ConfirmationWindow.ShowMessage(this, "Open GitHub repository",
                "Could not open the browser. Open this address manually:\n" + ExplorerInfo.RepositoryUrl);
        }
    }
    private void OpenCalculation(CalculationOperation operation) => OpenCalculation(operation, null);
    private void OpenCalculation(CalculationOperation operation, CalculationRequest? previous)
    {
        if (previous == null && operation.UsesPlot)
        {
            ViewModel.Equation.CommitEdit();
            ViewModel.FlushUpdate();
            if (ViewModel.HasInputError || ViewModel.HasIncompleteInput)
            {
                ConfirmationWindow.ShowMessage(this, "Check the equation", "Finish the curve equation before starting a calculation.");
                return;
            }
        }
        var dialog = new CalculationWindow(operation,
            previous?.Equation ?? CurveEquationText.Format(ViewModel.Snapshot.Curve), Workbench, previous) { Owner = this };
        dialog.RunRequested += RunCalculation;
        dialog.Show();
    }
    private async void RunCalculation(CalculationRequest request)
    {
        SetResultsVisible(true);
        await Workbench.RunAsync(request);
    }
    private void ResetEquationClick(object sender, RoutedEventArgs e)
    {
        if (confirmEquationReset()) ViewModel.ResetCommand.Execute(null);
    }
    private void MinimizeClick(object sender, RoutedEventArgs e) => SystemCommands.MinimizeWindow(this);
    private void MaximizeClick(object sender, RoutedEventArgs e)
    {
        if (WindowState == WindowState.Maximized) SystemCommands.RestoreWindow(this);
        else SystemCommands.MaximizeWindow(this);
    }

    private void CloseClick(object sender, RoutedEventArgs e) => Close();

    private void CoefficientEditFinished(object sender, KeyboardFocusChangedEventArgs e)
    {
        if (sender is FrameworkElement { DataContext: CoefficientViewModel coefficient }) coefficient.CommitEdit();
    }

    private void EquationEditFinished(object sender, KeyboardFocusChangedEventArgs e) => ViewModel.Equation.CommitEdit();

    private void EquationKeyDown(object sender, KeyEventArgs e)
    {
        if (e.Key != Key.Enter) return;
        ViewModel.Equation.CommitEdit();
        ViewModel.FlushUpdate();
        e.Handled = true;
    }

    private void CoefficientKeyDown(object sender, KeyEventArgs e)
    {
        if (e.Key != Key.Enter || sender is not FrameworkElement { DataContext: CoefficientViewModel coefficient }) return;
        coefficient.CommitEdit();
        ViewModel.FlushUpdate();
        e.Handled = true;
    }

    private void ResetSlider(object sender, MouseButtonEventArgs e)
    {
        if (sender is not FrameworkElement { DataContext: CoefficientViewModel coefficient }) return;
        coefficient.ResetCommand.Execute(null);
        ViewModel.FlushUpdate();
        e.Handled = true;
    }

    private async void CopyClick(object sender, RoutedEventArgs e)
    {
        if (ClipboardActions.CopyText(this, ViewModel.Snapshot.Summary, "Copy") && sender is Button button)
        {
            button.Content = "Copied ✓";
            await Task.Delay(1400);
            button.Content = "Copy";
        }
    }

}

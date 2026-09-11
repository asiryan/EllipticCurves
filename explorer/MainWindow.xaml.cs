using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Animation;
using System.Windows.Threading;
using EllipticCurves.Explorer.ViewModels;
using EllipticCurves.Explorer.Windowing;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

public partial class MainWindow : Window
{
    public MainViewModel ViewModel { get; } = new();
    public WorkbenchViewModel Workbench { get; } = new();
    private const double SidebarTabWidth = 32;
    private readonly SidebarState equationSidebar = new(238);
    private readonly SidebarState resultsSidebar = new(300);
    private readonly Func<bool> confirmEquationReset;

    private sealed class SidebarState(double minimumWidth)
    {
        public double MinimumWidth { get; } = minimumWidth;
        public double ExpandedWidth { get; set; } = minimumWidth;
        public bool IsVisible { get; set; } = true;
        public bool IsAnimating { get; set; }
        public int AnimationVersion { get; set; }
    }

    public MainWindow() : this(null) { }
    internal MainWindow(Func<bool>? confirmation)
    {
        confirmEquationReset = confirmation ?? (() => ConfirmationWindow.Confirm(this,
            "Reset equation?", "Restore the classic curve and recenter the plot. Your current equation will be replaced.",
            "Reset equation", "y^2 = x^3 - x"));
        InitializeComponent();
        DataContext = ViewModel;
        Results.DataContext = Workbench;
        ResultsTab.DataContext = Workbench;
        Explorer.OperationRequested += OpenCalculation;
        Results.HideRequested += () => SetResultsVisible(false);
        Results.RepeatRequested += request => OpenCalculation(CalculationCatalog.Get(request.OperationId), request);
        ViewModel.ViewResetRequested += ResetView;
        SourceInitialized += UpdateWindowInsets;
        StateChanged += UpdateWindowInsets;
        LocationChanged += UpdateWindowInsets;
        SizeChanged += UpdateWindowInsets;
        Workspace.SizeChanged += (_, _) => UpdateSidebarBounds();
        EquationHost.SizeChanged += (_, _) => UpdateSidebarBounds();
        ResultsHost.SizeChanged += (_, _) => UpdateSidebarBounds();
        DpiChanged += (_, _) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded,
            new Action(() => AppRoot.Margin = WindowWorkArea.GetContentMargin(this)));
    }

    private void UpdateWindowInsets(object? sender, EventArgs e)
    {
        AppRoot.Margin = WindowWorkArea.GetContentMargin(this);
        UpdateSidebarBounds();
    }
    private void UpdateSidebarBounds()
    {
        var available = Workspace.ActualWidth - PlotColumn.MinWidth
            - Workspace.ColumnDefinitions[1].Width.Value - Workspace.ColumnDefinitions[3].Width.Value;
        // Resolve both limits from requested widths so resizing the window cannot
        // make the panels repeatedly push each other's measured width back and forth.
        var equationWidth = Math.Min(EquationColumn.Width.Value,
            Math.Max(EquationColumn.MinWidth, available - ResultsColumn.MinWidth));
        ResultsColumn.MaxWidth = Math.Max(resultsSidebar.MinimumWidth, available - equationWidth);
        var resultsWidth = Math.Min(ResultsColumn.Width.Value, ResultsColumn.MaxWidth);
        EquationColumn.MaxWidth = Math.Max(equationSidebar.MinimumWidth, available - resultsWidth);
    }

    private void WindowLoaded(object sender, RoutedEventArgs e) => Plot.Fit();

    private void WindowClosed(object? sender, EventArgs e)
    {
        Workbench.Dispose();
        ViewModel.ViewResetRequested -= ResetView;
        ViewModel.Dispose();
    }

    private void RepositoryClick(object sender, RoutedEventArgs e)
    {
        const string repositoryUrl = "https://github.com/asiryan/EllipticCurves";
        try
        {
            using var browser = Process.Start(new ProcessStartInfo(repositoryUrl) { UseShellExecute = true });
        }
        catch (Exception error) when (error is Win32Exception or InvalidOperationException)
        {
            ConfirmationWindow.ShowMessage(this, "Open GitHub repository",
                "Could not open the browser. Open this address manually:\n" + repositoryUrl);
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
    private void ExpandResultsClick(object sender, RoutedEventArgs e) => SetResultsVisible(true);
    private void CollapseEquationClick(object sender, RoutedEventArgs e) => SetEquationVisible(false);
    private void ExpandEquationClick(object sender, RoutedEventArgs e) => SetEquationVisible(true);
    private void SetEquationVisible(bool visible) =>
        SetSidebarVisible(equationSidebar, visible, EquationColumn, EquationPanel, EquationTab, HorizontalAlignment.Right, EquationSplitter);
    private void SetResultsVisible(bool visible) =>
        SetSidebarVisible(resultsSidebar, visible, ResultsColumn, Results, ResultsTab, HorizontalAlignment.Left, ResultsSplitter);

    private void SetSidebarVisible(SidebarState state, bool visible, ColumnDefinition column,
        FrameworkElement panel, Button tab, HorizontalAlignment slideAlignment, GridSplitter? splitter = null)
    {
        if (state.IsVisible == visible) return;
        var from = column.ActualWidth > 0 ? column.ActualWidth : column.Width.Value;
        if (!visible && !state.IsAnimating) state.ExpandedWidth = Math.Max(state.MinimumWidth, from);
        var to = visible ? Math.Min(state.ExpandedWidth, column.MaxWidth) : SidebarTabWidth;
        var version = ++state.AnimationVersion;
        state.IsVisible = visible;
        column.BeginAnimation(ColumnDefinition.WidthProperty, null);
        column.MinWidth = SidebarTabWidth;
        column.Width = new GridLength(to);
        // Keep content at its expanded width while clipping it toward the outer edge.
        panel.Width = Math.Max(state.MinimumWidth, visible ? to : from);
        panel.HorizontalAlignment = slideAlignment;
        panel.Visibility = Visibility.Visible;
        tab.Visibility = Visibility.Collapsed;
        if (splitter != null) splitter.Visibility = Visibility.Collapsed;

        void Finish()
        {
            if (version != state.AnimationVersion) return;
            column.BeginAnimation(ColumnDefinition.WidthProperty, null);
            column.MinWidth = visible ? state.MinimumWidth : SidebarTabWidth;
            panel.Width = double.NaN;
            panel.HorizontalAlignment = HorizontalAlignment.Stretch;
            panel.Visibility = visible ? Visibility.Visible : Visibility.Collapsed;
            if (splitter != null) splitter.Visibility = panel.Visibility;
            tab.Visibility = visible ? Visibility.Collapsed : Visibility.Visible;
            state.IsAnimating = false;
        }

        if (!IsLoaded || !SystemParameters.ClientAreaAnimation)
        {
            Finish();
            return;
        }
        state.IsAnimating = true;
        var animation = new GridLengthAnimation { From = from, To = to, Duration = TimeSpan.FromMilliseconds(200) };
        animation.Completed += (_, _) => Finish();
        column.BeginAnimation(ColumnDefinition.WidthProperty, animation, HandoffBehavior.SnapshotAndReplace);
    }
    private void ResetView(object? sender, EventArgs e) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(Plot.Fit));
    private void ResetEquationClick(object sender, RoutedEventArgs e)
    {
        if (confirmEquationReset()) ViewModel.ResetCommand.Execute(null);
    }
    private void FitClick(object sender, RoutedEventArgs e) => Plot.Fit();
    private void ZoomInClick(object sender, RoutedEventArgs e) => Plot.Zoom(1 / 1.25);
    private void ZoomOutClick(object sender, RoutedEventArgs e) => Plot.Zoom(1.25);
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
        try
        {
            Clipboard.SetText(ViewModel.Snapshot.Summary);
            if (sender is Button button)
            {
                button.Content = "Copied ✓";
                await Task.Delay(1400);
                button.Content = "Copy";
            }
        }
        catch (ExternalException)
        {
            ConfirmationWindow.ShowMessage(this, "Copy", "The clipboard is busy. Please try again.");
        }
    }

    private void ExportClick(object sender, RoutedEventArgs e)
    {
        var dialog = new SaveFileDialog
        {
            Filter = "PNG image (*.png)|*.png",
            FileName = "elliptic-curve.png",
            Title = "Export the current real locus"
        };
        if (dialog.ShowDialog(this) != true) return;
        try
        {
            var bitmap = new RenderTargetBitmap((int)Math.Ceiling(Plot.ActualWidth * 2), (int)Math.Ceiling(Plot.ActualHeight * 2), 192, 192, PixelFormats.Pbgra32);
            bitmap.Render(Plot);
            var encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(bitmap));
            using var stream = File.Create(dialog.FileName);
            encoder.Save(stream);
        }
        catch (Exception error) when (error is IOException or UnauthorizedAccessException)
        {
            ConfirmationWindow.ShowMessage(this, "Export plot", "The image could not be saved. Check the destination and try again.");
        }
    }
}

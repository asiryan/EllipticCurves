using System.IO;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Threading;
using EllipticCurves.Visualizer.ViewModels;
using EllipticCurves.Visualizer.Windowing;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.Models;
using Microsoft.Win32;

namespace EllipticCurves.Visualizer;

public partial class MainWindow : Window
{
    public MainViewModel ViewModel { get; } = new();
    public WorkbenchViewModel Workbench { get; } = new();
    private bool resultsVisible = true;
    private double resultsWidth = 360;

    public MainWindow()
    {
        InitializeComponent();
        DataContext = ViewModel;
        Results.DataContext = Workbench;
        ResultsToggle.DataContext = Workbench;
        Explorer.OperationRequested += OpenCalculation;
        Results.HideRequested += () => SetResultsVisible(false);
        Results.RepeatRequested += request => OpenCalculation(CalculationCatalog.Get(request.OperationId), request);
        ViewModel.ViewResetRequested += ResetView;
        SourceInitialized += UpdateWindowInsets;
        StateChanged += UpdateWindowInsets;
        LocationChanged += UpdateWindowInsets;
        SizeChanged += UpdateWindowInsets;
        DpiChanged += (_, _) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded,
            new Action(() => AppRoot.Margin = WindowWorkArea.GetContentMargin(this)));
    }

    private void UpdateWindowInsets(object? sender, EventArgs e)
    {
        AppRoot.Margin = WindowWorkArea.GetContentMargin(this);
        ResultsColumn.MaxWidth = Math.Max(300, ActualWidth - 58 - 12 - 650);
    }

    private void WindowLoaded(object sender, RoutedEventArgs e) => Plot.Fit();
    private void WindowClosed(object? sender, EventArgs e) { Workbench.Dispose(); ViewModel.ViewResetRequested -= ResetView; ViewModel.Dispose(); }
    private void OpenCalculation(CalculationOperation operation) => OpenCalculation(operation, null);
    private void OpenCalculation(CalculationOperation operation, CalculationRequest? previous)
    {
        if (previous == null && operation.UsesPlot)
        {
            ViewModel.Equation.CommitEdit();
            ViewModel.FlushUpdate();
            if (ViewModel.HasInputError || ViewModel.HasIncompleteInput)
            {
                MessageBox.Show(this, "Finish the curve equation before starting a calculation.", "Explorer", MessageBoxButton.OK, MessageBoxImage.Information);
                return;
            }
        }
        var dialog = new CalculationWindow(operation, previous?.Equation ?? CurveEquationText.Format(ViewModel.Snapshot.Curve), Workbench, previous) { Owner = this };
        dialog.RunRequested += RunCalculation;
        dialog.Show();
    }
    private async void RunCalculation(CalculationRequest request)
    {
        SetResultsVisible(true);
        await Workbench.RunAsync(request);
    }
    private void ToggleResultsClick(object sender, RoutedEventArgs e) => SetResultsVisible(!resultsVisible);
    private void SetResultsVisible(bool visible)
    {
        if (resultsVisible && !visible) resultsWidth = ResultsColumn.ActualWidth;
        resultsVisible = visible;
        Results.Visibility = ResultsSplitter.Visibility = visible ? Visibility.Visible : Visibility.Collapsed;
        ResultsColumn.MinWidth = visible ? 300 : 0;
        ResultsColumn.Width = new GridLength(visible ? resultsWidth : 0);
        ResultsGap.Width = new GridLength(visible ? 12 : 0);
        EquationColumn.Width = new GridLength(visible ? 230 : 270);
        PanHint.Visibility = visible ? Visibility.Collapsed : Visibility.Visible;
    }
    private void ResetView(object? sender, EventArgs e) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(Plot.Fit));
    private void FitClick(object sender, RoutedEventArgs e) => Plot.Fit();
    private void ZoomInClick(object sender, RoutedEventArgs e) => Plot.Zoom(1 / 1.25);
    private void ZoomOutClick(object sender, RoutedEventArgs e) => Plot.Zoom(1.25);
    private void MinimizeClick(object sender, RoutedEventArgs e) => SystemCommands.MinimizeWindow(this);
    private void MaximizeClick(object sender, RoutedEventArgs e) { if (WindowState == WindowState.Maximized) SystemCommands.RestoreWindow(this); else SystemCommands.MaximizeWindow(this); }
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
        catch (ExternalException) { MessageBox.Show(this, "The clipboard is busy. Please try again.", "Copy", MessageBoxButton.OK, MessageBoxImage.Information); }
    }

    private void ExportClick(object sender, RoutedEventArgs e)
    {
        var dialog = new SaveFileDialog { Filter = "PNG image (*.png)|*.png", FileName = "elliptic-curve.png", Title = "Export the current real locus" };
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
        { MessageBox.Show(this, "The image could not be saved. Check the destination and try again.", "Export plot", MessageBoxButton.OK, MessageBoxImage.Warning); }
    }
}

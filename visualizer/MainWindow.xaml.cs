using System.IO;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Threading;
using EllipticCurves.Visualizer.ViewModels;
using EllipticCurves.Visualizer.Windowing;
using Microsoft.Win32;

namespace EllipticCurves.Visualizer;

public partial class MainWindow : Window
{
    public MainViewModel ViewModel { get; } = new();

    public MainWindow()
    {
        InitializeComponent();
        DataContext = ViewModel;
        ViewModel.ViewResetRequested += ResetView;
        SourceInitialized += UpdateWindowInsets;
        StateChanged += UpdateWindowInsets;
        LocationChanged += UpdateWindowInsets;
        SizeChanged += UpdateWindowInsets;
        DpiChanged += (_, _) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded,
            new Action(() => AppRoot.Margin = WindowWorkArea.GetContentMargin(this)));
    }

    private void UpdateWindowInsets(object? sender, EventArgs e) => AppRoot.Margin = WindowWorkArea.GetContentMargin(this);

    private void WindowLoaded(object sender, RoutedEventArgs e) => Plot.Fit();
    private void WindowClosed(object? sender, EventArgs e) { ViewModel.ViewResetRequested -= ResetView; ViewModel.Dispose(); }
    private void ResetView(object? sender, EventArgs e) => Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(Plot.Fit));
    private void FitClick(object sender, RoutedEventArgs e) => Plot.Fit();
    private void ZoomInClick(object sender, RoutedEventArgs e) => Plot.Zoom(1 / 1.25);
    private void ZoomOutClick(object sender, RoutedEventArgs e) => Plot.Zoom(1.25);
    private void MinimizeClick(object sender, RoutedEventArgs e) => SystemCommands.MinimizeWindow(this);
    private void MaximizeClick(object sender, RoutedEventArgs e) { if (WindowState == WindowState.Maximized) SystemCommands.RestoreWindow(this); else SystemCommands.MaximizeWindow(this); }
    private void CloseClick(object sender, RoutedEventArgs e) => Close();

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

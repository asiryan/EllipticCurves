using System.ComponentModel;
using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Imaging;
using System.Windows.Threading;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private bool realViewResetPending;

    private bool IsComplexView => ViewMode.SelectedIndex == 1;

    private void ViewModeChanged(object sender, SelectionChangedEventArgs e)
    {
        if (e.Source != sender || TorusView == null || FitViewButton == null) return;
        RealPlotHost.Visibility = IsComplexView ? Visibility.Collapsed : Visibility.Visible;
        TorusView.Visibility = IsComplexView ? Visibility.Visible : Visibility.Collapsed;
        PlotLegend.Visibility = RealPlotHost.Visibility;
        TorusLegend.Visibility = TorusView.Visibility;
        FitViewButton.ToolTip = IsComplexView ? "Reset torus camera (Ctrl+F)" : "Reset real plot view (Ctrl+F)";
        if (!IsComplexView && realViewResetPending) QueueRealViewReset();
        UpdateExportState();
        QueueSessionStatusRefresh();
        QueueHistory();
    }

    private void TorusStateChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName == nameof(ComplexTorusViewModel.HasLattice)) UpdateExportState();
    }

    private void UpdateExportState() => ExportPlotButton.IsEnabled = !IsComplexView || TorusView.Model.HasLattice;

    private void ResetCurveViews(object? sender, EventArgs e)
    {
        TorusView.Fit();
        realViewResetPending = true;
        QueueRealViewReset();
    }

    private void FitChangedCurve(CurveSnapshot previous, CurveSnapshot current)
    {
        if (!Plot.NeedsRefit(previous.Plot, current.Plot)) return;
        realViewResetPending = true;
        if (IsComplexView || Plot.ActualWidth <= 0 || Plot.ActualHeight <= 0) QueueRealViewReset();
        else
        {
            // Complete the visible fit before Enter/blur captures an undo checkpoint.
            // Memento restoration does not recalculate, so Undo keeps its saved camera.
            Plot.Fit(current.Plot);
            realViewResetPending = false;
        }
    }

    private void QueueRealViewReset()
    {
        var version = historyRestoreVersion;
        Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() =>
        {
            // A collapsed plot has no current layout. Keep the request until Real locus is shown.
            if (version != historyRestoreVersion || !realViewResetPending || IsComplexView) return;
            Plot.UpdateLayout();
            if (Plot.ActualWidth <= 0 || Plot.ActualHeight <= 0) return;
            Plot.Fit();
            realViewResetPending = false;
        }));
    }

    private void ResetView(object? sender, EventArgs e) => FitCurrentView();
    private void FitCurrentView()
    {
        if (IsComplexView) TorusView.Fit();
        else
        {
            realViewResetPending = true;
            QueueRealViewReset();
        }
    }
    private void FitClick(object sender, RoutedEventArgs e) => FitCurrentView();
    private void ZoomInClick(object sender, RoutedEventArgs e) => ZoomCurrentView(1 / 1.25);
    private void ZoomOutClick(object sender, RoutedEventArgs e) => ZoomCurrentView(1.25);
    private void ZoomCurrentView(double factor)
    {
        if (IsComplexView) TorusView.Zoom(factor);
        else Plot.Zoom(factor);
    }
    private void ExportClick(object sender, RoutedEventArgs e)
    {
        FrameworkElement target = IsComplexView ? TorusView : Plot;
        if (target.ActualWidth <= 0 || target.ActualHeight <= 0 || (IsComplexView && !TorusView.Model.HasLattice)) return;
        var dialog = new SaveFileDialog
        {
            Filter = "PNG image (*.png)|*.png",
            FileName = IsComplexView ? "elliptic-curve-torus.png" : "elliptic-curve.png",
            Title = IsComplexView ? "Export the period lattice and complex torus" : "Export the current real locus"
        };
        if (dialog.ShowDialog(this) != true) return;
        try
        {
            var bitmap = PlotImageExporter.Render(target);
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

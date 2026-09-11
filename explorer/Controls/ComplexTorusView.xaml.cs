using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer.Controls;

public partial class ComplexTorusView : UserControl, IDisposable
{
    public static readonly DependencyProperty SnapshotProperty = DependencyProperty.Register(nameof(Snapshot), typeof(CurveSnapshot), typeof(ComplexTorusView), new PropertyMetadata(null, InputChanged));
    public static readonly DependencyProperty SamplesProperty = DependencyProperty.Register(nameof(Samples), typeof(IReadOnlyList<EllipticCurvePoint>), typeof(ComplexTorusView), new PropertyMetadata(null, InputChanged));
    public static readonly DependencyProperty ShowGridProperty = DependencyProperty.Register(nameof(ShowGrid), typeof(bool), typeof(ComplexTorusView), new PropertyMetadata(true, GridChanged));
    private bool updateQueued, disposed;
    public ComplexTorusViewModel Model { get; } = new();
    public CurveSnapshot? Snapshot { get => (CurveSnapshot?)GetValue(SnapshotProperty); set => SetValue(SnapshotProperty, value); }
    public IReadOnlyList<EllipticCurvePoint>? Samples { get => (IReadOnlyList<EllipticCurvePoint>?)GetValue(SamplesProperty); set => SetValue(SamplesProperty, value); }
    public bool ShowGrid { get => (bool)GetValue(ShowGridProperty); set => SetValue(ShowGridProperty, value); }

    public ComplexTorusView()
    {
        InitializeComponent();
        ContentRoot.DataContext = Model;
        SizeChanged += (_, _) => ContentRoot.Height = Math.Max(380, ActualHeight - 10);
        Loaded += (_, _) => QueueUpdate();
        Unloaded += (_, _) => UpdateModel(false);
        IsVisibleChanged += (_, _) =>
        {
            if (!IsVisible) UpdateModel(false);
            else QueueUpdate();
        };
    }

    private static void InputChanged(DependencyObject sender, DependencyPropertyChangedEventArgs e) => ((ComplexTorusView)sender).QueueUpdate();
    private static void GridChanged(DependencyObject sender, DependencyPropertyChangedEventArgs e)
    {
        var view = (ComplexTorusView)sender;
        if (view.LatticePlot == null) return;
        view.LatticePlot.ShowGrid = view.TorusPlot.ShowGrid = view.ShowGrid;
    }

    private void QueueUpdate()
    {
        if (updateQueued || disposed) return;
        updateQueued = true;
        Dispatcher.BeginInvoke(DispatcherPriority.DataBind, new Action(() =>
        {
            updateQueued = false;
            if (!disposed) UpdateModel(IsLoaded && IsVisible);
        }));
    }

    private void UpdateModel(bool active)
    {
        if (!disposed && Snapshot != null) Model.Update(Snapshot.Curve, Samples ?? Array.Empty<EllipticCurvePoint>(), active);
    }

    public void Fit() => TorusPlot.Fit();
    public void Zoom(double factor) => TorusPlot.Zoom(factor);
    public void Dispose()
    {
        disposed = true;
        Model.Dispose();
    }
}

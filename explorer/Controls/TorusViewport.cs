using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Media3D;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Controls;

/// <summary>A rotatable torus whose markers share period coordinates with PeriodLatticePlot.</summary>
public sealed class TorusViewport : Grid
{
    public static readonly DependencyProperty PointsProperty = DependencyProperty.Register(nameof(Points), typeof(IReadOnlyList<TorusPoint>), typeof(TorusViewport), new PropertyMetadata(null, PointsChanged));
    public static readonly DependencyProperty SelectedPointProperty = DependencyProperty.Register(nameof(SelectedPoint), typeof(TorusPoint), typeof(TorusViewport), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.BindsTwoWayByDefault, SelectionChanged));
    public static readonly DependencyProperty ShowGridProperty = DependencyProperty.Register(nameof(ShowGrid), typeof(bool), typeof(TorusViewport), new PropertyMetadata(true, GridChanged));
    private static readonly MeshGeometry3D MarkerMesh = TorusMesh.Sphere();
    private static readonly Material Gold = TorusMesh.Material("#F7CC7C", true);
    private static readonly Material White = TorusMesh.Material("#EDF3F7", true);
    private readonly Viewport3D viewport = new();
    private readonly OrthographicCamera camera = new() { NearPlaneDistance = 0.1, FarPlaneDistance = 50 };
    private readonly Model3DGroup fineGrid = new();
    private readonly Model3DGroup markers = new();
    private readonly ModelVisual3D gridVisual = new();
    private readonly Dictionary<GeometryModel3D, (TorusPoint Point, ScaleTransform3D Scale)> markerPoints = new();
    private double azimuth = TorusCameraState.Default.Azimuth, elevation = TorusCameraState.Default.Elevation,
        span = TorusCameraState.Default.Span;
    private Point? pressedAt, previousMouse;
    private bool dragged;

    public IReadOnlyList<TorusPoint>? Points { get => (IReadOnlyList<TorusPoint>?)GetValue(PointsProperty); set => SetValue(PointsProperty, value); }
    public TorusPoint? SelectedPoint { get => (TorusPoint?)GetValue(SelectedPointProperty); set => SetValue(SelectedPointProperty, value); }
    public bool ShowGrid { get => (bool)GetValue(ShowGridProperty); set => SetValue(ShowGridProperty, value); }

    public TorusViewport()
    {
        Background = Brushes.Transparent;
        ClipToBounds = true;
        Focusable = true;
        Cursor = Cursors.Hand;
        viewport.Camera = camera;
        Children.Add(viewport);
        var lighting = new Model3DGroup();
        lighting.Children.Add(new AmbientLight(Color.FromRgb(95, 115, 125)));
        lighting.Children.Add(new DirectionalLight(Colors.White, new Vector3D(-2, -3, -4)));
        lighting.Children.Add(new DirectionalLight(Color.FromRgb(95, 160, 160), new Vector3D(2, 1, 3)));
        lighting.Freeze();
        viewport.Children.Add(new ModelVisual3D { Content = lighting });
        var surface = new GeometryModel3D(TorusMesh.Surface(), TorusMesh.Material("#245C59"));
        surface.Freeze();
        viewport.Children.Add(new ModelVisual3D { Content = surface });

        var gridMaterial = TorusMesh.Material("#41716F", true);
        for (var i = 1; i < 12; i++)
        {
            fineGrid.Children.Add(new GeometryModel3D(TorusMesh.Cycle(true, i / 12.0, 0.0009), gridMaterial));
            fineGrid.Children.Add(new GeometryModel3D(TorusMesh.Cycle(false, i / 12.0, 0.0009), gridMaterial));
        }
        fineGrid.Freeze();
        gridVisual.Content = fineGrid;
        viewport.Children.Add(gridVisual);
        var cycles = new Model3DGroup();
        cycles.Children.Add(new GeometryModel3D(TorusMesh.Cycle(true, 0, 0.003), TorusMesh.Material("#63E6CF", true)));
        cycles.Children.Add(new GeometryModel3D(TorusMesh.Cycle(false, 0, 0.002), TorusMesh.Material("#84AAFF", true)));
        cycles.Freeze();
        viewport.Children.Add(new ModelVisual3D { Content = cycles });
        viewport.Children.Add(new ModelVisual3D { Content = markers });
        SizeChanged += (_, _) => UpdateCamera();
        UpdateCamera();
    }

    public void Fit() => RestoreCamera(TorusCameraState.Default);

    public void Zoom(double factor)
    {
        span = Math.Clamp(span * factor, 3.5, 18);
        UpdateCamera();
    }

    private void UpdateCamera()
    {
        var a = azimuth * Math.PI / 180;
        var b = elevation * Math.PI / 180;
        var position = new Point3D(12 * Math.Cos(b) * Math.Sin(a), 12 * Math.Sin(b), 12 * Math.Cos(b) * Math.Cos(a));
        camera.Position = position;
        camera.LookDirection = new Vector3D(-position.X, -position.Y, -position.Z);
        camera.UpDirection = new Vector3D(0, 1, 0);
        camera.Width = span * Math.Max(1, ActualWidth / Math.Max(1, ActualHeight));
    }

    internal TorusCameraState CaptureCamera() => new(azimuth, elevation, span);
    internal void RestoreCamera(TorusCameraState state)
    {
        azimuth = state.Azimuth;
        elevation = state.Elevation;
        span = state.Span;
        UpdateCamera();
    }

    private static void GridChanged(DependencyObject sender, DependencyPropertyChangedEventArgs e)
    {
        var view = (TorusViewport)sender;
        view.gridVisual.Content = view.ShowGrid ? view.fineGrid : null;
    }

    private static void PointsChanged(DependencyObject sender, DependencyPropertyChangedEventArgs e) => ((TorusViewport)sender).RebuildMarkers();
    private static void SelectionChanged(DependencyObject sender, DependencyPropertyChangedEventArgs e) => ((TorusViewport)sender).UpdateSelection();

    private void RebuildMarkers()
    {
        markers.Children.Clear();
        markerPoints.Clear();
        if (Points == null) return;
        foreach (var point in Points)
        {
            var position = TorusMesh.Position(point.Coordinates.U, point.Coordinates.V, 0.045);
            var scale = new ScaleTransform3D(0.065, 0.065, 0.065);
            var transform = new Transform3DGroup();
            transform.Children.Add(scale);
            transform.Children.Add(new TranslateTransform3D(position.X, position.Y, position.Z));
            var marker = new GeometryModel3D(MarkerMesh, Gold) { Transform = transform };
            markers.Children.Add(marker);
            markerPoints.Add(marker, (point, scale));
        }
        UpdateSelection();
    }

    private void UpdateSelection()
    {
        foreach (var (marker, data) in markerPoints)
        {
            var selected = Equals(data.Point, SelectedPoint);
            marker.Material = selected ? White : Gold;
            data.Scale.ScaleX = data.Scale.ScaleY = data.Scale.ScaleZ = selected ? 0.105 : 0.065;
        }
    }

    private TorusPoint? Hit(Point position)
    {
        TorusPoint? point = null;
        VisualTreeHelper.HitTest(viewport, null, result =>
        {
            if (result is not RayMeshGeometry3DHitTestResult mesh) return HitTestResultBehavior.Continue;
            if (mesh.ModelHit is GeometryModel3D model && markerPoints.TryGetValue(model, out var hit)) point = hit.Point;
            // The surface occludes markers on its far side.
            return HitTestResultBehavior.Stop;
        }, new PointHitTestParameters(position));
        return point;
    }

    protected override void OnMouseLeftButtonDown(MouseButtonEventArgs e)
    {
        base.OnMouseLeftButtonDown(e);
        Focus();
        if (e.ClickCount == 2) Fit();
        pressedAt = previousMouse = e.GetPosition(this);
        dragged = false;
        CaptureMouse();
        e.Handled = true;
    }

    protected override void OnMouseMove(MouseEventArgs e)
    {
        base.OnMouseMove(e);
        var current = e.GetPosition(this);
        if (pressedAt is { } start && previousMouse is { } previous && IsMouseCaptured)
        {
            if ((current - start).Length > 4) dragged = true;
            if (dragged)
            {
                azimuth = (azimuth - (current.X - previous.X) * 0.5) % 360;
                elevation = Math.Clamp(elevation + (current.Y - previous.Y) * 0.5, -80, 80);
                UpdateCamera();
            }
            previousMouse = current;
            ToolTip = null;
            e.Handled = true;
        }
        else ToolTip = Hit(e.GetPosition(viewport))?.DisplayName;
    }

    protected override void OnMouseLeftButtonUp(MouseButtonEventArgs e)
    {
        base.OnMouseLeftButtonUp(e);
        if (pressedAt == null) return;
        if (!dragged && Hit(e.GetPosition(viewport)) is { } point) SetCurrentValue(SelectedPointProperty, point);
        ReleaseMouseCapture();
        pressedAt = previousMouse = null;
        e.Handled = true;
    }

    protected override void OnLostMouseCapture(MouseEventArgs e)
    {
        pressedAt = previousMouse = null;
        base.OnLostMouseCapture(e);
    }

    protected override void OnMouseWheel(MouseWheelEventArgs e)
    {
        Zoom(Math.Pow(1.2, -e.Delta / 120.0));
        e.Handled = true;
    }

    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        switch (e.Key)
        {
            case Key.Left: azimuth -= 8; break;
            case Key.Right: azimuth += 8; break;
            case Key.Up: elevation = Math.Min(80, elevation + 8); break;
            case Key.Down: elevation = Math.Max(-80, elevation - 8); break;
            case Key.Home: Fit(); break;
            case Key.Add: case Key.OemPlus: Zoom(1 / 1.2); break;
            case Key.Subtract: case Key.OemMinus: Zoom(1.2); break;
            default: return;
        }
        UpdateCamera();
        e.Handled = true;
    }
}

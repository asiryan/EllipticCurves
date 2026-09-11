using System.Globalization;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Controls;

/// <summary>The fundamental parallelogram of Λ/ω₁, with opposite edges identified.</summary>
public sealed class PeriodLatticePlot : FrameworkElement
{
    public static readonly DependencyProperty LatticeProperty = DependencyProperty.Register(nameof(Lattice), typeof(TorusLattice), typeof(PeriodLatticePlot), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.AffectsRender));
    public static readonly DependencyProperty PointsProperty = DependencyProperty.Register(nameof(Points), typeof(IReadOnlyList<TorusPoint>), typeof(PeriodLatticePlot), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.AffectsRender));
    public static readonly DependencyProperty SelectedPointProperty = DependencyProperty.Register(nameof(SelectedPoint), typeof(TorusPoint), typeof(PeriodLatticePlot), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.AffectsRender | FrameworkPropertyMetadataOptions.BindsTwoWayByDefault));
    public static readonly DependencyProperty ShowGridProperty = DependencyProperty.Register(nameof(ShowGrid), typeof(bool), typeof(PeriodLatticePlot), new FrameworkPropertyMetadata(true, FrameworkPropertyMetadataOptions.AffectsRender));
    private static readonly Brush Muted = Brush("#91A1B2"), Text = Brush("#EDF3F7"), First = Brush("#63E6CF"), Second = Brush("#84AAFF"), Gold = Brush("#F7CC7C");
    public TorusLattice? Lattice { get => (TorusLattice?)GetValue(LatticeProperty); set => SetValue(LatticeProperty, value); }
    public IReadOnlyList<TorusPoint>? Points { get => (IReadOnlyList<TorusPoint>?)GetValue(PointsProperty); set => SetValue(PointsProperty, value); }
    public TorusPoint? SelectedPoint { get => (TorusPoint?)GetValue(SelectedPointProperty); set => SetValue(SelectedPointProperty, value); }
    public bool ShowGrid { get => (bool)GetValue(ShowGridProperty); set => SetValue(ShowGridProperty, value); }

    public PeriodLatticePlot()
    {
        Focusable = true;
        Cursor = Cursors.Hand;
        ClipToBounds = true;
    }

    private Point Position(double u, double v)
    {
        var lattice = Lattice!;
        var width = Math.Max(1, ActualWidth - 78);
        var height = Math.Max(1, ActualHeight - 76);
        var span = 1 + lattice.TauReal;
        var scale = Math.Min(width / span, height / lattice.TauImaginary);
        var left = (ActualWidth - span * scale) / 2;
        var bottom = (ActualHeight + lattice.TauImaginary * scale) / 2;
        return new Point(left + (u + v * lattice.TauReal) * scale, bottom - v * lattice.TauImaginary * scale);
    }

    protected override void OnRender(DrawingContext dc)
    {
        dc.DrawRectangle(Brushes.Transparent, null, new Rect(RenderSize));
        if (Lattice == null || ActualWidth < 100 || ActualHeight < 90) return;
        var origin = Position(0, 0);
        var a = Position(1, 0);
        var b = Position(0, 1);
        var opposite = Position(1, 1);
        var fill = new StreamGeometry();
        using (var context = fill.Open())
        {
            context.BeginFigure(origin, true, true);
            context.PolyLineTo(new[] { a, opposite, b }, true, false);
        }
        dc.DrawGeometry(Brush("#1263E6CF"), null, fill);
        if (ShowGrid)
        {
            var grid = new Pen(Brush("#293B45"), 0.7);
            for (var i = 1; i < 8; i++)
            {
                var fraction = i / 8.0;
                dc.DrawLine(grid, Position(fraction, 0), Position(fraction, 1));
                dc.DrawLine(grid, Position(0, fraction), Position(1, fraction));
            }
        }
        var axis = new Pen(Brush("#40515E"), 1);
        dc.DrawLine(axis, new Point(12, origin.Y), new Point(ActualWidth - 12, origin.Y));
        dc.DrawLine(axis, new Point(origin.X, ActualHeight - 12), new Point(origin.X, 12));
        Label(dc, "Re", new Point(ActualWidth - 27, origin.Y + 7), Muted, 10);
        Label(dc, "Im", new Point(origin.X + 5, 8), Muted, 10);
        dc.DrawLine(new Pen(First, 2), origin, a);
        dc.DrawLine(new Pen(Second, 2), origin, b);
        dc.DrawLine(new Pen(First, 1.5) { DashStyle = DashStyles.Dash }, b, opposite);
        dc.DrawLine(new Pen(Second, 1.5) { DashStyle = DashStyles.Dash }, a, opposite);
        Label(dc, "1", a + new Vector(-3, 9), First);
        Label(dc, "τ", b + new Vector(-15, -23), Second);
        Label(dc, "1 + τ", opposite + new Vector(-19, -23), Muted, 10);
        if (Points == null) return;
        foreach (var point in Points.Where(point => !Equals(point, SelectedPoint))) DrawPoint(dc, point, false);
        if (SelectedPoint != null) DrawPoint(dc, SelectedPoint, true);
    }

    private IEnumerable<Point> Copies(TorusPoint point)
    {
        var (u, v) = point.Coordinates;
        yield return Position(u, v);
        if (u < 1e-9) yield return Position(1, v);
        if (v < 1e-9) yield return Position(u, 1);
        if (u < 1e-9 && v < 1e-9) yield return Position(1, 1);
    }

    private void DrawPoint(DrawingContext dc, TorusPoint point, bool selected)
    {
        foreach (var location in Copies(point))
        {
            if (selected) dc.DrawEllipse(Brush("#2663E6CF"), null, location, 10, 10);
            dc.DrawEllipse(Brush("#121920"), new Pen(selected ? Text : Gold, selected ? 2 : 1.5), location, selected ? 5 : 3.5, selected ? 5 : 3.5);
        }
        if (selected)
            Label(dc, point.Name, Position(point.Coordinates.U, point.Coordinates.V) + new Vector(9, -20), Text);
    }

    private TorusPoint? Hit(Point mouse) => Lattice == null ? null : Points?
        .Select(point => (Point: point, Distance: Copies(point).Min(position => (position - mouse).Length)))
        .Where(hit => hit.Distance <= 12).OrderBy(hit => hit.Distance).Select(hit => hit.Point).FirstOrDefault();

    protected override void OnMouseLeftButtonDown(MouseButtonEventArgs e)
    {
        base.OnMouseLeftButtonDown(e);
        Focus();
        if (Hit(e.GetPosition(this)) is { } point) SetCurrentValue(SelectedPointProperty, point);
        e.Handled = true;
    }

    protected override void OnMouseMove(MouseEventArgs e)
    {
        base.OnMouseMove(e);
        ToolTip = Hit(e.GetPosition(this))?.DisplayName;
    }

    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        if (Points == null || Points.Count == 0) return;
        var direction = e.Key is Key.Left or Key.Down ? -1 : e.Key is Key.Right or Key.Up ? 1 : 0;
        if (direction == 0 && e.Key != Key.Home) return;
        var index = Array.IndexOf(Points.ToArray(), SelectedPoint);
        SetCurrentValue(SelectedPointProperty, Points[e.Key == Key.Home ? 0 : (index + direction + Points.Count) % Points.Count]);
        e.Handled = true;
    }

    private void Label(DrawingContext dc, string text, Point location, Brush color, double size = 12) =>
        dc.DrawText(new FormattedText(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight,
            new Typeface("Consolas"), size, color, VisualTreeHelper.GetDpi(this).PixelsPerDip), location);

    private static Brush Brush(string color)
    {
        var brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color));
        brush.Freeze();
        return brush;
    }
}

using System.Globalization;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.Controls;

/// <summary>A real-locus renderer with adaptive axis scales, root-aware sampling and pointer navigation.</summary>
public sealed class CurvePlot : FrameworkElement
{
    public static readonly DependencyProperty SnapshotProperty = DependencyProperty.Register(nameof(Snapshot), typeof(CurveSnapshot), typeof(CurvePlot), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.AffectsRender, InvalidateGeometry));
    public static readonly DependencyProperty SamplesProperty = DependencyProperty.Register(nameof(Samples), typeof(IReadOnlyList<EllipticCurvePoint>), typeof(CurvePlot), new FrameworkPropertyMetadata(null, FrameworkPropertyMetadataOptions.AffectsRender));
    public static readonly DependencyProperty ShowGridProperty = DependencyProperty.Register(nameof(ShowGrid), typeof(bool), typeof(CurvePlot), new FrameworkPropertyMetadata(true, FrameworkPropertyMetadataOptions.AffectsRender));
    public static readonly DependencyProperty ShowPointsProperty = DependencyProperty.Register(nameof(ShowPoints), typeof(bool), typeof(CurvePlot), new FrameworkPropertyMetadata(true, FrameworkPropertyMetadataOptions.AffectsRender));
    private static readonly Brush CurveBrush = Brush("#63E6CF");
    private static readonly Brush PointBrush = Brush("#F7CC7C");
    private static readonly Brush LabelBrush = Brush("#EDF3F7");
    private static readonly Typeface LabelTypeface = new("Consolas");
    private readonly List<(StreamGeometry Upper, StreamGeometry Lower, StreamGeometry Fill)> geometry = new();
    private bool geometryDirty = true;
    private double centerX = PlotViewState.Default.CenterX, centerY = PlotViewState.Default.CenterY,
        verticalSpan = PlotViewState.Default.VerticalSpan, horizontalScaleRatio = PlotViewState.Default.HorizontalScaleRatio;
    private Point? dragPosition, pointer;
    internal event Action? ViewChanged;

    public CurveSnapshot? Snapshot { get => (CurveSnapshot?)GetValue(SnapshotProperty); set => SetValue(SnapshotProperty, value); }
    public IReadOnlyList<EllipticCurvePoint>? Samples { get => (IReadOnlyList<EllipticCurvePoint>?)GetValue(SamplesProperty); set => SetValue(SamplesProperty, value); }
    public bool ShowGrid { get => (bool)GetValue(ShowGridProperty); set => SetValue(ShowGridProperty, value); }
    public bool ShowPoints { get => (bool)GetValue(ShowPointsProperty); set => SetValue(ShowPointsProperty, value); }
    private Rect PlotBounds => new(48, 16, Math.Max(1, ActualWidth - 66), Math.Max(1, ActualHeight - 46));
    private double ScaleY => PlotBounds.Height / verticalSpan;
    private double ScaleX => ScaleY * horizontalScaleRatio;

    public CurvePlot()
    {
        Focusable = true;
        Cursor = Cursors.Cross;
        ClipToBounds = true;
        SizeChanged += (_, _) => RefreshGeometry();
    }

    public void Fit()
    {
        if (Snapshot?.Plot.IsDrawable != true) return;
        RestoreView(GetResetView());
    }

    internal PlotViewState GetResetView() => GetResetView(Snapshot?.Plot);

    internal void Fit(CurvePlotData data)
    {
        if (data.IsDrawable) RestoreView(GetResetView(data));
    }

    internal bool NeedsRefit(CurvePlotData previous, CurvePlotData current)
    {
        if (!current.IsDrawable) return false;
        if (!previous.IsDrawable) return true;
        // Compare the curves' natural fitted ranges, not the user's zoom/pan.
        // Local edits preserve a chosen view, including a deliberately empty one.
        var before = GetResetView(previous);
        var after = GetResetView(current);
        var widthBefore = before.VerticalSpan / before.HorizontalScaleRatio;
        var widthAfter = after.VerticalSpan / after.HorizontalScaleRatio;
        static bool DifferentScale(double a, double b) => a / b is < 0.125 or > 8;
        return DifferentScale(after.VerticalSpan, before.VerticalSpan) || DifferentScale(widthAfter, widthBefore)
            || Math.Abs(after.CenterX - before.CenterX) > 4 * Math.Max(widthBefore, widthAfter) * PlotBounds.Width / PlotBounds.Height
            || Math.Abs(after.CenterY - before.CenterY) > 4 * Math.Max(before.VerticalSpan, after.VerticalSpan);
    }

    private PlotViewState GetResetView(CurvePlotData? data)
    {
        if (data is null || !data.IsDrawable) return PlotViewState.Default;
        var left = data.Roots[0];
        var right = data.Roots[^1];
        // Include part of the unbounded branch even when there is only one real root.
        right += Math.Max(0.5, 0.25 * Math.Max(right - left, data.CharacteristicScale));
        if (!double.IsFinite(right) || !double.IsFinite(right - left)) return PlotViewState.Default;
        var fittedCenterX = left / 2 + right / 2;
        var ys = new List<double> { data.CenterY(left), data.CenterY(right) };
        if (ys.Any(y => !double.IsFinite(y))) return PlotViewState.Default;
        for (var i = 0; i <= 100; i++)
            if (data.TryEvaluate(left + (right - left) * i / 100, out var upper, out var lower)) { ys.Add(upper); ys.Add(lower); }
        var fittedCenterY = ys.Min() / 2 + ys.Max() / 2;
        var horizontalSpan = Math.Max(5, (right - left) * 1.8);
        var desiredVerticalSpan = Math.Max(3.4, (ys.Max() - ys.Min()) * 1.35);
        var aspect = PlotBounds.Height / PlotBounds.Width;
        // x and y can differ by many orders of magnitude (x~s², y~s³).
        // Keep equal units for ordinary views, but fit each axis when one would
        // otherwise collapse. Tick labels and point coordinates stay in the input model.
        var ratio = desiredVerticalSpan / horizontalSpan / aspect;
        if (ratio is >= 0.125 and <= 8) ratio = 1;
        var minimumSpan = MinimumVerticalSpan(fittedCenterX, fittedCenterY, ratio);
        if (!double.IsFinite(ratio) || ratio <= 0 || !double.IsFinite(minimumSpan)
            || minimumSpan > CurvePlotData.MaximumPlotMagnitude) return PlotViewState.Default;
        var fittedSpan = Math.Clamp(Math.Max(desiredVerticalSpan, horizontalSpan * aspect * ratio), minimumSpan, CurvePlotData.MaximumPlotMagnitude);
        return new(fittedCenterX, fittedCenterY, fittedSpan, ratio);
    }

    public void Zoom(double factor) => ZoomAt(factor, new Point(PlotBounds.Left + PlotBounds.Width / 2, PlotBounds.Top + PlotBounds.Height / 2));

    internal PlotViewState CaptureView() => new(centerX, centerY, verticalSpan, horizontalScaleRatio);
    internal void RestoreView(PlotViewState state)
    {
        centerX = state.CenterX;
        centerY = state.CenterY;
        verticalSpan = state.VerticalSpan;
        horizontalScaleRatio = state.HorizontalScaleRatio;
        pointer = dragPosition = null;
        RefreshGeometry();
    }

    public Point ToScreen(double x, double y) => new(PlotBounds.Left + PlotBounds.Width / 2 + (x - centerX) * ScaleX, PlotBounds.Top + PlotBounds.Height / 2 - (y - centerY) * ScaleY);
    public Point ToWorld(Point screen) => new(centerX + (screen.X - PlotBounds.Left - PlotBounds.Width / 2) / ScaleX, centerY - (screen.Y - PlotBounds.Top - PlotBounds.Height / 2) / ScaleY);

    protected override void OnRender(DrawingContext dc)
    {
        base.OnRender(dc);
        dc.DrawRectangle(Brush("#121920"), null, new Rect(RenderSize));
        if (ActualWidth < 90 || ActualHeight < 70 || Snapshot is null) return;
        var bounds = PlotBounds;
        DrawGrid(dc);
        dc.PushClip(new RectangleGeometry(bounds));
        if (geometryDirty) BuildGeometry();
        foreach (var piece in geometry)
        {
            dc.DrawGeometry(Brush("#0963E6CF"), null, piece.Fill);
            dc.DrawGeometry(null, new Pen(Brush("#1263E6CF"), 8), piece.Upper);
            dc.DrawGeometry(null, new Pen(Brush("#1263E6CF"), 8), piece.Lower);
            dc.DrawGeometry(null, new Pen(CurveBrush, 2), piece.Upper);
            dc.DrawGeometry(null, new Pen(CurveBrush, 2), piece.Lower);
        }
        if (Snapshot.IsSingular && Snapshot.Plot.IsDrawable)
            foreach (var root in Snapshot.Plot.Roots)
                dc.DrawEllipse(CurveBrush, null, ToScreen(root, Snapshot.Plot.CenterY(root)), 3, 3);
        if (ShowPoints && Samples is not null)
            foreach (var point in Samples)
            {
                var screen = ToScreen(CurvePlotData.ToDouble(point.X), CurvePlotData.ToDouble(point.Y));
                if (!bounds.Contains(screen)) continue;
                dc.DrawEllipse(Brush("#121920"), new Pen(PointBrush, 1.6), screen, 4, 4);
            }
        dc.Pop();
        DrawPointer(dc);
    }

    private void DrawGrid(DrawingContext dc)
    {
        var bounds = PlotBounds;
        var from = ToWorld(bounds.BottomLeft);
        var to = ToWorld(bounds.TopRight);
        var stepX = NiceStep(80 / ScaleX);
        var stepY = NiceStep(80 / ScaleY);
        if (ShowGrid)
        {
            var minorX = stepX / 5;
            var minorY = stepY / 5;
            var pen = new Pen(Brush("#1C2832"), 0.6);
            for (var x = Math.Ceiling(from.X / minorX) * minorX; x <= to.X; x += minorX)
            { var px = ToScreen(x, 0).X; dc.DrawLine(pen, new Point(px, bounds.Top), new Point(px, bounds.Bottom)); }
            for (var y = Math.Ceiling(from.Y / minorY) * minorY; y <= to.Y; y += minorY)
            { var py = ToScreen(0, y).Y; dc.DrawLine(pen, new Point(bounds.Left, py), new Point(bounds.Right, py)); }
        }
        var majorPen = new Pen(Brush("#2B3946"), 0.7);
        for (var x = Math.Ceiling(from.X / stepX) * stepX; x <= to.X; x += stepX)
        {
            var px = ToScreen(x, 0).X;
            if (ShowGrid) dc.DrawLine(majorPen, new Point(px, bounds.Top), new Point(px, bounds.Bottom));
            var text = Label(Number(x, stepX), 10);
            // Leave the axis name readable even when tick labels use exponents.
            if (px + text.Width / 2 < bounds.Right - 18)
                dc.DrawText(text, new Point(px - text.Width / 2, bounds.Bottom + 9));
        }
        for (var y = Math.Ceiling(from.Y / stepY) * stepY; y <= to.Y; y += stepY)
        {
            var py = ToScreen(0, y).Y;
            if (ShowGrid) dc.DrawLine(majorPen, new Point(bounds.Left, py), new Point(bounds.Right, py));
            var text = Label(Number(y, stepY), 10);
            dc.DrawText(text, new Point(bounds.Left - text.Width - 10, py - text.Height / 2));
        }
        var origin = ToScreen(0, 0);
        var axisPen = new Pen(Brush("#526171"), 1);
        if (origin.X >= bounds.Left && origin.X <= bounds.Right) dc.DrawLine(axisPen, new Point(origin.X, bounds.Top), new Point(origin.X, bounds.Bottom));
        if (origin.Y >= bounds.Top && origin.Y <= bounds.Bottom) dc.DrawLine(axisPen, new Point(bounds.Left, origin.Y), new Point(bounds.Right, origin.Y));
        dc.DrawText(Label("x", 12), new Point(bounds.Right - 6, bounds.Bottom + 9));
        dc.DrawText(Label("y", 12), new Point(bounds.Left - 20, bounds.Top - 15));
        if (horizontalScaleRatio != 1)
        {
            var note = Label("Independent axis scales", 10);
            dc.DrawText(note, new Point(bounds.Right - note.Width, bounds.Top - 15));
        }
    }

    private void BuildGeometry()
    {
        geometry.Clear();
        geometryDirty = false;
        if (Snapshot is null || !Snapshot.Plot.IsDrawable) return;
        var data = Snapshot.Plot;
        var bounds = PlotBounds;
        var left = ToWorld(bounds.TopLeft).X;
        var right = ToWorld(bounds.TopRight).X;
        var cuts = new[] { left }.Concat(data.Roots.Where(x => x > left && x < right)).Append(right).ToArray();
        for (var interval = 1; interval < cuts.Length; interval++)
        {
            var start = cuts[interval - 1];
            var end = cuts[interval];
            if (!data.TryEvaluate((start + end) / 2, out _, out _)) continue;
            var count = Math.Clamp((int)Math.Ceiling((end - start) * ScaleX), 32, 4096);
            var upper = new List<Point>();
            var lower = new List<Point>();
            for (var i = 0; i <= count; i++)
            {
                // Cosine spacing resolves vertical tangents at interval endpoints.
                var x = start + (end - start) * (1 - Math.Cos(Math.PI * i / count)) / 2;
                if (!data.TryEvaluate(x, out var high, out var low)) continue;
                upper.Add(BoundedScreen(x, high));
                lower.Add(BoundedScreen(x, low));
            }
            if (upper.Count < 2) continue;
            var fill = Path(upper.Concat(lower.AsEnumerable().Reverse()).ToArray(), true);
            geometry.Add((Path(upper, false), Path(lower, false), fill));
        }
    }

    private void DrawPointer(DrawingContext dc)
    {
        if (pointer is not Point mouse || dragPosition is not null || Snapshot is null || !PlotBounds.Contains(mouse)) return;
        var world = ToWorld(mouse);
        Point? selected = null;
        string? caption = null;
        var color = CurveBrush;
        if (ShowPoints && Samples is not null)
            foreach (var point in Samples)
            {
                var p = ToScreen(CurvePlotData.ToDouble(point.X), CurvePlotData.ToDouble(point.Y));
                if ((p - mouse).Length > 10) continue;
                selected = p;
                caption = $"Rational point\nx = {point.X}   y = {point.Y}";
                color = PointBrush;
                break;
            }
        if (selected is null && Snapshot.Plot.TryEvaluate(world.X, out var upper, out var lower))
        {
            var y = Math.Abs(world.Y - upper) < Math.Abs(world.Y - lower) ? upper : lower;
            var p = ToScreen(world.X, y);
            if ((p - mouse).Length < 32)
            {
                selected = p;
                caption = FormattableString.Invariant($"Real coordinates (approx.)\nx ≈ {world.X:G6}   y ≈ {y:G6}");
            }
        }
        if (selected is not Point screen || caption is null || !PlotBounds.Contains(screen)) return;
        var pen = new Pen(Brush("#596D7E"), 0.8) { DashStyle = DashStyles.Dash };
        dc.DrawLine(pen, new Point(screen.X, PlotBounds.Top), new Point(screen.X, PlotBounds.Bottom));
        dc.DrawLine(pen, new Point(PlotBounds.Left, screen.Y), new Point(PlotBounds.Right, screen.Y));
        dc.DrawEllipse(Brush("#121920"), new Pen(color, 2), screen, 5, 5);
        var label = Label(caption, 11, Brush("#EDF3F7"));
        var box = new Rect(Math.Clamp(screen.X + 14, PlotBounds.Left, Math.Max(PlotBounds.Left, PlotBounds.Right - label.Width - 24)),
            Math.Clamp(screen.Y - label.Height - 22, PlotBounds.Top, Math.Max(PlotBounds.Top, PlotBounds.Bottom - label.Height - 20)), label.Width + 24, label.Height + 16);
        dc.DrawRoundedRectangle(Brush("#222D37"), new Pen(Brush("#465564"), 1), box, 7, 7);
        dc.DrawText(label, new Point(box.X + 12, box.Y + 8));
    }

    protected override void OnMouseWheel(MouseWheelEventArgs e)
    {
        base.OnMouseWheel(e);
        if (!PlotBounds.Contains(e.GetPosition(this))) return;
        ZoomAt(Math.Pow(1.18, -e.Delta / 120.0), e.GetPosition(this));
        e.Handled = true;
    }

    protected override void OnMouseLeftButtonDown(MouseButtonEventArgs e)
    {
        base.OnMouseLeftButtonDown(e);
        Focus();
        if (e.ClickCount == 2) { Fit(); e.Handled = true; return; }
        if (!PlotBounds.Contains(e.GetPosition(this))) return;
        dragPosition = e.GetPosition(this);
        CaptureMouse();
        Cursor = Cursors.ScrollAll;
        e.Handled = true;
    }

    protected override void OnMouseMove(MouseEventArgs e)
    {
        base.OnMouseMove(e);
        var current = e.GetPosition(this);
        if (dragPosition is Point previous)
        {
            Pan(current.X - previous.X, current.Y - previous.Y);
            dragPosition = current;
        }
        pointer = current;
        InvalidateVisual();
    }

    protected override void OnMouseLeftButtonUp(MouseButtonEventArgs e) { base.OnMouseLeftButtonUp(e); ReleaseMouseCapture(); }
    protected override void OnLostMouseCapture(MouseEventArgs e) { base.OnLostMouseCapture(e); dragPosition = null; Cursor = Cursors.Cross; InvalidateVisual(); }
    protected override void OnMouseLeave(MouseEventArgs e) { base.OnMouseLeave(e); pointer = null; InvalidateVisual(); }
    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        if (e.Key == Key.Home) Fit();
        else if (e.Key is Key.Add or Key.OemPlus) Zoom(1 / 1.25);
        else if (e.Key is Key.Subtract or Key.OemMinus) Zoom(1.25);
        else return;
        e.Handled = true;
    }

    internal void Pan(double deltaX, double deltaY)
    {
        centerX -= deltaX / ScaleX;
        centerY += deltaY / ScaleY;
        RefreshGeometry();
    }

    internal void ZoomAt(double factor, Point anchor)
    {
        if (!double.IsFinite(factor) || factor <= 0) return;
        var before = ToWorld(anchor);
        var minimum = MinimumVerticalSpan(centerX, centerY, horizontalScaleRatio);
        if (!double.IsFinite(minimum) || minimum > CurvePlotData.MaximumPlotMagnitude) return;
        verticalSpan = Math.Clamp(verticalSpan * factor, minimum, CurvePlotData.MaximumPlotMagnitude);
        var after = ToWorld(anchor);
        centerX += before.X - after.X;
        centerY += before.Y - after.Y;
        RefreshGeometry();
    }

    private Point BoundedScreen(double x, double y) { var p = ToScreen(x, y); return new Point(Math.Clamp(p.X, -1e7, 1e7), Math.Clamp(p.Y, -1e7, 1e7)); }
    private double MinimumVerticalSpan(double x, double y, double ratio) =>
        Math.Max(Math.Max(1e-12, Math.Abs(y) * 1e-12), Math.Max(1e-12, Math.Abs(x) * 1e-12) * ratio * (PlotBounds.Height / PlotBounds.Width));
    private void RefreshGeometry() { geometryDirty = true; InvalidateVisual(); ViewChanged?.Invoke(); }
    private static void InvalidateGeometry(DependencyObject source, DependencyPropertyChangedEventArgs args) => ((CurvePlot)source).RefreshGeometry();
    private FormattedText Label(string text, double size, Brush? brush = null) => new(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, LabelTypeface, size, brush ?? LabelBrush, VisualTreeHelper.GetDpi(this).PixelsPerDip);
    private static string Number(double value, double step) => (Math.Abs(value) < step * 1e-7 ? 0 : value).ToString("G4", CultureInfo.InvariantCulture);
    private static double NiceStep(double rough) { var power = Math.Pow(10, Math.Floor(Math.Log10(rough))); var value = rough / power; return (value < 2 ? 2 : value < 5 ? 5 : 10) * power; }
    private static SolidColorBrush Brush(string color) { var result = new SolidColorBrush((Color)ColorConverter.ConvertFromString(color)); result.Freeze(); return result; }
    private static StreamGeometry Path(IReadOnlyList<Point> points, bool closed)
    {
        var result = new StreamGeometry();
        using (var context = result.Open()) { context.BeginFigure(points[0], closed, closed); context.PolyLineTo(points.Skip(1).ToArray(), true, false); }
        result.Freeze();
        return result;
    }
}

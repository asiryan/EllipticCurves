using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace EllipticCurves.Explorer.Controls;

internal static class PlotImageExporter
{
    public static RenderTargetBitmap Render(FrameworkElement target)
    {
        target.UpdateLayout();
        var size = target.RenderSize;
        const double scale = 2;
        var bitmap = new RenderTargetBitmap((int)Math.Ceiling(size.Width * scale), (int)Math.Ceiling(size.Height * scale),
            96 * scale, 96 * scale, PixelFormats.Pbgra32);
        var drawing = new DrawingVisual();
        using (var context = drawing.RenderOpen())
        {
            // Cover the entire bitmap, including a fractional layout size rounded up to a pixel.
            context.DrawRectangle(new SolidColorBrush(Color.FromRgb(0x12, 0x19, 0x20)), null,
                new Rect(0, 0, bitmap.PixelWidth / scale, bitmap.PixelHeight / scale));
            // WPF includes the target's offset within its parent when rendering a visual.
            // Map just its own bounds to the image origin, without changing the live layout.
            var offset = VisualTreeHelper.GetOffset(target);
            var brush = new VisualBrush(target)
            {
                ViewboxUnits = BrushMappingMode.Absolute,
                Viewbox = new Rect(new Point(offset.X, offset.Y), size),
                Stretch = Stretch.Fill
            };
            context.DrawRectangle(brush, null, new Rect(size));
        }
        bitmap.Render(drawing);
        return bitmap;
    }
}

using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Interop;
using System.Windows.Media;

namespace EllipticCurves.Visualizer.Windowing;

internal static class WindowWorkArea
{
    public static Thickness GetContentMargin(Window window)
    {
        var margin = new Thickness(1);
        if (window.WindowState != WindowState.Maximized) return margin;

        var handle = new WindowInteropHelper(window).Handle;
        if (handle == IntPtr.Zero || !GetWindowRect(handle, out var bounds)) return margin;
        var monitor = MonitorFromWindow(handle, 2); // MONITOR_DEFAULTTONEAREST
        var info = new MonitorInfo { Size = Marshal.SizeOf<MonitorInfo>() };
        if (!GetMonitorInfo(monitor, ref info)) return margin;

        // A maximized HWND extends its invisible resize frame beyond the work area.
        // WindowChrome draws into that frame, so keep content inside the visible area.
        // Measure the actual monitor and convert physical pixels to this window's DIPs.
        var dpi = VisualTreeHelper.GetDpi(window);
        margin.Left += Math.Max(0, info.WorkArea.Left - bounds.Left) / dpi.DpiScaleX;
        margin.Top += Math.Max(0, info.WorkArea.Top - bounds.Top) / dpi.DpiScaleY;
        margin.Right += Math.Max(0, bounds.Right - info.WorkArea.Right) / dpi.DpiScaleX;
        margin.Bottom += Math.Max(0, bounds.Bottom - info.WorkArea.Bottom) / dpi.DpiScaleY;
        return margin;
    }

    [DllImport("user32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool GetWindowRect(IntPtr window, out NativeRect rectangle);

    [DllImport("user32.dll")]
    private static extern IntPtr MonitorFromWindow(IntPtr window, uint flags);

    [DllImport("user32.dll", EntryPoint = "GetMonitorInfoW")]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool GetMonitorInfo(IntPtr monitor, ref MonitorInfo info);
}

using System.Windows;
using System.Windows.Controls.Primitives;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media;

namespace EllipticCurves.Explorer.Controls;

// Shared dismissal for title-bar menus. Automatic Popup dismissal would uncheck
// the toggle before its closing click, causing that click to reopen the menu.
internal sealed class TitleBarPopup
{
    private readonly FrameworkElement host;
    private readonly ToggleButton toggle;
    private readonly Popup popup;
    private Window? owner;
    private HwndSource? source;

    public TitleBarPopup(FrameworkElement host, ToggleButton toggle, Popup popup)
    {
        this.host = host;
        this.toggle = toggle;
        this.popup = popup;
        toggle.Checked += ToggleChanged;
        toggle.Unchecked += ToggleChanged;
        host.Unloaded += (_, _) => { Close(); Detach(); };
        popup.Child.PreviewKeyDown += MenuKeyDown;
    }

    public void Close() => toggle.IsChecked = false;

    private void ToggleChanged(object sender, RoutedEventArgs e)
    {
        Detach();
        if (toggle.IsChecked != true || Window.GetWindow(host) is not { } window) return;
        owner = window;
        owner.PreviewMouseDown += OwnerMouseDown;
        owner.PreviewKeyDown += MenuKeyDown;
        owner.Deactivated += CloseFromOwner;
        owner.LocationChanged += CloseFromOwner;
        owner.SizeChanged += CloseFromOwner;
        var handle = new WindowInteropHelper(owner).Handle;
        source = handle == IntPtr.Zero ? null : HwndSource.FromHwnd(handle);
        source?.AddHook(OwnerMessage);
    }

    private void Detach()
    {
        source?.RemoveHook(OwnerMessage);
        source = null;
        if (owner == null) return;
        owner.PreviewMouseDown -= OwnerMouseDown;
        owner.PreviewKeyDown -= MenuKeyDown;
        owner.Deactivated -= CloseFromOwner;
        owner.LocationChanged -= CloseFromOwner;
        owner.SizeChanged -= CloseFromOwner;
        owner = null;
    }

    private void OwnerMouseDown(object sender, MouseButtonEventArgs e)
    {
        if (e.OriginalSource is Visual visual &&
            (visual == toggle || toggle.IsAncestorOf(visual) ||
             visual == popup.Child || popup.Child.IsAncestorOf(visual))) return;
        Close();
    }

    private void MenuKeyDown(object sender, KeyEventArgs e)
    {
        if (e.Key != Key.Escape || toggle.IsChecked != true) return;
        Close();
        toggle.Focus();
        e.Handled = true;
    }

    private IntPtr OwnerMessage(IntPtr hwnd, int message, IntPtr wParam, IntPtr lParam, ref bool handled)
    {
        // WM_NC[L/R/M/X]BUTTONDOWN bypass WPF input. Leave them unhandled so
        // Windows can still drag/resize the window and open its system menu.
        if (message is 0x00A1 or 0x00A4 or 0x00A7 or 0x00AB) Close();
        return IntPtr.Zero;
    }

    private void CloseFromOwner(object? sender, EventArgs e) => Close();
}

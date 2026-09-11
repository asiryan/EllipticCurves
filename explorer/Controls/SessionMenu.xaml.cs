using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;

namespace EllipticCurves.Explorer.Controls;

public partial class SessionMenu : UserControl
{
    private readonly TitleBarPopup menu;
    public event Action? OpenRequested;
    public event Action? SaveRequested;
    public event Action? ExitRequested;

    public SessionMenu()
    {
        InitializeComponent();
        menu = new TitleBarPopup(this, Toggle, MenuPopup);
    }

    private void OpenClick(object sender, RoutedEventArgs e) { menu.Close(); OpenRequested?.Invoke(); }
    private void SaveClick(object sender, RoutedEventArgs e) { menu.Close(); SaveRequested?.Invoke(); }
    private void ExitClick(object sender, RoutedEventArgs e) { menu.Close(); ExitRequested?.Invoke(); }
    private void PopupOpened(object? sender, EventArgs e) =>
        Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() =>
        {
            if (MenuPopup.IsOpen) (OpenButton.IsEnabled ? OpenButton : SaveButton).Focus();
        }));
}

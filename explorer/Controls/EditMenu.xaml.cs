using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;

namespace EllipticCurves.Explorer.Controls;

public partial class EditMenu : UserControl
{
    private readonly TitleBarPopup menu;
    public event Action? UndoRequested;
    public event Action? RedoRequested;
    public EditMenu()
    {
        InitializeComponent();
        menu = new TitleBarPopup(this, Toggle, MenuPopup);
    }
    internal void Close() => menu.Close();
    private void UndoClick(object sender, RoutedEventArgs e) { Close(); UndoRequested?.Invoke(); }
    private void RedoClick(object sender, RoutedEventArgs e) { Close(); RedoRequested?.Invoke(); }
    private void PopupOpened(object? sender, EventArgs e) =>
        Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() =>
        {
            if (MenuPopup.IsOpen) (UndoButton.IsEnabled ? UndoButton : RedoButton.IsEnabled ? RedoButton : (Control)Toggle).Focus();
        }));
}

using System.Windows;
using System.Windows.Controls;
using System.Windows.Threading;

namespace EllipticCurves.Explorer.Controls;

public partial class HelpMenu : UserControl
{
    private readonly TitleBarPopup menu;
    public event Action? UserGuideRequested;
    public event Action? ShortcutsRequested;
    public event Action? LmfdbRequested;
    public event Action? ProjectRequested;
    public event Action? IssueRequested;
    public event Action? AboutRequested;

    public HelpMenu()
    {
        InitializeComponent();
        menu = new TitleBarPopup(this, Toggle, MenuPopup);
    }

    internal void Close() => menu.Close();
    private void Request(Action? action) { Close(); action?.Invoke(); }
    private void GuideClick(object sender, RoutedEventArgs e) => Request(UserGuideRequested);
    private void ShortcutsClick(object sender, RoutedEventArgs e) => Request(ShortcutsRequested);
    private void LmfdbClick(object sender, RoutedEventArgs e) => Request(LmfdbRequested);
    private void ProjectClick(object sender, RoutedEventArgs e) => Request(ProjectRequested);
    private void IssueClick(object sender, RoutedEventArgs e) => Request(IssueRequested);
    private void AboutClick(object sender, RoutedEventArgs e) => Request(AboutRequested);
    private void PopupOpened(object? sender, EventArgs e) =>
        Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() =>
        {
            if (MenuPopup.IsOpen) GuideButton.Focus();
        }));
}

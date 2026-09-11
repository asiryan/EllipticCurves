using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.ViewModels;

namespace EllipticCurves.Visualizer.Controls;

public partial class ExplorerMenu : UserControl
{
    public event Action<CalculationOperation>? OperationRequested;
    public ExplorerMenu() { InitializeComponent(); DataContext = new ExplorerMenuViewModel(); }
    private void PopupOpened(object? sender, EventArgs e) => Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() => SearchBox.Focus()));
    private void PopupKeyDown(object sender, KeyEventArgs e)
    { if (e.Key == Key.Escape) { Toggle.IsChecked = false; Toggle.Focus(); e.Handled = true; } }
    private void OperationClick(object sender, RoutedEventArgs e)
    {
        if (sender is not Button { Tag: CalculationOperation operation }) return;
        Toggle.IsChecked = false;
        OperationRequested?.Invoke(operation);
    }
}

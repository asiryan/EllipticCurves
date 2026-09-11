using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer.Controls;

public partial class ExplorerMenu : UserControl
{
    public event Action<CalculationOperation>? OperationRequested;
    public ExplorerMenu() { InitializeComponent(); DataContext = new ExplorerMenuViewModel(); }
    private void PopupOpened(object? sender, EventArgs e) => Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() => SearchBox.Focus()));
    private void PopupKeyDown(object sender, KeyEventArgs e)
    { if (e.Key == Key.Escape) { Toggle.IsChecked = false; Toggle.Focus(); e.Handled = true; } }
    private void CategorySelectionChanged(object sender, SelectionChangedEventArgs e)
    {
        if (e.Source != sender) return;
        // Reset after the category binding has replaced and laid out the items.
        Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() =>
        {
            OperationsList.ApplyTemplate();
            if (OperationsList.Template?.FindName("PART_ScrollViewer", OperationsList) is ScrollViewer scroll)
                scroll.ScrollToTop();
        }));
    }
    private void OperationMouseDown(object sender, MouseButtonEventArgs e)
    {
        if (e.ChangedButton != MouseButton.Left || e.OriginalSource is not DependencyObject source
            || ItemsControl.ContainerFromElement(OperationsList, source) is not ListBoxItem row) return;

        // Focusing a partly visible row on mouse-down must not scroll it away
        // before mouse-up. Intercept below the ScrollViewer, for this input turn
        // only; keyboard navigation and wheel/scrollbar commands still scroll.
        RequestBringIntoViewEventHandler suppressFocusScroll = (_, request) => request.Handled = true;
        row.RequestBringIntoView += suppressFocusScroll;
        Dispatcher.BeginInvoke(DispatcherPriority.Input,
            new Action(() => row.RequestBringIntoView -= suppressFocusScroll));
    }
    private void OperationClick(object sender, RoutedEventArgs e)
    {
        if (sender is not Button { Tag: CalculationOperation operation }) return;
        Toggle.IsChecked = false;
        OperationRequested?.Invoke(operation);
    }
}

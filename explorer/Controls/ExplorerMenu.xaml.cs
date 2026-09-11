using System.ComponentModel;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media;
using System.Windows.Threading;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer.Controls;

public partial class ExplorerMenu : UserControl
{
    private const int WmNcLButtonDown = 0x00A1;
    private const int WmNcRButtonDown = 0x00A4;
    private const int WmNcMButtonDown = 0x00A7;
    private const int WmNcXButtonDown = 0x00AB;
    private Window? popupOwner;
    private HwndSource? popupOwnerSource;
    public event Action<CalculationOperation>? OperationRequested;
    public ExplorerMenu()
    {
        InitializeComponent();
        var model = new ExplorerMenuViewModel();
        model.PropertyChanged += OperationsChanged;
        DataContext = model;
    }

    private void ToggleStateChanged(object sender, RoutedEventArgs e)
    {
        DetachPopupOwner();
        if (Toggle.IsChecked != true || Window.GetWindow(this) is not { } owner) return;
        popupOwner = owner;
        owner.PreviewMouseDown += OwnerMouseDown;
        owner.Deactivated += CloseFromOwner;
        owner.LocationChanged += CloseFromOwner;
        owner.SizeChanged += CloseFromOwner;
        var handle = new WindowInteropHelper(owner).Handle;
        popupOwnerSource = handle == IntPtr.Zero ? null : HwndSource.FromHwnd(handle);
        popupOwnerSource?.AddHook(OwnerWindowMessage);
    }

    private void DetachPopupOwner()
    {
        popupOwnerSource?.RemoveHook(OwnerWindowMessage);
        popupOwnerSource = null;
        if (popupOwner == null) return;
        popupOwner.PreviewMouseDown -= OwnerMouseDown;
        popupOwner.Deactivated -= CloseFromOwner;
        popupOwner.LocationChanged -= CloseFromOwner;
        popupOwner.SizeChanged -= CloseFromOwner;
        popupOwner = null;
    }

    private void OwnerMouseDown(object sender, MouseButtonEventArgs e)
    {
        // Let the toggle finish its own click without changing its state first.
        if (e.OriginalSource is Visual source &&
            (source == Toggle || Toggle.IsAncestorOf(source) ||
             source == MenuPopup.Child || MenuPopup.Child.IsAncestorOf(source))) return;
        Toggle.IsChecked = false;
    }

    private IntPtr OwnerWindowMessage(IntPtr hwnd, int message, IntPtr wParam, IntPtr lParam, ref bool handled)
    {
        // WindowChrome routes title-bar and resize-border clicks through native
        // messages instead of PreviewMouseDown. Dismiss without consuming the
        // message, so dragging, resizing and the system menu still work.
        if (message is WmNcLButtonDown or WmNcRButtonDown or WmNcMButtonDown or WmNcXButtonDown)
            Toggle.IsChecked = false;
        return IntPtr.Zero;
    }

    private void CloseFromOwner(object? sender, EventArgs e) => Toggle.IsChecked = false;

    private void MenuUnloaded(object sender, RoutedEventArgs e)
    {
        Toggle.IsChecked = false;
        DetachPopupOwner();
    }

    private void PopupOpened(object? sender, EventArgs e) =>
        Dispatcher.BeginInvoke(DispatcherPriority.Input, new Action(() =>
        {
            if (MenuPopup.IsOpen) SearchBox.Focus();
        }));

    private void PopupKeyDown(object sender, KeyEventArgs e)
    {
        if (e.Key == Key.Escape)
        {
            Toggle.IsChecked = false;
            Toggle.Focus();
            e.Handled = true;
        }
    }

    private void OperationsChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName != nameof(ExplorerMenuViewModel.Operations)) return;
        // Both category and search changes replace the list. Reset after binding/layout.
        Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(ScrollOperationsToTop));
    }

    private void ScrollOperationsToTop()
    {
        OperationsList.ApplyTemplate();
        if (OperationsList.Template?.FindName("PART_ScrollViewer", OperationsList) is ScrollViewer scroll)
            scroll.ScrollToTop();
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

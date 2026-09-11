using System.Windows;
using System.Windows.Automation;
using System.Windows.Input;

namespace EllipticCurves.Explorer;

public partial class ConfirmationWindow : Window
{
    internal bool Confirmed { get; private set; }

    internal ConfirmationWindow(string title, string message, string confirmText, string? detail = null)
    {
        InitializeComponent();
        Title = Heading.Text = title;
        Message.Text = message;
        ConfirmButton.Content = confirmText;
        if (!string.IsNullOrEmpty(detail))
        {
            Detail.Text = detail;
            DetailPanel.Visibility = Visibility.Visible;
        }
    }

    internal static bool Confirm(Window? owner, string title, string message, string confirmText, string? detail = null)
    {
        var dialog = new ConfirmationWindow(title, message, confirmText, detail) { Owner = owner };
        if (owner == null) dialog.WindowStartupLocation = WindowStartupLocation.CenterScreen;
        dialog.ShowDialog();
        return dialog.Confirmed;
    }

    internal static void ShowMessage(Window? owner, string title, string message)
    {
        var dialog = new ConfirmationWindow(title, message, string.Empty) { Owner = owner };
        dialog.ConfirmButton.Visibility = Visibility.Collapsed;
        dialog.CancelButton.Content = "OK";
        AutomationProperties.SetName(dialog.CloseButton, "Close message");
        if (owner == null) dialog.WindowStartupLocation = WindowStartupLocation.CenterScreen;
        dialog.ShowDialog();
    }

    private void WindowContentRendered(object? sender, EventArgs e) => CancelButton.Focus();
    private void ConfirmClick(object sender, RoutedEventArgs e) { Confirmed = true; Close(); }
    private void CancelClick(object sender, RoutedEventArgs e) => Close();
    private void HeaderMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
    {
        if (e.ButtonState == MouseButtonState.Pressed) { DragMove(); e.Handled = true; }
    }
}

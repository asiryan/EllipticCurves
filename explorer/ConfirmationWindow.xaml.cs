using System.Windows;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Input;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer;

internal enum SaveChangesChoice { Cancel, Save, Discard }
internal sealed record SaveChangesResult(SaveChangesChoice Choice, string? FileName = null);

public partial class ConfirmationWindow : Window
{
    internal bool Confirmed { get; private set; }
    internal SaveChangesChoice SaveChoice { get; private set; }
    internal SaveChangesResult SaveResult => new(SaveChoice,
        SaveChoice == SaveChangesChoice.Save ? SessionFileName.Text.Trim() : null);

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
        var dialog = new ConfirmationWindow(title, message, confirmText, detail);
        dialog.ShowForOwner(owner);
        return dialog.Confirmed;
    }

    internal static void ShowMessage(Window? owner, string title, string message)
    {
        var dialog = new ConfirmationWindow(title, message, string.Empty);
        dialog.ConfirmButton.Visibility = Visibility.Collapsed;
        dialog.CancelButton.Content = "OK";
        AutomationProperties.SetName(dialog.CloseButton, "Close message");
        dialog.ShowForOwner(owner);
    }

    internal static ConfirmationWindow CreateSaveChangesDialog(string sessionName)
    {
        var dialog = new ConfirmationWindow(SessionMessages.UnsavedChanges,
            "Save changes to this session before continuing? Choosing Discard will lose the unsaved changes.", "Save")
            { Width = 520 };
        dialog.DiscardButton.Visibility = Visibility.Visible;
        dialog.SessionNamePanel.Visibility = Visibility.Visible;
        dialog.SessionFileName.Text = sessionName;
        return dialog;
    }

    internal static SaveChangesResult AskToSaveChanges(Window owner, string sessionName)
    {
        var dialog = CreateSaveChangesDialog(sessionName);
        dialog.ShowForOwner(owner);
        return dialog.SaveResult;
    }

    private void SessionFileNameChanged(object sender, TextChangedEventArgs e)
    {
        if (SessionNamePanel.Visibility != Visibility.Visible) return;
        var valid = SessionFile.IsValidFileName(SessionFileName.Text);
        ConfirmButton.IsEnabled = valid;
        SessionNameError.Visibility = valid ? Visibility.Collapsed : Visibility.Visible;
    }

    private void ShowForOwner(Window? owner)
    {
        Owner = owner;
        if (owner == null) WindowStartupLocation = WindowStartupLocation.CenterScreen;
        ShowDialog();
    }

    private void WindowContentRendered(object? sender, EventArgs e) => DialogRoot.Focus();
    private void ConfirmClick(object sender, RoutedEventArgs e)
    {
        Confirmed = true;
        SaveChoice = SaveChangesChoice.Save;
        Close();
    }

    private void DiscardClick(object sender, RoutedEventArgs e)
    {
        SaveChoice = SaveChangesChoice.Discard;
        Close();
    }

    private void CancelClick(object sender, RoutedEventArgs e) => Close();

    private void HeaderMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
    {
        if (e.ButtonState == MouseButtonState.Pressed)
        {
            DragMove();
            e.Handled = true;
        }
    }
}

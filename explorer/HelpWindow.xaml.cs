using System.Windows;
using System.Windows.Input;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.Windowing;

namespace EllipticCurves.Explorer;

internal enum HelpPage { About, Shortcuts, License }
public sealed record HelpShortcut(string Description, string Keys);

public partial class HelpWindow : Window
{
    internal HelpPage Page { get; }
    internal event Action? LicenseRequested;

    internal HelpWindow(HelpPage page)
    {
        InitializeComponent();
        Page = page;
        switch (page)
        {
            case HelpPage.About:
                Title = "About EllipticCurves";
                SectionHeading.Text = "HELP · ABOUT";
                AboutContent.Visibility = CopyVersionButton.Visibility = Visibility.Visible;
                break;
            case HelpPage.Shortcuts:
                Title = "Keyboard Shortcuts";
                SectionHeading.Text = "HELP · SHORTCUTS";
                Width = 600;
                Height = 620;
                ShortcutsContent.Visibility = Visibility.Visible;
                ShortcutList.ItemsSource = new HelpShortcut[]
                {
                    new("User Guide", "F1"),
                    new("New session", "Ctrl + N"),
                    new("Open session", "Ctrl + O"),
                    new("Save changes", "Ctrl + S"),
                    new("Save as", "Ctrl + Shift + S"),
                    new("Undo", "Ctrl + Z"),
                    new("Redo", "Ctrl + Y"),
                    new("Reset the active plot view", "Ctrl + F"),
                    new("Finish an equation or coefficient edit", "Enter"),
                    new("Dismiss a menu, Help or import window", "Esc")
                };
                break;
            case HelpPage.License:
                Title = "MIT License";
                SectionHeading.Text = "HELP · ABOUT · LICENSE";
                Width = 760;
                Height = 620;
                LicenseText.Text = ExplorerInfo.LicenseText;
                LicenseContent.Visibility = Visibility.Visible;
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(page));
        }
        Loaded += (_, _) => CloseButton.Focus();
        PreviewKeyDown += (_, e) => { if (e.Key == Key.Escape) { e.Handled = true; Close(); } };
        SourceInitialized += UpdateInsets;
        StateChanged += UpdateInsets;
        SizeChanged += UpdateInsets;
        LocationChanged += UpdateInsets;
        DpiChanged += (_, _) => Dispatcher.BeginInvoke(new Action(() => HelpRoot.Margin = WindowWorkArea.GetContentMargin(this)));
    }

    private void UpdateInsets(object? sender, EventArgs e) => HelpRoot.Margin = WindowWorkArea.GetContentMargin(this);
    private void CloseClick(object sender, RoutedEventArgs e) => Close();
    private void LicenseClick(object sender, RoutedEventArgs e) => LicenseRequested?.Invoke();
    private async void CopyVersionClick(object sender, RoutedEventArgs e)
    {
        if (!ClipboardActions.CopyText(this, ExplorerInfo.VersionInfo, "Copy version info")) return;
        CopyVersionButton.Content = "Copied ✓";
        await Task.Delay(1400);
        CopyVersionButton.Content = "Copy version info";
    }
}

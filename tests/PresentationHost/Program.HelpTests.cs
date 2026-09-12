using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Interop;
using System.Windows.Media.Imaging;
using EllipticCurves;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Controls;

internal static partial class Program
{
    private static void CheckHelp()
    {
        var menu = new HelpMenu();
        var popup = (Popup)menu.FindName("MenuPopup");
        var toggle = (ToggleButton)menu.FindName("Toggle");
        BindingOperations.ClearBinding(popup, Popup.IsOpenProperty);
        var actions = new List<string>();
        menu.UserGuideRequested += () => actions.Add("User Guide");
        menu.ShortcutsRequested += () => actions.Add("Keyboard Shortcuts");
        menu.LmfdbRequested += () => actions.Add("LMFDB Website");
        menu.ProjectRequested += () => actions.Add("Project on GitHub");
        menu.IssueRequested += () => actions.Add("Report an Issue");
        menu.AboutRequested += () => actions.Add("About EllipticCurves");
        var buttons = Descendants(popup.Child).OfType<Button>().ToArray();
        foreach (var button in buttons)
        {
            toggle.IsChecked = true;
            button.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(toggle.IsChecked == false, "A Help action left its menu open.");
        }
        Require(actions.SequenceEqual(new[] { "User Guide", "Keyboard Shortcuts", "LMFDB Website", "Project on GitHub", "Report an Issue", "About EllipticCurves" }),
            "Help actions must dispatch exactly once and preserve their order.");
        Require(buttons.Select(button => (string)button.Content).SequenceEqual(actions), "Help labels differ from the actions they trigger.");
        var main = CreateMainWindow();
        try
        {
            Require(main.InputBindings.OfType<KeyBinding>().Any(binding => binding.Key == Key.F1 && binding.Command == ApplicationCommands.Help)
                && ApplicationCommands.Help.CanExecute(null, main), "F1 must open the user guide from the main window.");
            var help = (HelpMenu)main.FindName("Help");
            BindingOperations.ClearBinding((Popup)help.FindName("MenuPopup"), Popup.IsOpenProperty);
            var helpToggle = (ToggleButton)help.FindName("Toggle");
            var fileToggle = (ToggleButton)((SessionMenu)main.FindName("Session")).FindName("Toggle");
            helpToggle.IsChecked = true;
            fileToggle.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = Mouse.PreviewMouseDownEvent });
            Require(helpToggle.IsChecked == false, "Opening File must dismiss Help.");
        }
        finally { main.Close(); }

        Require(ExplorerInfo.LibraryVersion == typeof(EllipticCurveQ).Assembly.GetName().Version!.ToString(3)
            && ExplorerInfo.VersionInfo.Contains(ExplorerInfo.LibraryVersion), "About must use the loaded library version.");
        Require(ExplorerInfo.LicenseText.StartsWith("MIT License") && ExplorerInfo.LicenseText.Contains(ExplorerInfo.Author)
            && ExplorerInfo.LicenseText.Contains("SOFTWARE"), "The application must bundle the complete project license.");
        var menuRoot = (FrameworkElement)popup.Child;
        menuRoot.Measure(new Size(270, double.PositiveInfinity));
        RenderHelp(menuRoot, menuRoot.DesiredSize, "help-menu");
        foreach (var page in Enum.GetValues<HelpPage>())
        {
            var dialog = new HelpWindow(page);
            try
            {
                var root = (FrameworkElement)dialog.Content;
                if (page == HelpPage.About)
                {
                    var licenses = 0;
                    dialog.LicenseRequested += () => licenses++;
                    ((Button)dialog.FindName("LicenseButton")).RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                    Require(licenses == 1, "MIT License must open the bundled license.");
                }
                if (page == HelpPage.License)
                    Require(((TextBox)dialog.FindName("LicenseText")).Text == ExplorerInfo.LicenseText, "The license window must show the entire embedded license.");
                System.Windows.Documents.TextElement.SetFontSize(root, dialog.FontSize);
                System.Windows.Documents.TextElement.SetFontFamily(root, dialog.FontFamily);
                dialog.Content = null;
                foreach (var size in new[] { new Size(dialog.Width - 2, dialog.Height - 2), new Size(dialog.MinWidth - 2, dialog.MinHeight - 2) })
                {
                    RenderHelp(root, size, $"help-{page}-{size.Width}");
                    var close = (Button)dialog.FindName("CloseButton");
                    Require(new Rect(root.RenderSize).Contains(close.TransformToAncestor(root).TransformBounds(new Rect(close.RenderSize))),
                        "The Help close button must remain visible at the minimum size.");
                    var activeContent = new[] { "AboutContent", "ShortcutsContent", "LicenseContent" }.Select(name => (FrameworkElement)dialog.FindName(name))
                        .Single(content => content.Visibility == Visibility.Visible);
                    Require(activeContent.ActualHeight > 120, "Help content must remain readable at the minimum size.");
                }
            }
            finally { dialog.Close(); }
        }
    }

    private static void RenderHelp(FrameworkElement root, Size size, string name)
    {
        using var source = new HwndSource(new HwndSourceParameters("Help layout")
            { WindowStyle = 0, Width = (int)size.Width, Height = (int)size.Height }) { RootVisual = root };
        root.Measure(size);
        root.Arrange(new Rect(new Point(), size));
        root.UpdateLayout();
        if (Environment.GetEnvironmentVariable("ELLIPTIC_EXPLORER_PREVIEW_DIRECTORY") is not { Length: > 0 } directory) return;
        System.IO.Directory.CreateDirectory(directory);
        var encoder = new PngBitmapEncoder();
        encoder.Frames.Add(BitmapFrame.Create(PlotImageExporter.Render(root)));
        using var output = System.IO.File.Create(System.IO.Path.Combine(directory, name + ".png"));
        encoder.Save(output);
    }
}

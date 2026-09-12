using System.Windows;
using System.Windows.Input;
using EllipticCurves.Explorer.Controls;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private void InitializeHelp()
    {
        Help.UserGuideRequested += () => ApplicationCommands.Help.Execute(null, this);
        Help.ShortcutsRequested += () => OpenHelpWindow(HelpPage.Shortcuts);
        Help.LmfdbRequested += () => BrowserActions.Open(this, ExplorerInfo.LmfdbUrl, "LMFDB Website");
        Help.ProjectRequested += () => BrowserActions.Open(this, ExplorerInfo.RepositoryUrl, "Project on GitHub");
        Help.IssueRequested += () => BrowserActions.Open(this, ExplorerInfo.ReportIssueUrl, "Report an Issue");
        Help.AboutRequested += () => OpenHelpWindow(HelpPage.About);
    }

    private void HelpCommandCanExecute(object sender, CanExecuteRoutedEventArgs e)
    {
        e.CanExecute = true;
        e.Handled = true;
    }

    private void HelpCommandExecuted(object sender, ExecutedRoutedEventArgs e)
    {
        e.Handled = true;
        Session.Close();
        Edit.Close();
        Explorer.Close();
        Help.Close();
        BrowserActions.Open(this, ExplorerInfo.UserGuideUrl, "User Guide");
    }

    private void OpenHelpWindow(HelpPage page)
    {
        var existing = OwnedWindows.OfType<HelpWindow>().FirstOrDefault(window => window.Page == page);
        if (existing != null)
        {
            if (existing.WindowState == WindowState.Minimized) existing.WindowState = WindowState.Normal;
            existing.Activate();
            return;
        }
        var dialog = new HelpWindow(page) { Owner = this };
        dialog.LicenseRequested += () => OpenHelpWindow(HelpPage.License);
        dialog.Show();
    }
}

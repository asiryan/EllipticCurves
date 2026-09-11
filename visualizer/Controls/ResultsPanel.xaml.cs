using System.IO;
using System.Runtime.InteropServices;
using System.Windows;
using System.Windows.Controls;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.ViewModels;
using Microsoft.Win32;

namespace EllipticCurves.Visualizer.Controls;

public partial class ResultsPanel : UserControl
{
    public event Action? HideRequested;
    public event Action<CalculationRequest>? RepeatRequested;
    public ResultsPanel() => InitializeComponent();
    private CalculationJobViewModel? Selected => (DataContext as WorkbenchViewModel)?.Selected;
    private void HideClick(object sender, RoutedEventArgs e) => HideRequested?.Invoke();
    private void RepeatClick(object sender, RoutedEventArgs e) { if (Selected != null) RepeatRequested?.Invoke(Selected.Request); }
    private void HistoryContextMenuOpening(object sender, ContextMenuEventArgs e)
    {
        if (sender is not ComboBoxItem { DataContext: CalculationJobViewModel job, ContextMenu: { } menu }
            || DataContext is not WorkbenchViewModel workbench) return;
        // Capture the row under the pointer, independently of the displayed result.
        menu.DataContext = job;
        var delete = (MenuItem)menu.Items[0];
        delete.IsEnabled = workbench.CanDelete(job);
        // A Click handler nested in Setter.Value produces an invalid WPF BAML
        // connection. Wire it here, once even when this menu is reopened.
        delete.Click -= HistoryDeleteClick;
        delete.Click += HistoryDeleteClick;
    }
    private void HistoryDeleteClick(object sender, RoutedEventArgs e)
    {
        if (sender is MenuItem { DataContext: CalculationJobViewModel job }
            && DataContext is WorkbenchViewModel workbench) workbench.Delete(job);
    }
    private void CopyClick(object sender, RoutedEventArgs e)
    {
        if (Selected == null) return;
        try { Clipboard.SetText(Selected.Report); }
        catch (ExternalException) { MessageBox.Show(Window.GetWindow(this), "The clipboard is busy. Please try again.", "Copy result"); }
    }
    private void SaveClick(object sender, RoutedEventArgs e)
    {
        if (Selected == null) return;
        var dialog = new SaveFileDialog { Filter = "Text report (*.txt)|*.txt", FileName = "elliptic-calculation-" + Selected.StartedAt.ToString("yyyyMMdd-HHmmss") + ".txt", Title = "Save calculation report" };
        if (dialog.ShowDialog(Window.GetWindow(this)) != true) return;
        try { File.WriteAllText(dialog.FileName, Selected.Report); }
        catch (Exception error) when (error is IOException or UnauthorizedAccessException)
        { MessageBox.Show(Window.GetWindow(this), "Could not save the report. Check the destination and try again.", "Save result"); }
    }
}

using System.Windows;
using System.Windows.Input;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer;

public partial class LmfdbImportWindow : Window
{
    internal LmfdbImportViewModel Model { get; }
    public LmfdbCurveFormula? ImportedFormula { get; private set; }

    public LmfdbImportWindow() : this(new LmfdbImportViewModel()) { }
    internal LmfdbImportWindow(LmfdbImportViewModel model)
    {
        InitializeComponent();
        DataContext = Model = model;
        Loaded += (_, _) => ConductorInput.Focus();
        Closed += (_, _) => Model.Dispose();
        PreviewKeyDown += (_, e) => { if (e.Key == Key.Escape) { Close(); e.Handled = true; } };
    }

    private void CloseClick(object sender, RoutedEventArgs e) => Close();
    private void StopClick(object sender, RoutedEventArgs e) => Model.Cancel();
    private async void SearchClick(object sender, RoutedEventArgs e) => await Model.SearchAsync();
    private async void NextClick(object sender, RoutedEventArgs e) => await Model.NextAsync();
    private async void PreviousClick(object sender, RoutedEventArgs e) => await Model.PreviousAsync();
    private void ImportClick(object sender, RoutedEventArgs e)
    {
        if (!Model.CanImport) return;
        ImportedFormula = Model.Selected;
        DialogResult = true;
    }
}

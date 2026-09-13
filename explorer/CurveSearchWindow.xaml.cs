using System.IO;
using System.Text.Json;
using System.Windows;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;
using EllipticCurves.Explorer.Controls;
using Microsoft.Win32;

namespace EllipticCurves.Explorer;

public partial class CurveSearchWindow : Window
{
    public CurveSearchViewModel Model { get; }
    public event Action<CurveSearchCandidate>? OpenCurveRequested;
    public CurveSearchWindow(CurveSearchViewModel model)
    {
        Model = model;
        InitializeComponent();
        DataContext = model;
        Closed += (_, _) => Model.Pause();
    }
    private async void StartClick(object sender, RoutedEventArgs e) => await Model.StartAsync();
    private void PauseClick(object sender, RoutedEventArgs e) => Model.Pause();
    private void CloseClick(object sender, RoutedEventArgs e) => Close();
    private void SourceClick(object sender, RoutedEventArgs e) => BrowserActions.Open(this, ElkiesSearchFamily.Source, "Elkies family source");
    private void OpenCurveClick(object sender, RoutedEventArgs e)
    { if (Model.CanOpenCurve && Model.Selected is { } candidate) OpenCurveRequested?.Invoke(candidate); }
    private bool CanReplace() => !Model.IsDirty || ConfirmationWindow.Confirm(this, "Replace unsaved search?",
        "The current search checkpoint has not been saved. Save it first if you want to continue this range later.", "Replace search");
    private void NewClick(object sender, RoutedEventArgs e) { if (!Model.IsBusy && CanReplace()) Model.NewSearch(); }
    private void SaveClick(object sender, RoutedEventArgs e)
    {
        var dialog = new SaveFileDialog { Title = "Save Elkies search checkpoint", Filter = "Elkies search (*.ecsearch)|*.ecsearch", DefaultExt = ".ecsearch", FileName = "Elkies-search", AddExtension = true };
        if (dialog.ShowDialog(this) != true) return;
        try { Model.Snapshot.Save(dialog.FileName); Model.MarkSaved(); }
        catch (Exception error) when (error is IOException or UnauthorizedAccessException or ArgumentException or FormatException or OverflowException)
        { Model.ShowError("Could not save search: " + error.Message); }
    }
    private void LoadClick(object sender, RoutedEventArgs e)
    {
        if (Model.IsBusy || !CanReplace()) return;
        var dialog = new OpenFileDialog { Title = "Open Elkies search checkpoint", Filter = "Elkies search (*.ecsearch)|*.ecsearch", CheckFileExists = true };
        if (dialog.ShowDialog(this) != true) return;
        try { Model.Load(CurveSearchState.Load(dialog.FileName)); }
        catch (Exception error) when (error is IOException or UnauthorizedAccessException or ArgumentException or JsonException)
        { Model.ShowError("Could not open search: " + error.Message); }
    }
}

using System.Windows;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.ViewModels;

namespace EllipticCurves.Visualizer;

public partial class CalculationWindow : Window
{
    private readonly CalculationFormViewModel form;
    public event Action<CalculationRequest>? RunRequested;
    public CalculationWindow(CalculationOperation operation, string equation, WorkbenchViewModel workbench, CalculationRequest? previous = null)
    {
        InitializeComponent();
        DataContext = form = new(operation, equation, workbench, previous);
        Closed += (_, _) => form.Dispose();
    }
    private void CloseClick(object sender, RoutedEventArgs e) => Close();
    private void RunClick(object sender, RoutedEventArgs e)
    {
        if (!form.CanRun) return;
        var request = form.CreateRequest();
        Close();
        RunRequested?.Invoke(request);
    }
}

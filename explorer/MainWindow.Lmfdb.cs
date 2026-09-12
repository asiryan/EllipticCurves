using System.Windows.Threading;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private void OpenLmfdbImport()
    {
        if (sessionActionInProgress || !Workbench.CanRun)
        {
            ConfirmationWindow.ShowMessage(this, "Import curve", "Finish the current calculation or session operation before importing a curve.");
            return;
        }
        var dialog = new LmfdbImportWindow { Owner = this };
        dialog.ImportRequested += formula =>
        {
            if (sessionActionInProgress || !Workbench.CanRun)
            {
                ConfirmationWindow.ShowMessage(dialog, "Import curve", "Finish the current calculation or session operation before importing a curve.");
                return;
            }
            ImportLmfdbFormula(formula);
            dialog.Close();
        };
        dialog.Show();
    }

    internal void ImportLmfdbFormula(LmfdbCurveFormula formula)
    {
        if (!CurveEquationText.TryParse(formula.Equation, out var curve, out var error) || curve!.Discriminant.IsZero)
            throw new ArgumentException("Invalid imported formula: " + error, nameof(formula));
        CommitHistory();
        ViewModel.Equation.Text = formula.Equation;
        ViewModel.Equation.CommitEdit();
        ViewModel.FlushUpdate();
        ResetCurveViews(this, EventArgs.Empty);
        // The fit runs after layout, so record formula and viewport as one edit.
        var version = historyRestoreVersion;
        Dispatcher.BeginInvoke(DispatcherPriority.ContextIdle, new Action(() =>
        {
            if (version == historyRestoreVersion) CommitHistory();
        }));
    }
}

using System.Windows.Threading;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.ViewModels;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private CurveSearchWindow? searchWindow;
    private CurveSearchViewModel? searchModel;
    private void OpenCurveSearch()
    {
        if (searchWindow != null) { searchWindow.Activate(); return; }
        if (sessionActionInProgress) return;
        searchModel ??= new(Workbench);
        var dialog = new CurveSearchWindow(searchModel) { Owner = this };
        searchWindow = dialog;
        dialog.Closed += (_, _) => searchWindow = null;
        dialog.OpenCurveRequested += candidate =>
        {
            if (sessionActionInProgress || !Workbench.CanRun) return;
            OpenSearchCurve(candidate);
            dialog.Close();
            Activate();
        };
        dialog.Show();
    }

    internal void OpenSearchCurve(CurveSearchCandidate candidate)
    {
        // Recreate from the family parameter. A saved equation is display data.
        var (curve, _) = ElkiesSearchFamily.Create(candidate.Numerator, candidate.Denominator);
        if (curve.IsSingular) throw new ArgumentException("This parameter gives a singular curve.");
        CommitHistory();
        ViewModel.Equation.Text = CurveEquationText.Format(curve);
        ViewModel.Equation.CommitEdit();
        ViewModel.FlushUpdate();
        ResetCurveViews(this, EventArgs.Empty);
        var version = historyRestoreVersion;
        Dispatcher.BeginInvoke(DispatcherPriority.ContextIdle, new Action(() =>
        { if (version == historyRestoreVersion) CommitHistory(); }));
    }
}

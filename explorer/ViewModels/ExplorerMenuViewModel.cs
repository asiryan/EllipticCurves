#nullable enable
using EllipticCurves.Explorer.Computations;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class ExplorerMenuViewModel : ObservableObject
{
    public static CalculationOperation ImportCurve { get; } = new("explorer.import-curve", "Import curve from LMFDB",
        "LMFDB · internet", "Search LMFDB by conductor or conductor range and import a curve equation.", CalculationContext.Database);
    public static CalculationOperation SearchCurves { get; } = new(CurveSearchEngine.OperationId, "Search the Elkies family",
        "Ranks and arithmetic", "Generate curves from a built-in family with 17 known points. Score candidates, prove rank lower bounds, pause and resume.", CalculationContext.RationalCurve);
    // Keep the catalog's category order, and sort every action within its category.
    // Import opens a picker, so its menu entry is not part of the calculation catalog.
    private static readonly IReadOnlyList<CalculationOperation> MenuOperations = CalculationCatalog.All.Append(ImportCurve).Append(SearchCurves)
        .GroupBy(operation => operation.Group).SelectMany(category => category.OrderBy(operation => operation.Title)).ToArray();
    private string search = "", group = "All calculations";
    public IReadOnlyList<string> Groups { get; } = new[] { "All calculations" }.Concat(CalculationCatalog.All.Select(o => o.Group).Distinct()).ToArray();
    public string Search
    {
        get => search;
        set
        {
            search = value;
            OnPropertyChanged();
            NotifyOperationsChanged();
        }
    }

    public string Group
    {
        get => group;
        set
        {
            group = value;
            OnPropertyChanged();
            NotifyOperationsChanged();
        }
    }

    public IEnumerable<CalculationOperation> Operations => MenuOperations.Where(operation =>
        (Group == "All calculations" || operation.Group == Group) &&
        (operation.Title + " " + operation.Group + " " + operation.Description + " " + operation.Member?.Name).Contains(Search.Trim(), StringComparison.OrdinalIgnoreCase));
    public string Summary => Operations.Count() + " actions shown · select one to continue";

    private void NotifyOperationsChanged()
    {
        OnPropertyChanged(nameof(Operations));
        OnPropertyChanged(nameof(Summary));
    }
}

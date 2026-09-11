#nullable enable
using EllipticCurves.Visualizer.Computations;

namespace EllipticCurves.Visualizer.ViewModels;

public sealed class ExplorerMenuViewModel : ObservableObject
{
    private string search = "", group = "All calculations";
    public IReadOnlyList<string> Groups { get; } = new[] { "All calculations" }.Concat(CalculationCatalog.All.Select(o => o.Group).Distinct()).ToArray();
    public string Search { get => search; set { search = value; OnPropertyChanged(); Update(); } }
    public string Group { get => group; set { group = value; OnPropertyChanged(); Update(); } }
    public IEnumerable<CalculationOperation> Operations => CalculationCatalog.All.Where(operation =>
        (Group == "All calculations" || operation.Group == Group) &&
        (operation.Title + " " + operation.Group + " " + operation.Description + " " + operation.Member?.Name).Contains(Search.Trim(), StringComparison.OrdinalIgnoreCase));
    public string Summary => Operations.Count() + " calculations shown · select one to set its parameters";
    private void Update() { OnPropertyChanged(nameof(Operations)); OnPropertyChanged(nameof(Summary)); }
}

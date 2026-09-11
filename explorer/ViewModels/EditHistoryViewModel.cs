namespace EllipticCurves.Explorer.ViewModels;

public sealed class EditHistoryViewModel : ObservableObject
{
    public bool CanUndo { get; private set; }
    public bool CanRedo { get; private set; }
    internal void Update(bool canUndo, bool canRedo)
    {
        CanUndo = canUndo;
        CanRedo = canRedo;
        OnPropertyChanged(nameof(CanUndo));
        OnPropertyChanged(nameof(CanRedo));
    }
}

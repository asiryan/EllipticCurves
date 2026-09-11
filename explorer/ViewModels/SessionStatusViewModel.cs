#nullable enable
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class SessionStatusViewModel : ObservableObject
{
    private SessionState state;

    public string FileName => SessionFile.GetFileName(state.Path);
    public string DisplayName => FileName + (state.Modified ? " *" : "");
    public string FileLocation => state.Path ?? SessionMessages.NewFileLocation;
    public bool IsSaving => state.Saving;
    public bool NeedsSave => state.NeedsSave;
    public bool CanSave => state.CanSave;
    public bool CanSaveAs => state.CanSaveAs;
    public bool CanExit => state.CanExit;
    public bool CanReplace => state.CanReplace;
    public string Status => state.Saving ? SessionMessages.Saving : state.Failed ? SessionMessages.SaveFailed
        : state.Modified ? SessionMessages.UnsavedChanges : state.Path == null ? SessionMessages.NewSession : SessionMessages.Saved;

    internal void Update(SessionState next)
    {
        if (state == next) return;
        state = next;
        OnPropertyChanged(string.Empty);
    }
}

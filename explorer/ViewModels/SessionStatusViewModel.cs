#nullable enable
using System.IO;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class SessionStatusViewModel : ObservableObject
{
    private (string? Path, bool Modified, bool Saving, bool Failed, bool Busy, bool CanRun) state;

    public string FileName => state.Path == null ? "untitled.ec" : Path.GetFileName(state.Path);
    public string DisplayName => FileName + (state.Modified ? " *" : "");
    public string FileLocation => state.Path ?? "New session: not saved to a file yet.";
    public bool IsSaving => state.Saving;
    public bool NeedsSave => state.Path == null || state.Modified || state.Failed;
    public bool CanSave => !state.Busy && state.Path != null && NeedsSave;
    public bool CanSaveAs => !state.Busy;
    public bool CanExit => !state.Busy;
    public bool CanReplace => !state.Busy && state.CanRun;
    public string Status => state.Saving ? "Saving…" : state.Failed ? "Save failed"
        : state.Modified ? "Unsaved changes" : state.Path == null ? "New session" : "Saved";

    internal void Update(string? path, bool modified, bool saving, bool failed, bool busy, bool canRun)
    {
        var next = (path, modified, saving, failed, busy, canRun);
        if (state == next) return;
        state = next;
        OnPropertyChanged(string.Empty);
    }
}

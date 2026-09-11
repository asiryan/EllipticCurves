#nullable enable
namespace EllipticCurves.Explorer.Models;

// Shared by the menu, keyboard commands and save entry point. This describes
// the live session, independently of the snapshot serialized to disk.
public readonly record struct SessionState(string? Path, bool Modified, bool Saving,
    bool Failed, bool Busy, bool CanRun)
{
    public bool NeedsSave => Path == null || Modified || Failed;
    public bool CanSave => !Busy && Path != null && NeedsSave;
    public bool CanSaveAs => !Busy;
    public bool CanExit => !Busy;
    public bool CanReplace => !Busy && CanRun;
}

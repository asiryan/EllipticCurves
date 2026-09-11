#nullable enable
namespace EllipticCurves.Explorer.Models;

/// <summary>A bounded caretaker. Restoring a memento never executes an operation.</summary>
public sealed class MementoHistory<T>(int capacity = 50)
{
    private readonly int capacity = capacity > 0 ? capacity : throw new ArgumentOutOfRangeException(nameof(capacity));
    private readonly List<T> undo = new(), redo = new();
    public bool CanUndo => undo.Count != 0;
    public bool CanRedo => redo.Count != 0;

    public void Record(T previous)
    {
        undo.Add(previous);
        if (undo.Count > capacity) undo.RemoveAt(0);
        redo.Clear();
    }

    public T Undo(T current) => Move(undo, redo, current);
    public T Redo(T current) => Move(redo, undo, current);

    private static T Move(List<T> from, List<T> to, T current)
    {
        if (from.Count == 0) throw new InvalidOperationException("No history in this direction.");
        var state = from[^1];
        from.RemoveAt(from.Count - 1);
        to.Add(current);
        return state;
    }

    public void Clear() { undo.Clear(); redo.Clear(); }
}

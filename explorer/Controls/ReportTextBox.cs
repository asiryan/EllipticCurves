using System.Windows;
using System.Windows.Controls;

namespace EllipticCurves.Explorer.Controls;

public sealed class ReportTextBox : TextBox
{
    private object? displayedJob;

    protected override void OnPropertyChanged(DependencyPropertyChangedEventArgs e)
    {
        if (e.Property != TextProperty || !IsInitialized || !ReferenceEquals(displayedJob, DataContext))
        {
            base.OnPropertyChanged(e);
            if (e.Property == TextProperty) displayedJob = DataContext;
            return;
        }

        var start = SelectionStart;
        var length = SelectionLength;
        var offset = VerticalOffset;
        var oldText = (string?)e.OldValue ?? "";
        var newText = (string?)e.NewValue ?? "";
        if (newText.Length != oldText.Length)
        {
            // Longer status or timing text shifts the unchanged input/result text.
            var suffix = 0;
            while (suffix < Math.Min(oldText.Length, newText.Length)
                && oldText[^(suffix + 1)] == newText[^(suffix + 1)]) suffix++;
            if (start >= oldText.Length - suffix) start += newText.Length - oldText.Length;
        }

        base.OnPropertyChanged(e);
        start = Math.Clamp(start, 0, Text.Length);
        Select(start, Math.Min(length, Text.Length - start));
        ScrollToVerticalOffset(offset);
    }
}

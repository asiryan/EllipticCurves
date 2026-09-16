using System.Windows;
using System.Windows.Controls;

namespace EllipticCurves.Explorer.Controls;

public sealed class ReportTextBox : TextBox
{
    public static readonly DependencyProperty SourceTextProperty = DependencyProperty.Register(
        nameof(SourceText), typeof(string), typeof(ReportTextBox),
        new PropertyMetadata("", (owner, args) => ((ReportTextBox)owner).Display((string?)args.NewValue ?? "")));

    public string SourceText
    {
        get => (string?)GetValue(SourceTextProperty) ?? "";
        set => SetValue(SourceTextProperty, value);
    }

    private object? displayedJob;
    private string displayedText = "";

    public ReportTextBox() => DataObject.AddCopyingHandler(this, CopySourceText);

    private void Display(string source)
    {
        // WPF can wrap after leading spaces and leave a visually empty line
        // before a long number. Keep the indent with the first word for display.
        // This one-for-one substitution preserves all selection indices.
        char[]? characters = null;
        var leading = true;
        for (var i = 0; i < source.Length; i++)
        {
            if (leading && source[i] == ' ')
            {
                characters ??= source.ToCharArray();
                characters[i] = '\u00a0';
            }
            else leading = source[i] is '\r' or '\n';
        }
        displayedText = characters == null ? source : new string(characters);
        SetCurrentValue(TextProperty, displayedText);
    }

    private void CopySourceText(object sender, DataObjectCopyingEventArgs e)
    {
        if (Text != displayedText) return;
        // Keyboard/context-menu copy must return the same original characters
        // as the report's Copy and Export buttons, including real NBSPs in data.
        var selected = SourceText.Substring(SelectionStart, SelectionLength);
        e.DataObject.SetData(DataFormats.UnicodeText, selected);
        e.DataObject.SetData(DataFormats.Text, selected);
        e.DataObject.SetData(DataFormats.StringFormat, selected);
    }

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

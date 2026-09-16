using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using EllipticCurves.Explorer.Controls;

internal static partial class Program
{
    private static void CheckReportWrapping()
    {
        const string factors = "    [21447358407501679739, 1]\r\n"
            + "    [229195053659697773915482188781052624679234107, 1]";
        var report = new ReportTextBox
        {
            Style = (Style)Application.Current.FindResource("ReportText"),
            IsReadOnly = true,
            AcceptsReturn = true,
            TextWrapping = TextWrapping.Wrap,
            FontFamily = new FontFamily("Cascadia Mono, Consolas"),
            FontSize = 12,
            VerticalContentAlignment = VerticalAlignment.Top,
            SourceText = factors
        };
        foreach (var width in new[] { 200, 266, 400, 640 })
        {
            report.Measure(new Size(width, 400));
            report.Arrange(new Rect(0, 0, width, 400));
            report.UpdateLayout();
            Require(report.LineCount >= 2, "The report wrapping check needs both factor rows.");
            for (var i = 0; i < report.LineCount; i++)
                Require(!string.IsNullOrWhiteSpace(report.GetLineText(i)),
                    $"An indented factor created an empty visual row at width {width}.");
            Require(report.SourceText == factors, "Wrapping changed the source report used for Copy and Export.");
            report.SelectAll();
            RequireReportCopy(report, factors);
        }

        // A selection can start inside an indent and end inside a wrapped number.
        var start = factors.IndexOf('\n') + 3;
        report.Select(start, 35);
        RequireReportCopy(report, factors.Substring(start, 35));

        // Keep intentional blank lines and nonbreaking spaces in the original data.
        const string spaced = "Heading\r\n\r\n  Value: a\u00a0b\r\n";
        report.SourceText = spaced;
        Require(report.Text.Contains("\r\n\r\n"), "The display removed an intentional blank line.");
        report.SelectAll();
        RequireReportCopy(report, spaced);
    }

    private static void RequireReportCopy(ReportTextBox report, string expected)
    {
        var copied = false;
        void CheckCopy(object sender, DataObjectCopyingEventArgs e)
        {
            try
            {
                copied = true;
                foreach (var format in new[] { DataFormats.UnicodeText, DataFormats.Text, DataFormats.StringFormat })
                    Require(Equals(e.DataObject.GetData(format), expected),
                        $"Copying a report selection changed its original characters in {format}.");
            }
            finally { e.CancelCommand(); } // Leave the user's system clipboard untouched.
        }
        DataObject.AddCopyingHandler(report, CheckCopy);
        try
        {
            Require(ApplicationCommands.Copy.CanExecute(null, report), "Copy must be available for a report selection.");
            ApplicationCommands.Copy.Execute(null, report);
            Require(copied, "The WPF Copy command did not prepare the selected report text.");
        }
        finally { DataObject.RemoveCopyingHandler(report, CheckCopy); }
    }
}

using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Automation;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using EllipticCurves.Visualizer.Models;
using EllipticCurves.Visualizer.Windowing;
using EllipticCurves.Visualizer;
using EllipticCurves.Visualizer.Computations;
using EllipticCurves.Visualizer.Controls;
using EllipticCurves.Visualizer.ViewModels;

internal static class Program
{
    [STAThread]
    private static int Main()
    {
        try
        {
            // Load the actual compiled BAML and theme, without Show() or an app event loop.
            var app = new App();
            app.InitializeComponent();
            using var workbench = new WorkbenchViewModel();
            var request = new CalculationRequest("Q.TorsionStructure", "y^2 = x^3 - x", new());
            var selected = new CalculationJobViewModel(request, "Displayed");
            var other = new CalculationJobViewModel(request, "Right-clicked");
            workbench.Jobs.Add(selected);
            workbench.Jobs.Add(other);
            workbench.Selected = selected;
            var panel = new ResultsPanel { DataContext = workbench };
            panel.Measure(new Size(300, 800));
            panel.Arrange(new Rect(0, 0, 300, 800));
            panel.UpdateLayout();
            var history = Descendants(panel).OfType<ComboBox>().Single();
            var style = history.ItemContainerStyle;
            var opening = style.Setters.OfType<EventSetter>().Single(s => s.Event == FrameworkElement.ContextMenuOpeningEvent);
            var onOpening = (ContextMenuEventHandler)opening.Handler;
            var row = new ComboBoxItem { DataContext = other, Content = other, Style = style };
            var menu = row.ContextMenu ?? throw new Exception("History row has no context menu.");
            Require(menu.Items.Count == 1, "Expected one Delete action.");
            var delete = (MenuItem)menu.Items[0];
            Require((string)delete.Header == "Delete", "Delete action is missing.");

            // Invoke the same EventSetter that WPF raises on right-click. Opening twice
            // must not attach duplicate handlers or act on a stale row.
            onOpening(row, null!);
            onOpening(row, null!);
            var deletions = 0;
            workbench.Jobs.CollectionChanged += (_, _) => deletions++;
            delete.RaiseEvent(new RoutedEventArgs(MenuItem.ClickEvent, delete));
            Require(deletions == 1 && !workbench.Jobs.Contains(other), "Delete must remove the context row exactly once.");
            Require(workbench.Selected == selected, "Delete changed an unrelated displayed result.");
            row.DataContext = selected;
            onOpening(row, null!);
            delete.RaiseEvent(new RoutedEventArgs(MenuItem.ClickEvent, delete));
            Require(workbench.Jobs.Count == 0 && workbench.Selected == null, "Last-result deletion did not clear the panel.");
            Require(Descendants(panel).OfType<Button>().Any(b => Equals(b.Content, "Repeat")), "Repeat button is missing.");
            CheckWorkspaceLayout(request);
            CheckDockAnimation();
            CheckPlotRendering();
            Console.WriteLine("PASS: compiled XAML loads; history deletion, Repeat, title indicator, PNG rendering, navigation placement and both full-height sidebars checked. No windows opened.");
            app.Shutdown();
            return 0;
        }
        catch (Exception error) { Console.Error.WriteLine(error); return 1; }
    }

    private static void Require(bool condition, string message) { if (!condition) throw new Exception(message); }
    private static void CheckPlotRendering()
    {
        foreach (var equation in new[] { "y^2=x^3-x", "y^2+xy+y=x^3-x", "y^2=x^3", "y^2=x^3-3*x-2", "y^2=x^3+1e400*x" })
        {
            Require(CurveEquationText.TryParse(equation, out var curve, out var error), error);
            var plot = new CurvePlot { Snapshot = new CurveSnapshot(curve!) };
            plot.Measure(new Size(800, 500));
            plot.Arrange(new Rect(0, 0, 800, 500));
            plot.Fit();
            plot.Zoom(0.8);
            plot.Zoom(1.25);
            plot.UpdateLayout();
            var bitmap = new RenderTargetBitmap(1600, 1000, 192, 192, PixelFormats.Pbgra32);
            bitmap.Render(plot);
            var encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(bitmap));
            using var stream = new System.IO.MemoryStream();
            encoder.Save(stream);
            stream.Position = 0;
            var decoded = BitmapFrame.Create(stream, BitmapCreateOptions.None, BitmapCacheOption.OnLoad);
            Require(decoded.PixelWidth == 1600 && decoded.PixelHeight == 1000, "PNG export lost its dimensions.");
            var pixels = new byte[1600 * 1000 * 4];
            bitmap.CopyPixels(pixels, 1600 * 4, 0);
            // Antialiased strokes can round composited alpha down by one; check the
            // untouched background and ensure no part of the bitmap is left empty.
            Require(pixels[3] == 255 && pixels[^1] == 255 && pixels.Where((_, i) => i % 4 == 3).All(alpha => alpha != 0),
                "Plot export has a transparent background.");
            if (equation == "y^2=x^3-x")
            {
                var hasCurve = Enumerable.Range(0, 1600 * 1000).Any(i => pixels[4 * i] == 207 && pixels[4 * i + 1] == 230 && pixels[4 * i + 2] == 99);
                Require(hasCurve, "PNG export omitted the real locus.");
            }
        }
    }
    private static void CheckDockAnimation()
    {
        var column = new ColumnDefinition { Width = new GridLength(320) };
        var clock = new GridLengthAnimation { From = 320, To = 32, Duration = TimeSpan.FromMilliseconds(200) }.CreateClock();
        column.ApplyAnimationClock(ColumnDefinition.WidthProperty, clock);
        clock.Controller!.Begin();
        clock.Controller.SeekAlignedToLastTick(TimeSpan.FromMilliseconds(100), TimeSeekOrigin.BeginTime);
        Require(column.Width.Value > 32 && column.Width.Value < 320, "Dock animation did not interpolate the column width.");
        clock.Controller.SeekAlignedToLastTick(TimeSpan.FromMilliseconds(200), TimeSeekOrigin.BeginTime);
        Require(Math.Abs(column.Width.Value - 32) < 0.01, "Dock animation did not finish at the tab width.");
        column.ApplyAnimationClock(ColumnDefinition.WidthProperty, null);
    }
    private static void CheckWorkspaceLayout(CalculationRequest request)
    {
        var window = new MainWindow();
        try
        {
            var root = (FrameworkElement)window.Content;
            var size = new Size(1438, 918);
            void Layout()
            {
                root.Measure(size);
                root.Arrange(new Rect(new Point(), size));
                root.UpdateLayout();
                CheckPlotControls();
                CheckSidebarDividers();
            }
            Rect Bounds(FrameworkElement element) => element.TransformToAncestor(root).TransformBounds(new Rect(element.RenderSize));
            void CheckSidebarDividers()
            {
                var left = (FrameworkElement)window.FindName("EquationPanel");
                var right = (FrameworkElement)window.FindName("Results");
                var leftDivider = (GridSplitter)window.FindName("EquationSplitter");
                var rightDivider = (GridSplitter)window.FindName("ResultsSplitter");
                Require(leftDivider.Visibility == left.Visibility && rightDivider.Visibility == right.Visibility,
                    "A sidebar divider must fold and reopen with its panel.");
                if (left.Visibility == Visibility.Collapsed) left = (FrameworkElement)window.FindName("EquationTab");
                if (right.Visibility == Visibility.Collapsed) right = (FrameworkElement)window.FindName("ResultsTab");
                var plot = Bounds((FrameworkElement)window.FindName("PlotCard"));
                var leftGap = plot.Left - Bounds(left).Right;
                var rightGap = Bounds(right).Left - plot.Right;
                Require(Math.Abs(leftGap - rightGap) < 1, "Left and right sidebar gaps differ.");
                if (leftDivider.Visibility == Visibility.Visible && rightDivider.Visibility == Visibility.Visible)
                    Require(leftDivider.RenderSize == rightDivider.RenderSize,
                        "The sidebar dividers must have the same width and height.");
            }
            void CheckPlotControls()
            {
                var plotBounds = Bounds((FrameworkElement)window.FindName("Plot"));
                var cardBounds = Bounds((FrameworkElement)window.FindName("PlotCard"));
                var navigation = (Panel)window.FindName("PlotNavigation");
                var legend = (Panel)window.FindName("PlotLegend");
                foreach (var button in navigation.Children.OfType<Button>())
                {
                    var bounds = Bounds(button);
                    Require(bounds.Top >= plotBounds.Bottom, "Plot navigation overlaps the graph or axis labels.");
                    Require(cardBounds.Contains(bounds), "Plot navigation extends outside the plot card.");
                    foreach (var item in legend.Children.OfType<FrameworkElement>())
                        Require(!Bounds(item).IntersectsWith(bounds), "Plot navigation overlaps the legend.");
                }
            }
            void CheckRail(Rect expanded, FrameworkElement rail, bool right)
            {
                var folded = Bounds(rail);
                Require(Math.Abs(folded.Width - 32) < 1, "A folded sidebar must be 32 pixels wide.");
                Require(Math.Abs(folded.Top - expanded.Top) < 1 && Math.Abs(folded.Bottom - expanded.Bottom) < 1,
                    "Folding changed the sidebar's top, bottom or height.");
                Require(Math.Abs(right ? folded.Right - expanded.Right : folded.Left - expanded.Left) < 1,
                    "Folding moved the sidebar's outer edge.");
            }
            void Click(Button button)
            {
                button.RaiseEvent(new RoutedEventArgs(Button.ClickEvent, button));
                Layout();
            }
            Layout();
            var results = (ResultsPanel)window.FindName("Results");
            var tab = (Button)window.FindName("ResultsTab");
            var column = (ColumnDefinition)window.FindName("ResultsColumn");
            var plot = (CurvePlot)window.FindName("Plot");
            var plotCard = (Border)window.FindName("PlotCard");
            var export = (Button)window.FindName("ExportPlotButton");
            var equation = (Border)window.FindName("EquationPanel");
            var equationColumn = (ColumnDefinition)window.FindName("EquationColumn");
            var equationTab = (Button)window.FindName("EquationTab");
            var equationCollapse = (Button)window.FindName("EquationCollapseButton");
            var coefficients = Descendants(equation).OfType<Expander>().Single();
            var equationScroll = (ScrollViewer)window.FindName("EquationScroll");
            var equationIntro = (TextBlock)window.FindName("EquationIntro");
            Require(equationScroll.ScrollableHeight == 0, "The collapsed controls should fit at the default window size.");
            var textBounds = Bounds(equationIntro);
            var viewportWidth = equationScroll.ViewportWidth;
            equationScroll.Height = 280;
            Layout();
            Require(equationScroll.ScrollableHeight > 0, "The overflow layout did not enable scrolling.");
            Require(Math.Abs(equationScroll.ViewportWidth - viewportWidth) < 1 && equationIntro.RenderSize == textBounds.Size,
                "The sidebar scrollbar narrowed the text and changed its wrapping.");
            var scrollBar = Descendants(equationScroll).OfType<System.Windows.Controls.Primitives.ScrollBar>().Single();
            Require(Bounds(scrollBar).Left >= textBounds.Right,
                $"The sidebar scrollbar overlaps the text column: bar {Bounds(scrollBar)}, text {textBounds}, padding {equationScroll.Padding}.");
            equationScroll.ScrollToBottom();
            Layout();
            Require(equationScroll.VerticalOffset > 0, "The sidebar content no longer scrolls.");
            equationScroll.Height = double.NaN;
            equationScroll.ScrollToTop();
            Layout();
            coefficients.IsExpanded = true;
            Layout();
            Require(equationIntro.RenderSize == textBounds.Size, "Expanding Coefficients changed the help text wrapping.");
            Require(Descendants(plotCard).Contains(export), "Export must be on the plot panel.");
            Require(window.FindName("NativeIndicator") != null, "Native status indicator is missing.");
            Require(window.FindName("ResultsToggle") == null, "The old Results toolbar button remains.");
            Require(!Descendants(root).OfType<TextBlock>().Any(t => t.Text.Contains("A little change")), "The slogan remains.");
            Require(Math.Abs(results.ActualWidth - 300) < 1, "Results must start compact.");
            Require(Math.Abs(equation.ActualWidth - 238) < 1, "Equation must start compact.");

            var saved = new CalculationJobViewModel(request, "Preserved while folded");
            window.Workbench.Jobs.Add(saved);
            window.Workbench.Selected = saved;
            column.Width = new GridLength(470);
            equationColumn.Width = new GridLength(320);
            Layout();
            var expandedResults = Bounds(results);
            var expandedEquation = Bounds(equation);
            var expandedPlotWidth = plot.ActualWidth;
            var collapse = Descendants(results).OfType<Button>().Single(b => AutomationProperties.GetName(b) == "Collapse results panel");
            Click(collapse);
            Require(results.Visibility == Visibility.Collapsed && tab.Visibility == Visibility.Visible, "Results did not fold into the edge tab.");
            Require(column.ActualWidth <= 40 && plot.ActualWidth > expandedPlotWidth, "Folding must give space back to the plot.");
            CheckRail(expandedResults, tab, right: true);
            var oneFoldedPlotWidth = plot.ActualWidth;
            Click(equationCollapse);
            Require(equation.Visibility == Visibility.Collapsed && equationTab.Visibility == Visibility.Visible, "Equation did not fold into its side tab.");
            CheckRail(expandedEquation, equationTab, right: false);
            CheckRail(expandedResults, tab, right: true);
            Require(plot.ActualWidth > oneFoldedPlotWidth, "Folding Equation must give more space to the plot.");
            Click(tab);
            Require(results.Visibility == Visibility.Visible && tab.Visibility == Visibility.Collapsed, "The side tab did not reopen results.");
            Require(Math.Abs(column.ActualWidth - 470) < 1, "The resized results width was not retained.");
            Click(equationTab);
            Require(equation.Visibility == Visibility.Visible && equationTab.Visibility == Visibility.Collapsed, "The side tab did not reopen Equation.");
            Require(Math.Abs(equation.ActualWidth - 320) < 1 && coefficients.IsExpanded, "Folding lost the resized Equation width or Coefficients state.");
            Require(window.Workbench.Selected == saved, "Folding lost the selected result.");

            column.Width = new GridLength(300);
            equationColumn.Width = new GridLength(238);
            foreach (var windowSize in new[] { new Size(1118, 758), new Size(1438, 918), new Size(1918, 1078) })
            {
                size = windowSize;
                Layout();
                expandedResults = Bounds(results);
                expandedEquation = Bounds(equation);
                Require(Math.Abs(expandedResults.Top - expandedEquation.Top) < 1 && Math.Abs(expandedResults.Bottom - expandedEquation.Bottom) < 1,
                    "The two sidebars must share the same vertical bounds.");
                Click(collapse);
                CheckRail(expandedResults, tab, right: true);
                Click(equationCollapse);
                CheckRail(expandedResults, tab, right: true);
                CheckRail(expandedEquation, equationTab, right: false);
                Click(tab);
                CheckRail(expandedEquation, equationTab, right: false);
                Click(equationTab);
                Require(Bounds(results) == expandedResults && Bounds(equation) == expandedEquation,
                    "Reopening both sidebars did not restore their original bounds.");
            }
            column.Width = new GridLength(column.MaxWidth);
            Layout();
            Require(Math.Abs(plotCard.ActualWidth - 400) < 1, "The narrowest plot layout was not checked.");
            equationColumn.Width = new GridLength(450);
            column.Width = new GridLength(900);
            Layout();
            size = new Size(1118, 758);
            Layout();
            var workspace = (FrameworkElement)window.FindName("Workspace");
            Require(plotCard.ActualWidth >= 400 && Bounds(workspace).Contains(Bounds(results)) && Bounds(workspace).Contains(Bounds(equation)),
                "Shrinking the window after widening both sidebars pushed content outside the workspace.");
            size = new Size(1918, 1078);
            Layout();
            Require(Math.Abs(equation.ActualWidth - 450) < 1 && Math.Abs(results.ActualWidth - 900) < 1,
                "Growing the window did not restore the requested sidebar widths.");
        }
        finally { window.Close(); window.Workbench.Dispose(); window.ViewModel.Dispose(); }
    }
    private static IEnumerable<DependencyObject> Descendants(DependencyObject root)
    {
        for (var i = 0; i < VisualTreeHelper.GetChildrenCount(root); i++)
        {
            var child = VisualTreeHelper.GetChild(root, i);
            yield return child;
            foreach (var descendant in Descendants(child)) yield return descendant;
        }
    }
}

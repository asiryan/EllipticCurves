using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Automation;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using System.Windows.Threading;
using System.Windows.Interop;
using System.Runtime.InteropServices;
using EllipticCurves;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.Windowing;
using EllipticCurves.Explorer;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Controls;
using EllipticCurves.Explorer.ViewModels;

internal static partial class Program
{
    [STAThread]
    private static int Main()
    {
        try
        {
            // Load the actual compiled controls and theme. App.OnStartup calls
            // MainWindow.Show(), so use a resource-only Application for tests
            // that pump the dispatcher, including its queued Startup event.
            var app = new Application { ShutdownMode = ShutdownMode.OnExplicitShutdown };
            app.Resources.MergedDictionaries.Add(new ResourceDictionary
                { Source = new Uri("/EllipticCurves.Explorer;component/Themes/Theme.xaml", UriKind.Relative) });
            CheckTitleBarToggle(new ExplorerMenu());
            CheckTitleBarToggle(new SessionMenu());
            CheckTitleBarToggle(new EditMenu());
            CheckTitleBarToggle(new HelpMenu());
            CheckCaptionDismissal("Explorer");
            CheckCaptionDismissal("Session");
            CheckCaptionDismissal("Edit");
            CheckCaptionDismissal("Help");
            CheckHelp();
            CheckEditHistory();
            CheckSessionMenu();
            CheckSessionRestore();
            CheckSessionLifecycle();
            CheckSessionNavigation();
            CheckSessionSaveName();
            CheckSessionShortcuts();
            CheckSessionSaveStatus();
            CheckSessionSaveAsOverwrite();
            CheckSessionSavingConcurrency();
            CheckSaveChangesDialog();
            CheckExplorerSelection();
            CheckLmfdbImport();
            Require(app.MainWindow == null, "The presentation host must not launch the application window.");
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
            Require(Descendants(panel).OfType<TextBlock>().Any(text => text.Text == "SESSION HISTORY · LAST 50"),
                "The history heading must show the session's report limit.");
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
            CheckConfirmations(request);
            CheckConfirmationDialogs();
            CheckDockAnimation();
            CheckExportBounds();
            CheckPlotRendering();
            CheckExtremePlotViews();
            CheckLargeCurvePlotScales();
            CheckViewSwitching();
            CheckTorusCycleColors();
            CheckComplexTorusView();
            Console.WriteLine("PASS: compiled XAML loads; Help actions, bundled license and layouts, LMFDB formula import, Edit Undo/Redo, graph and result mementos, File save/load and shared title-bar menus, Tools click/focus scrolling, history deletion, themed Clear/Reset dialogs and confirmation paths, Repeat, PNG rendering, navigation placement and both full-height sidebars checked. No windows shown.");
            app.Shutdown();
            return 0;
        }
        catch (Exception error) { Console.Error.WriteLine(error); return 1; }
    }

    private static void Require(bool condition, string message) { if (!condition) throw new Exception(message); }

    private static MainWindow CreateMainWindow(Func<bool>? confirmation = null) => new(confirmation,
        new SessionDialogs(_ => new(SaveChangesChoice.Discard), _ => throw new Exception("Unexpected save dialog."),
            () => throw new Exception("Unexpected open dialog."), (title, message) => throw new Exception(title + ": " + message)));

    private static void CheckViewSwitching()
    {
        var window = CreateMainWindow();
        try
        {
            window.ViewModel.ShowPoints = false;
            var root = (FrameworkElement)window.Content;
            var mode = (ComboBox)window.FindName("ViewMode");
            var plot = (CurvePlot)window.FindName("Plot");
            void Layout()
            {
                root.Measure(new Size(1438, 918));
                root.Arrange(new Rect(0, 0, 1438, 918));
                root.UpdateLayout();
                Dispatcher.CurrentDispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(root.UpdateLayout));
            }
            (Point, Point) Projection() => (plot.ToScreen(0, 0), plot.ToScreen(1, 1));
            void RequireFitted()
            {
                var viewport = new Rect(plot.RenderSize);
                foreach (var x in plot.Snapshot!.Plot.Roots)
                    Require(viewport.Contains(plot.ToScreen(x, plot.Snapshot.Plot.CenterY(x))),
                        "A preset selected in Complex torus left real branch points outside the viewport.");
                var before = Projection();
                plot.Fit();
                Require(Projection() == before, "The deferred real view does not match Reset view at its actual size.");
            }

            // Select a curve before the real plot has ever been laid out.
            mode.SelectedIndex = 1;
            window.ViewModel.ApplyPreset(CurvePreset.All.Single(preset => preset.Name == "48.a3"));
            Layout();
            mode.SelectedIndex = 0;
            Layout();
            RequireFitted();

            plot.Zoom(0.4);
            var zoomed = Projection();
            mode.SelectedIndex = 1;
            Layout();
            window.ViewModel.FitCommand.Execute(null); // Ctrl+F resets only the torus camera.
            Layout();
            mode.SelectedIndex = 0;
            Layout();
            Require(Projection() == zoomed, "Switching modes or resetting the torus discarded the real viewport.");

            mode.SelectedIndex = 1;
            window.ViewModel.Equation.Text = "y^2 = x^3 - 20*x + 10";
            window.ViewModel.FlushUpdate();
            Layout();
            mode.SelectedIndex = 0;
            Layout();
            Require(Projection() == zoomed, "Editing coefficients in Complex torus reset the real viewport.");

            // Cover a return both before and after the queued reset has run while hidden.
            foreach (var flushWhileHidden in new[] { false, true })
            {
                mode.SelectedIndex = 1;
                window.ViewModel.ResetCommand.Execute(null);
                if (flushWhileHidden) Layout();
                mode.SelectedIndex = 0;
                Layout();
                RequireFitted();
                plot.Zoom(0.2);
            }
        }
        finally { window.Close(); }
    }

    private static void CheckTorusCycleColors()
    {
        var torus = new TorusViewport { ShowGrid = false };
        torus.Zoom(0.75);
        var host = new Border { Background = new SolidColorBrush(Color.FromRgb(0x12, 0x19, 0x20)), Child = torus };
        host.Measure(new Size(600, 440));
        host.Arrange(new Rect(0, 0, 600, 440));
        host.UpdateLayout();
        var bitmap = PlotImageExporter.Render(host);
        var pixels = new byte[bitmap.PixelWidth * bitmap.PixelHeight * 4];
        bitmap.CopyPixels(pixels, bitmap.PixelWidth * 4, 0);
        foreach (var color in new[] { Color.FromRgb(0x63, 0xE6, 0xCF), Color.FromRgb(0x84, 0xAA, 0xFF) })
        {
            var matches = 0;
            for (var i = 0; i < pixels.Length; i += 4)
                if (Math.Abs(pixels[i] - color.B) <= 12 && Math.Abs(pixels[i + 1] - color.G) <= 12
                    && Math.Abs(pixels[i + 2] - color.R) <= 12) matches++;
            Require(matches > 200, "A torus cycle lost its legend color against the surface: " + color);
        }
    }

    private static void CheckComplexTorusView()
    {
        using var view = new ComplexTorusView();
        var host = new Border { Padding = new Thickness(8, 59, 8, 0), Child = view };
        var curve = new EllipticCurveQ(0, 0, 0, -1, 0);
        var dispatcher = Dispatcher.CurrentDispatcher;
        dispatcher.Invoke(() => view.Model.Update(curve, new[]
        {
            new EllipticCurvePoint(-1, 0), new EllipticCurvePoint(0, 0), new EllipticCurvePoint(1, 0)
        }, true));
        var preparation = view.Model.PendingUpdate;
        var frame = new DispatcherFrame();
        var watch = System.Diagnostics.Stopwatch.StartNew();
        var timer = new DispatcherTimer(TimeSpan.FromMilliseconds(20), DispatcherPriority.Background,
            (_, _) => { if (preparation.IsCompleted || watch.Elapsed.TotalSeconds > 20) frame.Continue = false; }, dispatcher);
        try { Dispatcher.PushFrame(frame); }
        finally { timer.Stop(); }
        Require(preparation.IsCompleted, "Complex torus preparation did not finish.");
        preparation.GetAwaiter().GetResult();
        Require(view.Model.HasLattice && view.Model.Points.Count == 4, "Complex view lost the classic curve's half-period points.");

        void Layout(double width, double height)
        {
            host.Measure(new Size(width + 16, height + 59));
            host.Arrange(new Rect(0, 0, width + 16, height + 59));
            host.UpdateLayout();
            dispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(host.UpdateLayout));
        }
        Layout(800, 540);
        var lattice = (PeriodLatticePlot)view.FindName("LatticePlot");
        var torus = (TorusViewport)view.FindName("TorusPlot");
        var picker = Descendants(view).OfType<ComboBox>().Single();
        picker.SelectedIndex = 2;
        dispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(view.UpdateLayout));
        Require(Equals(lattice.SelectedPoint, view.Model.Points[2]) && Equals(torus.SelectedPoint, lattice.SelectedPoint),
            "The point selector must highlight the same point on the lattice and torus.");
        torus.SetCurrentValue(TorusViewport.SelectedPointProperty, view.Model.Points[3]);
        dispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(view.UpdateLayout));
        Require(Equals(picker.SelectedItem, view.Model.Points[3]) && Equals(lattice.SelectedPoint, picker.SelectedItem),
            "Selecting a torus marker did not update the other views.");
        picker.SelectedIndex = 0;
        dispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(view.UpdateLayout));
        view.Model.RestoreSelection(view.Model.Points[2].Point.ToString());
        using (var source = new HwndSource(new HwndSourceParameters("Session selection check") { ParentWindow = new IntPtr(-3) }))
            lattice.RaiseEvent(new KeyEventArgs(Keyboard.PrimaryDevice, source, Environment.TickCount, Key.Home)
                { RoutedEvent = Keyboard.KeyDownEvent });
        Require(view.Model.SessionSelection == EllipticCurvePoint.Infinity.ToString(),
            "Home on an already selected origin must cancel the deferred saved selection through the WPF binding.");
        view.ShowGrid = false;
        Require(!lattice.ShowGrid && !torus.ShowGrid, "Complex grid visibility is not shared.");
        view.Zoom(0.8);
        view.Fit();
        var bitmap = PlotImageExporter.Render(view);
        var pixels = new byte[bitmap.PixelWidth * bitmap.PixelHeight * 4];
        bitmap.CopyPixels(pixels, bitmap.PixelWidth * 4, 0);
        Require(bitmap.PixelWidth == 1600 && bitmap.PixelHeight == 1080
            && pixels.Where((_, i) => i % 4 == 3).All(alpha => alpha != 0),
            "The embedded complex view exported with an offset or transparent background.");
        Layout(400, 240);
        Require(lattice.ActualWidth >= 100 && lattice.ActualHeight >= 90 && torus.ActualHeight >= 90,
            "The complex diagrams collapsed in a small viewport.");
        var scroll = Descendants(view).OfType<ScrollViewer>().First();
        Require(scroll.ScrollableHeight > 0, "The compact complex view must scroll instead of crushing the diagrams.");
    }

    private static void CheckSessionMenu()
    {
        var owner = CreateMainWindow();
        try
        {
            var session = (SessionMenu)owner.FindName("Session");
            var explorer = (ExplorerMenu)owner.FindName("Explorer");
            var edit = (EditMenu)owner.FindName("Edit");
            var help = (HelpMenu)owner.FindName("Help");
            var header = (Panel)session.Parent;
            Require(header.Children.IndexOf(session) + 1 == header.Children.IndexOf(edit)
                && header.Children.IndexOf(edit) + 1 == header.Children.IndexOf(explorer)
                && header.Children.IndexOf(explorer) + 1 == header.Children.IndexOf(help),
                "Title-bar menus must appear in File, Edit, Tools, Help order.");
            var toggle = (ToggleButton)session.FindName("Toggle");
            var otherToggle = (ToggleButton)explorer.FindName("Toggle");
            Require(ReferenceEquals(toggle.Style, otherToggle.Style), "Title-bar menu styles differ.");
            foreach (var control in new UserControl[] { session, explorer })
                System.Windows.Data.BindingOperations.ClearBinding((Popup)control.FindName("MenuPopup"), Popup.IsOpenProperty);
            toggle.IsChecked = true;
            otherToggle.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = Mouse.PreviewMouseDownEvent });
            Require(toggle.IsChecked == false, "Opening Explorer did not dismiss Session.");
            otherToggle.IsChecked = true;
            toggle.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = Mouse.PreviewMouseDownEvent });
            Require(otherToggle.IsChecked == false, "Opening Session did not dismiss Explorer.");

            // Use an unconnected menu so actions can be exercised without file dialogs.
            var menu = new SessionMenu();
            var popup = (Popup)menu.FindName("MenuPopup");
            var menuToggle = (ToggleButton)menu.FindName("Toggle");
            System.Windows.Data.BindingOperations.ClearBinding(popup, Popup.IsOpenProperty);
            var actions = new List<string>();
            menu.NewRequested += () => actions.Add("New");
            menu.OpenRequested += () => actions.Add("Open");
            menu.SaveRequested += () => actions.Add("Save");
            menu.SaveAsRequested += () => actions.Add("Save as");
            menu.ExitRequested += () => actions.Add("Exit");
            var buttons = Descendants(popup.Child).OfType<Button>().ToArray();
            Require(buttons.Select(button => button.Content).SequenceEqual(new[] { "New", "Open", "Save", "Save as", "Exit" }),
                "Session must contain New, Open, Save, Save as and Exit, in English.");
            var menuRoot = (FrameworkElement)popup.Child;
            menuRoot.Measure(new Size(220, double.PositiveInfinity));
            menuRoot.Arrange(new Rect(menuRoot.DesiredSize));
            menuRoot.UpdateLayout();
            var shortcutRights = new List<double>();
            foreach (var (button, shortcut) in buttons.Take(4).Zip(new[] { "Ctrl+N", "Ctrl+O", "Ctrl+S", "Ctrl+Shift+S" }))
            {
                var texts = Descendants(button).OfType<TextBlock>().ToArray();
                var label = texts.Single(text => text.Text == (string)button.Content);
                var hint = texts.Single(text => text.Text == shortcut);
                Rect Bounds(FrameworkElement element) => element.TransformToAncestor(menuRoot).TransformBounds(new Rect(element.RenderSize));
                Require(Bounds(hint).Left > Bounds(label).Right && new Rect(menuRoot.RenderSize).Contains(Bounds(hint)),
                    "A session shortcut overlaps its label or extends outside the menu.");
                Require(AutomationProperties.GetAcceleratorKey(button) == shortcut, "The session shortcut is missing from accessibility information.");
                shortcutRights.Add(Bounds(hint).Right);
            }
            Require(shortcutRights.Max() - shortcutRights.Min() < 1, "Session shortcuts must align at the right edge.");
            foreach (var button in buttons)
            {
                menuToggle.IsChecked = true;
                button.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                Require(menuToggle.IsChecked == false, "A session action left its popup open.");
            }
            Require(actions.SequenceEqual(new[] { "New", "Open", "Save", "Save as", "Exit" }), "Session actions were not dispatched exactly once.");

            if (Environment.GetEnvironmentVariable("EC_SESSION_PREVIEW") is { Length: > 0 } preview)
            {
                var root = (FrameworkElement)owner.Content;
                root.Measure(new Size(1438, 918));
                root.Arrange(new Rect(0, 0, 1438, 918));
                root.UpdateLayout();
                var sessionPopup = (Popup)session.FindName("MenuPopup");
                System.Windows.Data.BindingOperations.ClearBinding(sessionPopup, Popup.IsOpenProperty);
                toggle.IsChecked = true;
                var dropdown = (FrameworkElement)sessionPopup.Child;
                dropdown.Measure(new Size(220, double.PositiveInfinity));
                dropdown.Arrange(new Rect(dropdown.DesiredSize));
                dropdown.UpdateLayout();
                var drawing = new DrawingVisual();
                using (var context = drawing.RenderOpen())
                {
                    context.DrawRectangle(new VisualBrush(root), null, new Rect(0, 0, 1438, 918));
                    var location = session.TranslatePoint(new Point(0, session.ActualHeight), root);
                    context.DrawRectangle(new VisualBrush(dropdown), null, new Rect(location, dropdown.RenderSize));
                }
                var bitmap = new RenderTargetBitmap(1438, 918, 96, 96, PixelFormats.Pbgra32);
                bitmap.Render(drawing);
                var encoder = new PngBitmapEncoder();
                encoder.Frames.Add(BitmapFrame.Create(bitmap));
                using var stream = System.IO.File.Create(preview);
                encoder.Save(stream);
                toggle.IsChecked = false;
            }
        }
        finally { owner.Close(); }
    }

    private static void CheckSessionRestore()
    {
        var original = CreateMainWindow();
        var restored = CreateMainWindow();
        var path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), Guid.NewGuid() + ".ec");
        try
        {
            void Layout(Window window)
            {
                var root = (FrameworkElement)window.Content;
                root.Measure(new Size(1438, 918));
                root.Arrange(new Rect(0, 0, 1438, 918));
                root.UpdateLayout();
            }
            Layout(original);
            original.ViewModel.Equation.Text = "y^2 + x*y + y = x^3 - 1/3*x + 2/7";
            original.ViewModel.Step.Text = "1/7";
            original.ViewModel.FlushUpdate();
            original.ViewModel.ActiveCoefficients[3].SliderOffset = 4;
            original.ViewModel.ShowGrid = false;
            original.ViewModel.ShowPoints = false;
            ((Expander)original.FindName("CoefficientsExpander")).IsExpanded = true;
            ((ColumnDefinition)original.FindName("EquationColumn")).Width = new GridLength(320);
            ((ColumnDefinition)original.FindName("ResultsColumn")).Width = new GridLength(410);
            Layout(original);
            ((Button)original.FindName("EquationCollapseButton")).RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            ((CurvePlot)original.FindName("Plot")).RestoreView(new PlotViewState(2.3, -4.5, 7.8));
            ((ComplexTorusView)original.FindName("TorusView")).RestoreCamera(new TorusCameraState(-71, 23, 8.9));
            ((ComboBox)original.FindName("ViewMode")).SelectedIndex = 1;
            original.Workbench.Jobs.Add(CalculationJobViewModel.FromSession(new CalculationSession(
                new CalculationRequest("Q.TorsionStructure", "y^2 = x^3 - x", new(), 45, 99),
                DateTime.Now, "Completed", "Done", TimeSpan.FromSeconds(3), 100, "Saved exact result")));
            original.Workbench.Selected = original.Workbench.Jobs[0];
            var saved = original.CaptureSession();
            SessionFile.Save(path, saved);
            restored.RestoreSession(SessionFile.Load(path));
            Layout(restored);
            Dispatcher.CurrentDispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(() => Layout(restored)));
            var actual = restored.CaptureSession();
            Require(actual.Equation == saved.Equation && restored.ViewModel.Step.Text == ExplorerSession.DefaultSliderStep
                && restored.ViewModel.ActiveCoefficients.All(coefficient => coefficient.SliderOffset == 0)
                && original.ViewModel.Step.Text == "1/7", "Opening must recover the equation with default sliders; saving must leave local settings unchanged.");
            Require(((ComboBox)restored.FindName("ViewMode")).SelectedIndex == 0
                && ((ComplexTorusView)restored.FindName("TorusView")).CaptureCamera() == TorusCameraState.Default,
                "The session must reopen in Real locus with default cameras.");
            Require(restored.ViewModel.ShowGrid && restored.ViewModel.ShowPoints
                && !((Expander)restored.FindName("CoefficientsExpander")).IsExpanded,
                "The document must not restore the previous window's visual settings.");
            Require(restored.Workbench.Selected?.Report == original.Workbench.Selected.Report,
                "The session lost calculation results or repeat parameters.");
            var before = System.Text.Json.JsonSerializer.Serialize(actual);
            saved.History.Add(null!);
            try { restored.RestoreSession(saved); throw new Exception("A malformed session was accepted."); }
            catch (System.IO.InvalidDataException) { }
            Require(System.Text.Json.JsonSerializer.Serialize(restored.CaptureSession()) == before,
                "A malformed session partially replaced the current workspace.");
            restored.ViewModel.Equation.Text = "y^2 =";
            try { restored.CaptureSession(); throw new Exception("An incomplete equation was silently saved."); }
            catch (InvalidOperationException) { }
        }
        finally { original.Close(); restored.Close(); System.IO.File.Delete(path); }
    }

    private static void CheckCaptionDismissal(string menuName)
    {
        var owner = CreateMainWindow();
        var menu = (UserControl)owner.FindName(menuName);
        var toggle = (ToggleButton)menu.FindName("Toggle");
        var popup = (Popup)menu.FindName("MenuPopup");
        // Keep both windows hidden, but send real native messages through the
        // same HWND and WindowChrome hooks used by the application.
        System.Windows.Data.BindingOperations.ClearBinding(popup, Popup.IsOpenProperty);
        try
        {
            var handle = new WindowInteropHelper(owner).EnsureHandle();
            var source = HwndSource.FromHwnd(handle);
            var expectedMessage = 0;
            var forwarded = false;
            var openWhenForwarded = false;
            IntPtr ObserveMessage(IntPtr hwnd, int message, IntPtr wParam, IntPtr lParam, ref bool handled)
            {
                if (message != expectedMessage) return IntPtr.Zero;
                forwarded = true;
                openWhenForwarded = toggle.IsChecked == true;
                // Stop only in the test, before Windows starts a modal move or
                // resize loop. The menu's hook must pass the message through.
                handled = true;
                return IntPtr.Zero;
            }
            source.AddHook(ObserveMessage);
            try
            {
                void Send(int message)
                {
                    expectedMessage = message;
                    forwarded = false;
                    // HTCAPTION, plus XBUTTON1 for the extended-button message.
                    var hitTest = message == 0x00AB ? 0x00010002 : 2;
                    SendMessage(handle, message, new IntPtr(hitTest), IntPtr.Zero);
                    Require(forwarded, "Explorer swallowed a native message needed by WindowChrome.");
                }

                foreach (var message in new[] { 0x00A1, 0x00A4, 0x00A7, 0x00AB })
                {
                    toggle.IsChecked = true;
                    Send(0x00A0); // WM_NCMOUSEMOVE
                    Require(openWhenForwarded, "Hovering over the title bar dismissed Explorer.");
                    Send(message); // WM_NC[L/R/M/X]BUTTONDOWN
                    Require(!openWhenForwarded && toggle.IsChecked == false,
                        $"Title-bar click 0x{message:X} did not dismiss Explorer before native window handling.");
                }
                Require(!owner.IsVisible && !popup.IsOpen, "Native input checks must not show windows.");
            }
            finally { source.RemoveHook(ObserveMessage); }
        }
        finally { owner.Close(); }
    }

    [DllImport("user32.dll", EntryPoint = "SendMessageW")]
    private static extern IntPtr SendMessage(IntPtr window, int message, IntPtr wParam, IntPtr lParam);

    private static void CheckTitleBarToggle(UserControl menu)
    {
        var toggle = (ToggleButton)menu.FindName("Toggle");
        var popup = (Popup)menu.FindName("MenuPopup");
        Require(popup.StaysOpen, "Automatic popup dismissal must not preempt the Explorer toggle's closing click.");
        // Exercise routed input and the real ToggleButton click without showing
        // either the owner window or its popup.
        System.Windows.Data.BindingOperations.ClearBinding(popup, Popup.IsOpenProperty);
        var outside = new Button { Content = "Outside Explorer" };
        var panel = new StackPanel { Orientation = Orientation.Horizontal, Children = { menu, outside } };
        var owner = new Window { Content = panel };
        const System.Reflection.BindingFlags methods = System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic;
        var click = typeof(ToggleButton).GetMethod("OnClick", methods)!;
        MouseButtonEventArgs MouseEvent(UIElement target, RoutedEvent routedEvent)
        {
            var input = new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = routedEvent };
            target.RaiseEvent(input);
            return input;
        }
        try
        {
            using var source = new HwndSource(new HwndSourceParameters("Explorer toggle headless layout")
                { WindowStyle = 0, Width = 300, Height = 44 }) { RootVisual = panel };
            panel.Measure(new Size(300, 44));
            panel.Arrange(new Rect(0, 0, 300, 44));
            panel.UpdateLayout();
            var chevron = (FrameworkElement)toggle.Template.FindName("Chevron", toggle);
            var points = new[]
            {
                toggle.TranslatePoint(new Point(3, toggle.ActualHeight / 2), menu),
                toggle.TranslatePoint(new Point(toggle.ActualWidth / 2, toggle.ActualHeight / 2), menu),
                chevron.TranslatePoint(new Point(chevron.ActualWidth / 2, chevron.ActualHeight / 2), menu)
            };
            foreach (var point in points)
            {
                var target = (UIElement)menu.InputHitTest(point);
                for (var i = 0; i < 4; i++)
                {
                    var wasOpen = toggle.IsChecked == true;
                    MouseEvent(target, Mouse.PreviewMouseDownEvent);
                    Require((toggle.IsChecked == true) == wasOpen, "Explorer closed before its header click completed.");
                    MouseEvent(target, Mouse.PreviewMouseUpEvent);
                    click.Invoke(toggle, null);
                    Require((toggle.IsChecked == true) != wasOpen, "Repeated Explorer clicks must alternate open and closed.");
                }
            }
            toggle.IsChecked = true;
            MouseEvent(popup.Child, Mouse.PreviewMouseDownEvent);
            Require(toggle.IsChecked == true, "Clicking inside Explorer must not dismiss it.");
            var outsideClick = MouseEvent(outside, Mouse.PreviewMouseDownEvent);
            Require(toggle.IsChecked == false && !outsideClick.Handled, "An outside click must close Explorer without swallowing the click.");
            foreach (var ownerEvent in new[] { "OnDeactivated", "OnLocationChanged" })
            {
                toggle.IsChecked = true;
                typeof(Window).GetMethod(ownerEvent, methods)!.Invoke(owner, new object[] { EventArgs.Empty });
                Require(toggle.IsChecked == false, $"Explorer remained open after {ownerEvent}.");
            }
            toggle.IsChecked = true;
            Require(toggle.Focusable && toggle.IsTabStop && toggle.IsEnabled,
                "The open Explorer toggle lost keyboard support.");
            var escape = new KeyEventArgs(Keyboard.PrimaryDevice, source, Environment.TickCount, Key.Escape)
                { RoutedEvent = Keyboard.PreviewKeyDownEvent };
            popup.Child.RaiseEvent(escape);
            Require(escape.Handled && toggle.IsChecked == false && toggle.IsHitTestVisible,
                "Escape must close Explorer and restore its opening button.");
            toggle.IsChecked = true;
            var toggleEscape = new KeyEventArgs(Keyboard.PrimaryDevice, source, Environment.TickCount, Key.Escape)
                { RoutedEvent = Keyboard.PreviewKeyDownEvent };
            toggle.RaiseEvent(toggleEscape);
            Require(toggleEscape.Handled && toggle.IsChecked == false,
                "Escape must also close the menu when focus remains on its toggle, including during a save.");
            toggle.IsChecked = true;
            menu.RaiseEvent(new RoutedEventArgs(FrameworkElement.UnloadedEvent));
            Require(toggle.IsChecked == false, "Unloading Explorer must close the popup.");
        }
        finally
        {
            owner.Content = null;
            owner.Close();
        }
    }

    private static void CheckExplorerSelection()
    {
        var menu = new ExplorerMenu();
        var popup = (Popup)menu.FindName("MenuPopup");
        var root = (FrameworkElement)popup.Child;
        popup.Child = null;
        root.DataContext = menu.DataContext;
        // Attach a presentation source so WPF loads the virtualized list's item
        // sizes. Without one, BringIntoView cannot reproduce the scrolling bug.
        // WS_VISIBLE is absent: this host is never shown or activated.
        using var source = new HwndSource(new HwndSourceParameters("Explorer headless layout")
            { WindowStyle = 0, Width = 680, Height = 540 }) { RootVisual = root };
        void Layout()
        {
            root.Measure(new Size(680, double.PositiveInfinity));
            root.Arrange(new Rect(new Point(), root.DesiredSize));
            root.UpdateLayout();
        }
        void Flush()
        {
            Dispatcher.CurrentDispatcher.Invoke(DispatcherPriority.ContextIdle, new Action(Layout));
        }
        Layout();
        Flush();
        var list = Descendants(root).OfType<ListBox>().Single(b => AutomationProperties.GetName(b) == "Available calculations");
        var scroll = Descendants(list).OfType<ScrollViewer>().Single();
        var viewport = Descendants(scroll).OfType<ScrollContentPresenter>().Single();
        Rect Bounds(FrameworkElement element) => element.TransformToAncestor(viewport).TransformBounds(new Rect(element.RenderSize));
        var requested = new List<CalculationOperation>();
        menu.OperationRequested += requested.Add;
        foreach (var clickTarget in new[] { "button", "text", "row padding" })
        {
            scroll.ScrollToTop();
            Flush();
            var button = Descendants(viewport).OfType<Button>().First(b => Bounds(b).Top < viewport.ActualHeight && Bounds(b).Bottom > viewport.ActualHeight);
            var row = (ListBoxItem)ItemsControl.ContainerFromElement(list, button);
            var target = clickTarget == "button" ? (UIElement)button : clickTarget == "text"
                ? Descendants(button).OfType<TextBlock>().First() : row;
            var offset = scroll.VerticalOffset;
            var bounds = Bounds(button);
            target.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = Mouse.PreviewMouseDownEvent });
            // Exercise both immediate and deferred focus requests before release.
            button.BringIntoView();
            row.Dispatcher.BeginInvoke(DispatcherPriority.Loaded, new Action(() => row.BringIntoView()));
            Flush();
            Require(scroll.VerticalOffset == offset && Bounds(button) == bounds,
                $"Clicking Explorer {clickTarget} moved the operation before mouse-up: offset {offset} -> {scroll.VerticalOffset}.");

            target.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
                { RoutedEvent = Mouse.PreviewMouseUpEvent });
            if (clickTarget != "row padding")
            {
                var count = requested.Count;
                button.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                Require(requested.Count == count + 1 && requested[^1] == button.Tag,
                    "An Explorer click must request the original operation exactly once.");
            }
            // Keyboard/programmatic focus must still reveal a clipped operation
            // after the mouse input turn has finished.
            Require(button.Focusable && button.IsTabStop, "Explorer operations lost keyboard focus support.");
            button.BringIntoView();
            Flush();
            Require(scroll.VerticalOffset > offset, "Explorer focus scrolling remained blocked after the click.");
        }

        scroll.ScrollToTop();
        Flush();
        var wheel = new MouseWheelEventArgs(Mouse.PrimaryDevice, Environment.TickCount, -120) { RoutedEvent = Mouse.MouseWheelEvent };
        viewport.RaiseEvent(wheel);
        Flush();
        Require(wheel.Handled && scroll.VerticalOffset > 0, "Explorer mouse-wheel scrolling is broken.");
        var wheelOffset = scroll.VerticalOffset;
        var scrollbar = Descendants(scroll).OfType<ScrollBar>().Single(b => b.Orientation == Orientation.Vertical);
        scrollbar.RaiseEvent(new MouseButtonEventArgs(Mouse.PrimaryDevice, Environment.TickCount, MouseButton.Left)
            { RoutedEvent = Mouse.PreviewMouseDownEvent });
        ScrollBar.PageDownCommand.Execute(null, scroll);
        Flush();
        Require(scroll.VerticalOffset > wheelOffset, "Explorer scrollbar clicks no longer scroll.");

        var model = (ExplorerMenuViewModel)menu.DataContext;
        var categories = Descendants(root).OfType<ListBox>().Single(b => AutomationProperties.GetName(b) == "Calculation categories");
        foreach (var filter in new Action[]
        {
            () => model.Search = "curve",
            () => model.Search = "",
            () => categories.SelectedItem = "Ranks and arithmetic",
            () => categories.SelectedIndex = 0
        })
        {
            scroll.ScrollToEnd();
            Flush();
            Require(scroll.VerticalOffset > 0, "The filter regression must start with a scrolled list.");
            filter();
            Flush();
            Require(scroll.VerticalOffset == 0, "Changing or clearing an Explorer filter retained the old scroll position.");
        }
    }
    private static void CheckConfirmationDialogs()
    {
        foreach (var reset in new[] { false, true })
        foreach (var action in new[] { "cancel", "confirm", "close" })
        {
            var dialog = new ConfirmationWindow(reset ? "Reset equation?" : "Clear history?",
                reset ? "Restore the classic curve and recenter the plot. Your current equation will be replaced."
                    : "Remove all calculation results from this session? You can restore them with Edit → Undo.",
                reset ? "Reset equation" : "Clear history", reset ? "y^2 = x^3 - x" : null);
            try
            {
                Require(dialog.WindowStyle == WindowStyle.None && dialog.AllowsTransparency && !dialog.ShowInTaskbar,
                    "Confirmation must use the custom window frame.");
                Require(dialog.WindowStartupLocation == WindowStartupLocation.CenterOwner, "Confirmation must center on its owner.");
                var root = (FrameworkElement)dialog.Content;
                root.Measure(new Size(dialog.Width, double.PositiveInfinity));
                root.Arrange(new Rect(new Point(), root.DesiredSize));
                root.UpdateLayout();
                var cancel = (Button)dialog.FindName("CancelButton");
                var confirm = (Button)dialog.FindName("ConfirmButton");
                var detail = (Border)dialog.FindName("DetailPanel");
                Require(cancel.IsDefault && cancel.IsCancel && !confirm.IsDefault, "Enter and Escape must default to cancellation.");
                Require(FocusManager.GetFocusedElement(dialog) == root, "Initial confirmation focus must not highlight an action button.");
                Require((detail.Visibility == Visibility.Visible) == reset, "Only Reset should show the equation preview.");
                Rect Bounds(FrameworkElement element) => element.TransformToAncestor(root).TransformBounds(new Rect(element.RenderSize));
                Require(Bounds(confirm).Right < Bounds(cancel).Left && new Rect(root.RenderSize).Contains(Bounds(cancel))
                    && new Rect(root.RenderSize).Contains(Bounds(confirm)),
                    "Cancel must follow the confirmation action, with both buttons inside the window.");
                if (action == "cancel")
                {
                    var bitmap = new RenderTargetBitmap((int)root.ActualWidth, (int)root.ActualHeight, 96, 96, PixelFormats.Pbgra32);
                    bitmap.Render(root);
                    var encoder = new PngBitmapEncoder();
                    encoder.Frames.Add(BitmapFrame.Create(bitmap));
                    using var stream = System.IO.File.Create(System.IO.Path.Combine(AppContext.BaseDirectory, reset ? "confirmation-reset.png" : "confirmation-clear.png"));
                    encoder.Save(stream);
                }
                var button = action == "confirm" ? confirm : action == "cancel" ? cancel
                    : Descendants(root).OfType<Button>().Single(b => AutomationProperties.GetName(b) == "Close confirmation");
                button.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
                Require(dialog.Confirmed == (action == "confirm"), "The custom dialog returned the wrong confirmation result.");
            }
            finally { dialog.Close(); }
        }
    }
    private static void CheckConfirmations(CalculationRequest request)
    {
        using var workbench = new WorkbenchViewModel();
        var first = new CalculationJobViewModel(request, "First");
        var second = new CalculationJobViewModel(request, "Second");
        workbench.Jobs.Add(first);
        workbench.Jobs.Add(second);
        workbench.Selected = second;
        var clearAnswer = false;
        var clearConfirmations = 0;
        var panel = new ResultsPanel(() => { clearConfirmations++; return clearAnswer; }) { DataContext = workbench };
        panel.Measure(new Size(300, 800));
        panel.Arrange(new Rect(0, 0, 300, 800));
        panel.UpdateLayout();
        var clear = Descendants(panel).OfType<Button>().Single(b => Equals(b.Content, "Clear"));
        var collapse = Descendants(panel).OfType<Button>().Single(b => AutomationProperties.GetName(b) == "Collapse results panel");
        var title = Descendants(panel).OfType<TextBlock>().Single(t => t.Text == "Results");
        Rect Bounds(FrameworkElement element) => element.TransformToAncestor(panel).TransformBounds(new Rect(element.RenderSize));
        Require(Bounds(title).Right < Bounds(clear).Left && Bounds(clear).Right < Bounds(collapse).Left,
            "Clear overlaps the Results title or collapse button at minimum width.");
        Require(clear.IsEnabled, "Clear must be enabled when history is available.");
        clear.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
        Require(clearConfirmations == 1 && workbench.Jobs.Count == 2 && workbench.Selected == second,
            "Rejecting Clear changed the history or selected result.");
        clearAnswer = true;
        clear.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
        panel.UpdateLayout();
        Require(clearConfirmations == 2 && !workbench.HasResults && !workbench.HasSelection && !clear.IsEnabled,
            "Accepting Clear did not empty and disable the history.");
        clear.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
        Require(clearConfirmations == 2, "An empty history must not request confirmation.");

        var resetAnswer = false;
        var resetConfirmations = 0;
        var window = CreateMainWindow(() => { resetConfirmations++; return resetAnswer; });
        try
        {
            var root = (FrameworkElement)window.Content;
            root.Measure(new Size(1438, 918));
            root.Arrange(new Rect(0, 0, 1438, 918));
            root.UpdateLayout();
            var plotted = Descendants(root).OfType<TextBlock>().Single(t => t.Inlines.OfType<Run>().Any(r => r.Text.StartsWith("EQUATION OVER ℚ", StringComparison.Ordinal)));
            var equationPreview = plotted.Inlines.OfType<Run>().Last();
            Require(equationPreview.Text == window.ViewModel.Snapshot.Equation, "The initial plotted equation is missing.");
            window.ViewModel.ShowPoints = false;
            window.ViewModel.Equation.Text = "y^2 + y = x^3 - x";
            window.ViewModel.FlushUpdate();
            root.UpdateLayout();
            Require(equationPreview.Text == window.ViewModel.Snapshot.Equation, "The plotted equation did not update after editing.");
            window.Workbench.Jobs.Add(first);
            window.Workbench.Selected = first;
            var snapshot = window.ViewModel.Snapshot;
            var resets = 0;
            window.ViewModel.CurveResetRequested += (_, _) => resets++;
            var reset = Descendants(root).OfType<Button>().Single(b => Equals(b.Content, "Reset"));
            Require(reset.Command == null, "Reset must not bypass confirmation through a command binding.");
            reset.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(resetConfirmations == 1 && window.ViewModel.Snapshot == snapshot && resets == 0,
                "Rejecting Reset changed the equation or viewport.");
            resetAnswer = true;
            reset.RaiseEvent(new RoutedEventArgs(Button.ClickEvent));
            Require(resetConfirmations == 2 && window.ViewModel.Equation.Text == "y^2 = x^3 - x" && resets == 1,
                "Accepting Reset did not restore and recenter the classic curve.");
            root.UpdateLayout();
            Require(equationPreview.Text == window.ViewModel.Snapshot.Equation, "The plotted equation did not update after Reset.");
            Require(window.Workbench.Selected == first && window.Workbench.Jobs.Count == 1,
                "Reset must preserve calculation history.");
        }
        finally { window.Close(); window.Workbench.Dispose(); window.ViewModel.Dispose(); }
    }
    private static void CheckExportBounds()
    {
        // Keep the source attached, with both a row offset and a margin, as in MainWindow.
        // An origin-only rendering check misses the transparent strips and cropped content.
        var target = new Grid { Width = 200.25, Height = 140.25, Margin = new Thickness(8, 0, 0, 0) };
        target.Children.Add(new Border { Width = 12, Height = 12, Background = Brushes.Lime,
            HorizontalAlignment = HorizontalAlignment.Left, VerticalAlignment = VerticalAlignment.Top });
        target.Children.Add(new Border { Width = 12, Height = 12, Background = Brushes.Gold,
            HorizontalAlignment = HorizontalAlignment.Right, VerticalAlignment = VerticalAlignment.Bottom });
        var host = new Grid();
        host.RowDefinitions.Add(new RowDefinition { Height = new GridLength(59) });
        host.RowDefinitions.Add(new RowDefinition());
        Grid.SetRow(target, 1);
        host.Children.Add(target);
        host.Measure(new Size(220, 220));
        host.Arrange(new Rect(0, 0, 220, 220));
        host.UpdateLayout();
        var originalOffset = VisualTreeHelper.GetOffset(target);
        var originalSize = target.RenderSize;
        var bitmap = PlotImageExporter.Render(target);
        Require(bitmap.PixelWidth == 401 && bitmap.PixelHeight == 281, "Export must round fractional dimensions up.");
        var pixels = new byte[bitmap.PixelWidth * bitmap.PixelHeight * 4];
        bitmap.CopyPixels(pixels, bitmap.PixelWidth * 4, 0);
        Color Pixel(int x, int y)
        {
            var i = (y * bitmap.PixelWidth + x) * 4;
            return Color.FromArgb(pixels[i + 3], pixels[i + 2], pixels[i + 1], pixels[i]);
        }
        Require(Pixel(4, 4) == Colors.Lime && Pixel(390, 270) == Colors.Gold,
            "Export shifted or cropped the source's corner markers.");
        Require(Pixel(200, 140) == Color.FromRgb(0x12, 0x19, 0x20)
            && pixels.Where((_, i) => i % 4 == 3).All(alpha => alpha != 0),
            "The export canvas must have a dark background, including the rounded edge pixels.");
        Require(VisualTreeHelper.GetParent(target) == host && VisualTreeHelper.GetOffset(target) == originalOffset
            && target.RenderSize == originalSize, "Export must not move or resize the live view.");
    }

    private static void CheckPlotRendering()
    {
        foreach (var equation in new[] { "y^2=x^3-x", "y^2+xy+y=x^3-x", "y^2=x^3", "y^2=x^3-3*x-2", "y^2=x^3+1e400*x" })
        {
            Require(CurveEquationText.TryParse(equation, out var curve, out var error), error);
            var plot = new CurvePlot { Snapshot = new CurveSnapshot(curve!) };
            var host = new Border { Padding = new Thickness(8, 59, 8, 0), Child = plot };
            host.Measure(new Size(816, 559));
            host.Arrange(new Rect(0, 0, 816, 559));
            plot.Fit();
            plot.Zoom(0.8);
            plot.Zoom(1.25);
            plot.UpdateLayout();
            var bitmap = PlotImageExporter.Render(plot);
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

    private static void CheckExtremePlotViews()
    {
        var window = CreateMainWindow();
        try
        {
            foreach (var equation in new[]
            {
                "y^2 + 1e308*x*y = x^3 - 2.5e615*x^2 - 3",
                "y^2 - 1e308*x*y = x^3 - 2.5e615*x^2 - 3",
                "y^2 = x^3 + 1e308*x^2"
            })
            {
                window.RestoreSession(ExplorerSession.New() with { Equation = equation });
                SettleSession(window);
                var plot = (CurvePlot)window.FindName("Plot");
                Require(window.ViewModel.Snapshot.IsPlotUnavailable, "Out-of-range curves must show the plot precision message.");
                Require(!string.IsNullOrEmpty(window.ViewModel.Snapshot.Discriminant), "Unavailable plots must retain exact invariants.");
                window.ViewModel.FitCommand.Execute(null);
                SettleSession(window);
                plot.Fit();
                plot.Zoom(0.8);
                plot.Zoom(1.25);
                var state = plot.CaptureView();
                Require(double.IsFinite(state.CenterX) && double.IsFinite(state.CenterY)
                    && double.IsFinite(state.VerticalSpan) && state.VerticalSpan > 0, "Fit or Zoom produced an invalid viewport.");
                var image = PlotImageExporter.Render(plot);
                Require(image.PixelWidth > 0 && image.PixelHeight > 0, "Export failed after an out-of-range equation.");
            }
            window.RestoreSession(ExplorerSession.New());
            SettleSession(window);
            Require(!window.ViewModel.Snapshot.IsPlotUnavailable, "A normal curve must recover after an unavailable plot.");
        }
        finally { window.Close(); }
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
        var window = CreateMainWindow();
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

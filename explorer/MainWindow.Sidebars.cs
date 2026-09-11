using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Animation;
using EllipticCurves.Explorer.Models;
using EllipticCurves.Explorer.Windowing;

namespace EllipticCurves.Explorer;

public partial class MainWindow
{
    private const double SidebarTabWidth = 32;
    private readonly SidebarState equationSidebar = new(238);
    private readonly SidebarState resultsSidebar = new(300);

    private sealed class SidebarState(double minimumWidth)
    {
        public double MinimumWidth { get; } = minimumWidth;
        public double ExpandedWidth { get; set; } = minimumWidth;
        public bool IsVisible { get; set; } = true;
        public bool IsAnimating { get; set; }
        public int AnimationVersion { get; set; }
    }

    private void UpdateSidebarBounds()
    {
        var available = Workspace.ActualWidth - PlotColumn.MinWidth
            - Workspace.ColumnDefinitions[1].Width.Value - Workspace.ColumnDefinitions[3].Width.Value;
        // Resolve both limits from requested widths so resizing the window cannot
        // make the panels repeatedly push each other's measured width back and forth.
        var equationWidth = Math.Min(EquationColumn.Width.Value,
            Math.Max(EquationColumn.MinWidth, available - ResultsColumn.MinWidth));
        ResultsColumn.MaxWidth = Math.Max(resultsSidebar.MinimumWidth, available - equationWidth);
        var resultsWidth = Math.Min(ResultsColumn.Width.Value, ResultsColumn.MaxWidth);
        EquationColumn.MaxWidth = Math.Max(equationSidebar.MinimumWidth, available - resultsWidth);
        QueueSessionStatusRefresh();
    }

    private void ExpandResultsClick(object sender, RoutedEventArgs e) => SetResultsVisible(true);
    private void CollapseEquationClick(object sender, RoutedEventArgs e) => SetEquationVisible(false);
    private void ExpandEquationClick(object sender, RoutedEventArgs e) => SetEquationVisible(true);
    private void SetEquationVisible(bool visible) =>
        SetSidebarVisible(equationSidebar, visible, EquationColumn, EquationPanel, EquationTab, HorizontalAlignment.Right, EquationSplitter);
    private void SetResultsVisible(bool visible) =>
        SetSidebarVisible(resultsSidebar, visible, ResultsColumn, Results, ResultsTab, HorizontalAlignment.Left, ResultsSplitter);

    private void SetSidebarVisible(SidebarState state, bool visible, ColumnDefinition column,
        FrameworkElement panel, Button tab, HorizontalAlignment slideAlignment, GridSplitter? splitter = null)
    {
        if (state.IsVisible == visible) return;
        var from = column.ActualWidth > 0 ? column.ActualWidth : column.Width.Value;
        if (!visible && !state.IsAnimating) state.ExpandedWidth = Math.Max(state.MinimumWidth, from);
        var to = visible ? Math.Min(state.ExpandedWidth, column.MaxWidth) : SidebarTabWidth;
        var version = ++state.AnimationVersion;
        state.IsVisible = visible;
        QueueSessionStatusRefresh();
        column.BeginAnimation(ColumnDefinition.WidthProperty, null);
        column.MinWidth = SidebarTabWidth;
        column.Width = new GridLength(to);
        // Keep content at its expanded width while clipping it toward the outer edge.
        panel.Width = Math.Max(state.MinimumWidth, visible ? to : from);
        panel.HorizontalAlignment = slideAlignment;
        panel.Visibility = Visibility.Visible;
        tab.Visibility = Visibility.Collapsed;
        if (splitter != null) splitter.Visibility = Visibility.Collapsed;

        void Finish()
        {
            if (version != state.AnimationVersion) return;
            column.BeginAnimation(ColumnDefinition.WidthProperty, null);
            column.MinWidth = visible ? state.MinimumWidth : SidebarTabWidth;
            panel.Width = double.NaN;
            panel.HorizontalAlignment = HorizontalAlignment.Stretch;
            panel.Visibility = visible ? Visibility.Visible : Visibility.Collapsed;
            if (splitter != null) splitter.Visibility = panel.Visibility;
            tab.Visibility = visible ? Visibility.Collapsed : Visibility.Visible;
            state.IsAnimating = false;
        }

        if (!IsLoaded || !SystemParameters.ClientAreaAnimation)
        {
            Finish();
            return;
        }
        state.IsAnimating = true;
        var animation = new GridLengthAnimation { From = from, To = to, Duration = TimeSpan.FromMilliseconds(200) };
        animation.Completed += (_, _) => Finish();
        column.BeginAnimation(ColumnDefinition.WidthProperty, animation, HandoffBehavior.SnapshotAndReplace);
    }
    private static SidebarSession CaptureSidebar(SidebarState state, ColumnDefinition column) => new(state.IsVisible,
        Math.Max(state.MinimumWidth, state.IsVisible && !state.IsAnimating ? column.Width.Value : state.ExpandedWidth));

    private static void RestoreSidebar(SidebarState state, SidebarSession saved, ColumnDefinition column,
        FrameworkElement panel, Button tab, GridSplitter splitter)
    {
        state.AnimationVersion++;
        state.IsAnimating = false;
        state.IsVisible = saved.Visible;
        state.ExpandedWidth = saved.Width;
        column.BeginAnimation(ColumnDefinition.WidthProperty, null);
        column.MinWidth = saved.Visible ? state.MinimumWidth : SidebarTabWidth;
        column.MaxWidth = double.PositiveInfinity;
        column.Width = new GridLength(saved.Visible ? saved.Width : SidebarTabWidth);
        panel.Width = double.NaN;
        panel.HorizontalAlignment = HorizontalAlignment.Stretch;
        panel.Visibility = splitter.Visibility = saved.Visible ? Visibility.Visible : Visibility.Collapsed;
        tab.Visibility = saved.Visible ? Visibility.Collapsed : Visibility.Visible;
    }
}

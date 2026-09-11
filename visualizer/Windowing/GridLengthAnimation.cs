using System.Windows;
using System.Windows.Media.Animation;

namespace EllipticCurves.Visualizer.Windowing;

/// <summary>Animates a pixel-sized dock column while its contents slide inside a clip.</summary>
public sealed class GridLengthAnimation : AnimationTimeline
{
    public static readonly DependencyProperty FromProperty = DependencyProperty.Register(nameof(From), typeof(double), typeof(GridLengthAnimation));
    public static readonly DependencyProperty ToProperty = DependencyProperty.Register(nameof(To), typeof(double), typeof(GridLengthAnimation));
    public double From { get => (double)GetValue(FromProperty); set => SetValue(FromProperty, value); }
    public double To { get => (double)GetValue(ToProperty); set => SetValue(ToProperty, value); }
    public override Type TargetPropertyType => typeof(GridLength);
    protected override Freezable CreateInstanceCore() => new GridLengthAnimation();
    public override object GetCurrentValue(object defaultOriginValue, object defaultDestinationValue, AnimationClock animationClock)
    {
        var progress = animationClock.CurrentProgress ?? 0;
        var eased = 1 - Math.Pow(1 - progress, 3);
        return new GridLength(Math.Max(0, From + (To - From) * eased));
    }
}

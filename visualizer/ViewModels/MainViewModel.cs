#nullable enable
using EllipticCurves.Visualizer.Models;

namespace EllipticCurves.Visualizer.ViewModels;

public sealed class MainViewModel : ObservableObject, IDisposable
{
    private bool updating, disposed;
    private bool showPoints = true;
    private bool showGrid = true;
    private CancellationTokenSource? sampleCancellation;
    private CurvePreset? selectedPreset;
    private IReadOnlyList<EllipticCurvePoint> samples = Array.Empty<EllipticCurvePoint>();
    private string sampleStatus = "";
    public IReadOnlyList<CoefficientViewModel> Coefficients { get; }
    public IReadOnlyList<CurvePreset> Presets => CurvePreset.All;
    public CurveSnapshot Snapshot { get; private set; } = new(new EllipticCurveQ(0, 0, 0, -1, 0));
    public bool HasInputError => Coefficients.Any(value => !value.IsValid);
    public string InputStatus => HasInputError ? "Incomplete input · showing the last valid curve" : "Exact arithmetic · coefficients in ℚ";
    public string PresetDescription => selectedPreset?.Description ?? "Your own coefficients";
    public Task PendingSamples { get; private set; } = Task.CompletedTask;
    public event EventHandler? ViewResetRequested;
    public RelayCommand ResetCommand { get; }
    public RelayCommand FitCommand { get; }

    public MainViewModel()
    {
        Coefficients = new[]
        {
            new CoefficientViewModel("a₁", "xy", CoefficientsChanged),
            new CoefficientViewModel("a₂", "x²", CoefficientsChanged),
            new CoefficientViewModel("a₃", "y", CoefficientsChanged),
            new CoefficientViewModel("a₄", "x", CoefficientsChanged),
            new CoefficientViewModel("a₆", "constant", CoefficientsChanged)
        };
        ResetCommand = new RelayCommand(_ => ApplyPreset(CurvePreset.All[0]));
        FitCommand = new RelayCommand(_ => ViewResetRequested?.Invoke(this, EventArgs.Empty));
        ApplyPreset(CurvePreset.All[0]);
    }

    public CurvePreset? SelectedPreset
    {
        get => selectedPreset;
        set { if (value is not null && value != selectedPreset) ApplyPreset(value); }
    }

    public bool ShowGrid
    {
        get => showGrid;
        set { showGrid = value; OnPropertyChanged(); }
    }

    public bool ShowPoints
    {
        get => showPoints;
        set { showPoints = value; OnPropertyChanged(); RefreshSamples(); }
    }

    public IReadOnlyList<EllipticCurvePoint> Samples
    {
        get => samples;
        private set { samples = value; OnPropertyChanged(); }
    }

    public string SampleStatus
    {
        get => sampleStatus;
        private set { sampleStatus = value; OnPropertyChanged(); }
    }

    public void ApplyPreset(CurvePreset preset)
    {
        updating = true;
        try
        {
            var values = new[] { preset.A1, preset.A2, preset.A3, preset.A4, preset.A6 };
            for (var i = 0; i < values.Length; i++) Coefficients[i].Value = values[i];
        }
        finally { updating = false; }
        Recalculate();
        selectedPreset = preset;
        OnPropertyChanged(nameof(SelectedPreset));
        OnPropertyChanged(nameof(PresetDescription));
        ViewResetRequested?.Invoke(this, EventArgs.Empty);
    }

    private void CoefficientsChanged()
    {
        if (updating || disposed) return;
        selectedPreset = null;
        OnPropertyChanged(nameof(SelectedPreset));
        OnPropertyChanged(nameof(PresetDescription));
        Recalculate();
    }

    private void Recalculate()
    {
        OnPropertyChanged(nameof(HasInputError));
        OnPropertyChanged(nameof(InputStatus));
        if (HasInputError) return;
        Snapshot = new CurveSnapshot(new EllipticCurveQ(Coefficients[0].ExactValue, Coefficients[1].ExactValue,
            Coefficients[2].ExactValue, Coefficients[3].ExactValue, Coefficients[4].ExactValue));
        OnPropertyChanged(nameof(Snapshot));
        RefreshSamples();
    }

    private void RefreshSamples()
    {
        sampleCancellation?.Cancel();
        sampleCancellation?.Dispose();
        sampleCancellation = null;
        Samples = Array.Empty<EllipticCurvePoint>();
        if (disposed || !showPoints || Snapshot.IsSingular)
        {
            SampleStatus = Snapshot.IsSingular ? "Singular cubic · rational samples disabled" : "Rational samples hidden";
            PendingSamples = Task.CompletedTask;
            return;
        }
        sampleCancellation = new CancellationTokenSource();
        SampleStatus = "Finding rational samples…";
        PendingSamples = FindSamplesAsync(Snapshot.Curve, sampleCancellation.Token);
    }

    private async Task FindSamplesAsync(EllipticCurveQ curve, CancellationToken token)
    {
        try
        {
            // A short debounce keeps rapid slider changes from queuing point searches.
            await Task.Delay(120, token);
            var points = await Task.Run(() =>
            {
                var result = new List<EllipticCurvePoint>();
                foreach (var point in curve.RationalPoints(12, 4))
                {
                    token.ThrowIfCancellationRequested();
                    if (!point.IsInfinity) result.Add(point);
                }
                return result.ToArray();
            }, token);
            token.ThrowIfCancellationRequested();
            Samples = points;
            SampleStatus = $"{points.Length} affine rational samples · x = m/n, |m| ≤ 12, 1 ≤ n ≤ 4";
        }
        catch (OperationCanceledException) when (token.IsCancellationRequested) { }
        catch (Exception) when (!token.IsCancellationRequested)
        {
            Samples = Array.Empty<EllipticCurvePoint>();
            SampleStatus = "Rational sample search unavailable";
        }
    }

    public void Dispose()
    {
        if (disposed) return;
        disposed = true;
        sampleCancellation?.Cancel();
        sampleCancellation?.Dispose();
        sampleCancellation = null;
    }
}

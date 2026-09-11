#nullable enable
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class MainViewModel : ObservableObject, IDisposable
{
    private bool updating, disposed;
    private bool showPoints = true;
    private bool showGrid = true;
    private bool isSimpleForm = true, isUpdatePending;
    private CancellationTokenSource? updateCancellation;
    private CancellationTokenSource? sampleCancellation;
    private BigRational? appliedSliderStep;
    private CurvePreset? selectedPreset;
    private IReadOnlyList<EllipticCurvePoint> samples = Array.Empty<EllipticCurvePoint>();
    private string sampleStatus = "";
    public IReadOnlyList<CoefficientViewModel> Coefficients { get; }
    public IReadOnlyList<CoefficientViewModel> SimpleCoefficients { get; }
    public IReadOnlyList<CoefficientViewModel> ActiveCoefficients => IsSimpleForm ? SimpleCoefficients : Coefficients;
    public CoefficientViewModel Step { get; }
    public EquationViewModel Equation { get; }
    public IReadOnlyList<CurvePreset> Presets => CurvePreset.All;
    public CurveSnapshot Snapshot { get; private set; } = new(new EllipticCurveQ(0, 0, 0, -1, 0));
    public bool HasInputError => Equation.Error.Length != 0 || ActiveCoefficients.Any(value => value.Error.Length != 0);
    public bool HasIncompleteInput => !Equation.IsValid || ActiveCoefficients.Any(value => !value.IsValid);
    public string InputStatus => HasInputError ? "Check the equation · showing the last valid curve"
        : HasIncompleteInput ? "Editing · showing the last valid curve"
        : isUpdatePending ? "Updating curve…" : "Exact arithmetic · coefficients in ℚ";
    public string PresetDescription => selectedPreset?.Description ?? "Your own coefficients";
    public Task PendingSamples { get; private set; } = Task.CompletedTask;
    public Task PendingUpdate { get; private set; } = Task.CompletedTask;
    public event EventHandler? CurveResetRequested;
    public event EventHandler? ViewResetRequested;
    public RelayCommand ResetCommand { get; }
    public RelayCommand FitCommand { get; }
    public RelayCommand SetStepCommand { get; }

    public MainViewModel()
    {
        Equation = new EquationViewModel(EquationChanged);
        Step = new CoefficientViewModel("Step", "", StepChanged, requirePositive: true);
        BigRational? CurrentStep() => Step.IsValid ? Step.ExactValue : null;
        Coefficients = new[]
        {
            new CoefficientViewModel("a₁", "xy", CoefficientsChanged, CurrentStep),
            new CoefficientViewModel("a₂", "x²", CoefficientsChanged, CurrentStep),
            new CoefficientViewModel("a₃", "y", CoefficientsChanged, CurrentStep),
            new CoefficientViewModel("a₄", "x", CoefficientsChanged, CurrentStep),
            new CoefficientViewModel("a₆", "constant", CoefficientsChanged, CurrentStep)
        };
        SimpleCoefficients = new[]
        {
            new CoefficientViewModel("A", "x", CoefficientsChanged, CurrentStep),
            new CoefficientViewModel("B", "constant", CoefficientsChanged, CurrentStep)
        };
        Step.SetExact(new BigRational(1, 100));
        updating = true;
        Coefficients[3].SetExact(-1);
        updating = false;
        ResetCommand = new RelayCommand(_ => ApplyPreset(CurvePreset.All[0]));
        FitCommand = new RelayCommand(_ => ViewResetRequested?.Invoke(this, EventArgs.Empty));
        SetStepCommand = new RelayCommand(value => { Step.Text = value?.ToString() ?? ""; Step.CommitEdit(); });
        ApplyPreset(CurvePreset.All[0]);
        // The initial snapshot already contains the classic curve.
        RefreshSamples();
    }

    public bool IsSimpleForm => isSimpleForm;
    public bool IsGeneralForm => !IsSimpleForm;
    public string FormDescription => IsSimpleForm ? "Simple Weierstrass form · A, B" : "General Weierstrass form · a₁, a₂, a₃, a₄, a₆";

    private void NotifyForm()
    {
        OnPropertyChanged(nameof(IsSimpleForm));
        OnPropertyChanged(nameof(IsGeneralForm));
        OnPropertyChanged(nameof(ActiveCoefficients));
        OnPropertyChanged(nameof(FormDescription));
    }

    private void EquationChanged()
    {
        if (updating || disposed) return;
        if (Equation.IsValid)
        {
            var curve = Equation.Curve;
            updating = true;
            try
            {
                isSimpleForm = curve.A1.IsZero && curve.A2.IsZero && curve.A3.IsZero;
                var values = IsSimpleForm ? new[] { curve.A4, curve.A6 }
                    : new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 };
                for (var i = 0; i < values.Length; i++)
                    if (!ActiveCoefficients[i].IsValid || ActiveCoefficients[i].ExactValue != values[i])
                        ActiveCoefficients[i].SetExact(values[i]);
            }
            finally { updating = false; }
            NotifyForm();
        }
        ScheduleUpdate();
    }

    private void SynchronizeEquation()
    {
        if (ActiveCoefficients.Any(value => !value.IsValid)) return;
        updating = true;
        try
        {
            Equation.SetCurve(IsSimpleForm
                ? new EllipticCurveQ(0, 0, 0, SimpleCoefficients[0].ExactValue, SimpleCoefficients[1].ExactValue)
                : new EllipticCurveQ(Coefficients[0].ExactValue, Coefficients[1].ExactValue,
                    Coefficients[2].ExactValue, Coefficients[3].ExactValue, Coefficients[4].ExactValue));
        }
        finally { updating = false; }
    }

    private void StepChanged()
    {
        var nextStep = Step.IsValid ? Step.ExactValue : (BigRational?)null;
        if (appliedSliderStep == nextStep) return;
        appliedSliderStep = nextStep;
        foreach (var coefficient in Coefficients.Concat(SimpleCoefficients)) coefficient.RecenterSlider();
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
        set
        {
            if (showPoints == value) return;
            showPoints = value;
            OnPropertyChanged();
            RefreshSamples();
        }
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
        CancelUpdate();
        updating = true;
        try
        {
            isSimpleForm = preset.A1 == 0 && preset.A2 == 0 && preset.A3 == 0;
            var values = new[] { preset.A1, preset.A2, preset.A3, preset.A4, preset.A6 };
            if (IsSimpleForm)
            {
                SimpleCoefficients[0].SetExact(preset.A4);
                SimpleCoefficients[1].SetExact(preset.A6);
            }
            else for (var i = 0; i < values.Length; i++) Coefficients[i].SetExact(values[i]);
        }
        finally { updating = false; }
        SynchronizeEquation();
        NotifyForm();
        Recalculate();
        selectedPreset = preset;
        OnPropertyChanged(nameof(SelectedPreset));
        OnPropertyChanged(nameof(PresetDescription));
        CurveResetRequested?.Invoke(this, EventArgs.Empty);
    }

    private void CoefficientsChanged()
    {
        if (updating || disposed) return;
        SynchronizeEquation();
        ScheduleUpdate();
    }

    private void ScheduleUpdate()
    {
        CancelUpdate();
        NotifyInput();
        // Validation and equivalent text edits do not change the plotted curve.
        if (!HasIncompleteInput && Snapshot.Curve.Equals(Equation.Curve)) return;
        selectedPreset = null;
        OnPropertyChanged(nameof(SelectedPreset));
        OnPropertyChanged(nameof(PresetDescription));
        if (HasIncompleteInput) return;
        updateCancellation = new CancellationTokenSource();
        isUpdatePending = true;
        OnPropertyChanged(nameof(InputStatus));
        PendingUpdate = UpdateAfterPauseAsync(updateCancellation.Token);
    }

    private async Task UpdateAfterPauseAsync(CancellationToken token)
    {
        try
        {
            await Task.Delay(300, token);
            token.ThrowIfCancellationRequested();
            isUpdatePending = false;
            Recalculate();
        }
        catch (OperationCanceledException) when (token.IsCancellationRequested) { }
    }

    public void FlushUpdate()
    {
        CancelUpdate();
        if (!disposed) Recalculate();
    }

    private void CancelUpdate()
    {
        updateCancellation?.Cancel();
        updateCancellation?.Dispose();
        updateCancellation = null;
        isUpdatePending = false;
        PendingUpdate = Task.CompletedTask;
    }

    private void NotifyInput()
    {
        OnPropertyChanged(nameof(HasInputError));
        OnPropertyChanged(nameof(HasIncompleteInput));
        OnPropertyChanged(nameof(InputStatus));
    }

    private void Recalculate()
    {
        NotifyInput();
        if (HasIncompleteInput || disposed || Snapshot.Curve.Equals(Equation.Curve)) return;
        Snapshot = new CurveSnapshot(Equation.Curve);
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
            // Keep rapid edits from queuing point searches.
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
        CancelUpdate();
        sampleCancellation?.Cancel();
        sampleCancellation?.Dispose();
        sampleCancellation = null;
    }
}

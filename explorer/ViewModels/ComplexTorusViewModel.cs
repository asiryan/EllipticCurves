#nullable enable
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

/// <summary>Lazy, cancellable preparation of a torus view; late results never replace a newer curve.</summary>
public sealed class ComplexTorusViewModel : ObservableObject, IDisposable
{
    private readonly SemaphoreSlim workerGate = new(1, 1);
    private readonly Func<EllipticCurveQ, CancellationToken, TorusLattice> prepare;
    private CancellationTokenSource? pending;
    private EllipticCurveQ? curve;
    private IReadOnlyList<EllipticCurvePoint> samples = Array.Empty<EllipticCurvePoint>();
    private bool active, disposed, samplesMapped;
    private int version;
    private TorusPoint? selectedPoint;
    private string? restoredSelection;
    private bool restoringSelection;
    public string? SessionSelection => restoredSelection ?? SelectedPoint?.Point.ToString();

    public void RestoreSelection(string? point)
    {
        restoredSelection = point ?? EllipticCurvePoint.Infinity.ToString();
    }
    public TorusLattice? Lattice { get; private set; }
    public IReadOnlyList<TorusPoint> Points { get; private set; } = Array.Empty<TorusPoint>();
    public bool IsBusy { get; private set; }
    public bool HasLattice => Lattice != null;
    public bool IsUnavailable => !HasLattice;
    public string Status { get; private set; } = "Choose Complex torus to compute the period lattice.";
    public Task PendingUpdate { get; private set; } = Task.CompletedTask;

    public TorusPoint? SelectedPoint
    {
        get => selectedPoint;
        set
        {
            if (Equals(selectedPoint, value)) return;
            if (!restoringSelection) restoredSelection = null;
            selectedPoint = value;
            OnPropertyChanged();
        }
    }

    public ComplexTorusViewModel() : this(TorusLattice.Create) { }
    internal ComplexTorusViewModel(Func<EllipticCurveQ, CancellationToken, TorusLattice> prepare) => this.prepare = prepare;

    public void Update(EllipticCurveQ nextCurve, IReadOnlyList<EllipticCurvePoint> nextSamples, bool isActive)
    {
        if (disposed) return;
        var curveChanged = !ReferenceEquals(curve, nextCurve);
        var pointsChanged = !samples.SequenceEqual(nextSamples);
        if (!curveChanged && !pointsChanged && active == isActive)
        {
            if (restoredSelection != null) SetPoints(Points);
            return;
        }
        CancelPending();
        curve = nextCurve;
        samples = nextSamples.ToArray();
        active = isActive;
        if (curveChanged)
        {
            samplesMapped = false;
            Lattice = null;
            SetPoints(Array.Empty<TorusPoint>());
        }
        else if (pointsChanged)
        {
            samplesMapped = false;
            SetPoints(Lattice == null ? Array.Empty<TorusPoint>() : new[] { TorusPoint.Origin });
        }
        if (!active)
        {
            NotifyState();
            return;
        }
        if (curve.IsSingular)
        {
            Status = "Δ = 0 · a singular cubic has no smooth complex torus.";
            NotifyState();
            return;
        }
        if (HasLattice && samplesMapped)
        {
            NotifyState();
            return;
        }
        StartUpdate();
    }

    private void StartUpdate()
    {
        var request = ++version;
        pending = new CancellationTokenSource();
        IsBusy = true;
        Status = HasLattice ? "Mapping rational samples…" : "Computing period lattice…";
        NotifyState();
        PendingUpdate = PrepareAsync(curve!, samples, request, pending.Token);
    }

    private async Task PrepareAsync(EllipticCurveQ input, IReadOnlyList<EllipticCurvePoint> inputSamples,
        int request, CancellationToken token)
    {
        try
        {
            await Task.Delay(180, token);
            // One worker per view, including when a previous calculation is still observing cancellation.
            await workerGate.WaitAsync(token);
            try
            {
                using var deadline = CancellationTokenSource.CreateLinkedTokenSource(token);
                deadline.CancelAfter(TimeSpan.FromSeconds(15));
                if (Lattice == null)
                {
                    var lattice = await Task.Run(() => prepare(input, deadline.Token), deadline.Token);
                    if (!IsCurrent(request, token)) return;
                    Lattice = lattice;
                    SetPoints(new[] { TorusPoint.Origin });
                    Status = "Mapping rational samples…";
                    NotifyState();
                }
                var basis = Lattice!;
                var mapped = await Task.Run(() => basis.MapPoints(inputSamples, deadline.Token), deadline.Token);
                if (!IsCurrent(request, token)) return;
                SetPoints(mapped.Points);
                samplesMapped = true;
                Status = mapped.Summary;
            }
            finally { workerGate.Release(); }
        }
        catch (OperationCanceledException) when (IsCurrent(request, token))
        {
            Status = HasLattice ? "Point mapping reached the 15 s limit. The period lattice is available."
                : "Period calculation reached the 15 s limit. Try a simpler curve.";
        }
        catch (OperationCanceledException) { }
        catch (Exception error) when (error is ArithmeticException or ArgumentException or InvalidOperationException)
        {
            if (IsCurrent(request, token)) Status = "Complex view unavailable: " + error.Message;
        }
        finally
        {
            if (IsCurrent(request, token))
            {
                IsBusy = false;
                NotifyState();
            }
        }
    }

    private bool IsCurrent(int request, CancellationToken token) => !disposed && request == version && !token.IsCancellationRequested;

    private void SetPoints(IReadOnlyList<TorusPoint> points)
    {
        var previousPoint = SelectedPoint?.Point;
        restoringSelection = true;
        try
        {
            Points = points;
            // The saved selection may arrive before background point mapping.
            OnPropertyChanged(nameof(Points));
            var restored = Points.FirstOrDefault(point => point.Point.ToString() == restoredSelection);
            SelectedPoint = restored ?? Points.FirstOrDefault(point => previousPoint.HasValue && point.Point.Equals(previousPoint.Value))
                ?? Points.FirstOrDefault();
            if (restored != null) restoredSelection = null;
        }
        finally { restoringSelection = false; }
    }

    private void NotifyState()
    {
        OnPropertyChanged(nameof(Lattice));
        OnPropertyChanged(nameof(HasLattice));
        OnPropertyChanged(nameof(IsUnavailable));
        OnPropertyChanged(nameof(IsBusy));
        OnPropertyChanged(nameof(Status));
    }

    private void CancelPending()
    {
        version++;
        pending?.Cancel();
        pending?.Dispose();
        pending = null;
        IsBusy = false;
    }

    public void Dispose()
    {
        if (disposed) return;
        disposed = true;
        CancelPending();
    }
}

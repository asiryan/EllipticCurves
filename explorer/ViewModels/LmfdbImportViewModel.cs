#nullable enable
using System.IO;
using System.Net.Http;
using System.Text.Json;
using EllipticCurves.Explorer.Computations;
using EllipticCurves.Explorer.Models;

namespace EllipticCurves.Explorer.ViewModels;

public sealed class LmfdbImportViewModel : ObservableObject, IDisposable
{
    private readonly Func<LmfdbConductorRange, LmfdbCurveFormula?, CancellationToken, Task<IReadOnlyList<LmfdbCurveFormula>>> search;
    private readonly List<LmfdbCurveFormula?> pageStarts = new();
    private CancellationTokenSource? running;
    private string query = "37";
    private bool disposed, hasNext;
    private int pageIndex;
    private LmfdbCurveFormula? selected;

    public LmfdbImportViewModel(Func<LmfdbConductorRange, LmfdbCurveFormula?, CancellationToken, Task<IReadOnlyList<LmfdbCurveFormula>>>? search = null)
        => this.search = search ?? new LmfdbCurveSearch().SearchAsync;

    public string Query
    {
        get => query;
        set
        {
            if (query == value) return;
            query = value;
            running?.Cancel();
            running = null;
            Curves = Array.Empty<LmfdbCurveFormula>();
            Selected = null;
            pageStarts.Clear();
            hasNext = false;
            pageIndex = 0;
            Error = "";
            Status = "Choose Search to load formulas from LMFDB.";
            NotifyState();
            OnPropertyChanged();
        }
    }

    public string InputError => LmfdbConductorRange.TryParse(Query, out _) ? "" : "Enter a conductor or range from 1 to 500000, such as 37 or 11-100.";
    public IReadOnlyList<LmfdbCurveFormula> Curves { get; private set; } = Array.Empty<LmfdbCurveFormula>();
    public LmfdbCurveFormula? Selected
    {
        get => selected;
        set { selected = value; OnPropertyChanged(); OnPropertyChanged(nameof(CanImport)); }
    }
    public bool IsBusy => running != null;
    public bool CanSearch => !disposed && !IsBusy && InputError.Length == 0;
    public bool CanImport => !disposed && !IsBusy && Selected != null && Curves.Contains(Selected);
    public bool CanPrevious => CanSearch && pageIndex > 0;
    public bool CanNext => CanSearch && hasNext;
    public string Status { get; private set; } = "Search by conductor, then select a formula to import.";
    public string Error { get; private set; } = "";
    public string PageInfo => Curves.Count == 0 ? "" : $"Page {pageIndex + 1} · {Curves.Count} formulas";

    public Task SearchAsync() => CanSearch ? LoadAsync(0, null, true) : Task.CompletedTask;
    public Task NextAsync() => CanNext ? LoadAsync(pageIndex + 1, Curves[^1], false) : Task.CompletedTask;
    public Task PreviousAsync() => CanPrevious ? LoadAsync(pageIndex - 1, pageStarts[pageIndex - 1], false) : Task.CompletedTask;
    public void Cancel() => running?.Cancel();

    private async Task LoadAsync(int targetPage, LmfdbCurveFormula? after, bool reset)
    {
        if (!LmfdbConductorRange.TryParse(Query, out var range)) return;
        using var cancellation = new CancellationTokenSource();
        running = cancellation;
        Error = "";
        Status = "Loading formulas from LMFDB…";
        NotifyState();
        try
        {
            var curves = await search(range!, after, cancellation.Token);
            cancellation.Token.ThrowIfCancellationRequested();
            if (disposed || running != cancellation) return;
            if (!reset && targetPage > pageIndex && curves.Count == 0)
            {
                hasNext = false;
                Status = "No more formulas in this conductor range.";
                return;
            }
            if (reset) pageStarts.Clear();
            if (targetPage == pageStarts.Count) pageStarts.Add(after);
            else pageStarts[targetPage] = after;
            pageIndex = targetPage;
            Curves = curves;
            Selected = curves.FirstOrDefault();
            hasNext = curves.Count == LmfdbCurveSearch.PageSize;
            Status = curves.Count == 0 ? $"No curves found for N = {range}." : $"N = {range} · select a formula to import.";
        }
        catch (OperationCanceledException)
        {
            if (running == cancellation) Status = cancellation.IsCancellationRequested
                ? "Search cancelled." : "LMFDB did not respond within 30 seconds. Try again.";
        }
        catch (Exception error) when (error is HttpRequestException or IOException or JsonException or FormatException or InvalidOperationException or KeyNotFoundException or OverflowException)
        {
            if (running == cancellation)
            {
                Status = "The formulas could not be loaded. You can retry the search.";
                Error = error is HttpRequestException or IOException ? "Could not connect to LMFDB. Check your connection and try again."
                    : error is FormatException ? error.Message : "LMFDB returned an unexpected response. Try again later.";
            }
        }
        finally
        {
            if (running == cancellation)
            {
                running = null;
                NotifyState();
            }
        }
    }

    private void NotifyState()
    {
        foreach (var name in new[] { nameof(Curves), nameof(IsBusy), nameof(CanSearch), nameof(CanImport), nameof(CanPrevious),
            nameof(CanNext), nameof(Status), nameof(Error), nameof(InputError), nameof(PageInfo) }) OnPropertyChanged(name);
    }

    public void Dispose()
    {
        disposed = true;
        running?.Cancel();
        running = null;
        NotifyState();
    }
}

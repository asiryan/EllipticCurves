# Certified factorization performance

Environment: 2026-09-16, Windows x64, .NET 8.0.31, 24 logical CPUs.
PARI/GP 2.19.0 development 31236-59c418aa0c uses one thread,
`factor_proven=1` and a 512 MB stack growth limit. Both implementations prove
primality; neither receives known factors. GP uses a fresh `ellinit` per run.

Times are wall-clock measurements. The first two tables use up to four native
workers and medians of three runs after the first run. Other tables contain
individual samples. Benchmarks ran without concurrent benchmark/test processes.

## Reported curve

```
y^2 = x^3 + x^2 - 221556180740323405132844117936*x
      + 35386140191724122461245294467670188433973860
```

| Conductor computation | First run | Median of subsequent runs |
| --- | ---: | ---: |
| Native ECM + SIQS | 731.5 ms | 572.2 ms |
| PARI/GP, proven factors | 828 ms | 827 ms |

These times exclude Explorer worker startup and UI updates.

## Independent inputs

The products of the prime pairs in `tests/Fixtures/factorization.csv` were
generated independently using a fixed GP seed. The generator proves each prime.
Cases include balanced semiprimes distinct from the reported curve's factors.

| Case | Decimal digits | Native median | GP median |
| --- | ---: | ---: | ---: |
| 1 | 38 | 19.4 ms | 8 ms |
| 2 | 42 | 54.4 ms | 25 ms |
| 3 | 46 | 58.2 ms | 77 ms |
| 4 | 50 | 162.0 ms | 196 ms |
| 5 | 54 | 337.9 ms | 477 ms |
| 6 | 57 | 768.3 ms | 1227 ms |
| 7 | 62 | 1890.9 ms | 3276 ms |
| 8 | 66 | 5889.4 ms | 8586 ms |

Factorization cost depends on the input, not just its decimal size.

## CPU scaling on the larger curve

```
y^2 + x*y = x^3 - 20820207864197471248300179976626*x
            + 36732936589138673862895758597955508398047757956
```

After removing small prime powers, the discriminant has a 75-digit composite
cofactor. The worker limit applies to recursive factorization and primality
proofs without nested worker groups.

Single-run measurements of the **complete conductor calculation**, in Release,
on the same machine with 24 logical CPUs:

| Explicit worker limit | Elapsed |
| --- | ---: |
| 4 | 80.77 s |
| 8 | 52.33 s |
| 12 | 37.05 s |
| 24 | 57.91 s |

The library's automatic mode (zero), which uses up to four workers for smaller
recursive subproblems, took 68.06 s in a separate run. The direct GP calculation
with one thread and `factor_proven=1` took 80.66 s. Every result was checked against
the same exact conductor; no known factors were provided to either factorizer.

Explorer defaults to four workers, bounded by available CPUs. More workers need
not be faster; 37.05 s is a single measurement, not an expected completion time.

To compare complete calculations at explicit worker limits:

```powershell
dotnet run --project tests/FactorizationBenchmark/FactorizationBenchmark.csproj -c Release -- --conductor-cpu 4 8 12 24
```

This prints wall time and total process CPU time, with a five-minute limit per
calculation. Add `0` for automatic worker selection; repeat a number for repeated
runs in the same process, for example `--conductor-cpu 12 12 12`.

## Explorer worker measurements

Follow-up runs on the larger curve, all with a limit of 12 workers:

| Execution | Elapsed |
| --- | ---: |
| Direct library, Release, two successive calls | 32.54 s; 53.76 s |
| Explorer worker, Debug, fresh process | 78.74 s |
| Explorer worker, Release, fresh process | 83.62 s |

The 12-worker setting reaches the factorizer correctly. Worker startup took
about 0.10 s and formatting/exit about 0.007 s. These overheads do not explain
the difference; its cause remains unresolved. Debug configuration alone does
not explain the observations either.

To measure the actual Explorer worker on Windows:

```powershell
dotnet build explorer/EllipticCurves.Explorer.csproj -c Release
dotnet run --project tests/FactorizationBenchmark/FactorizationBenchmark.csproj -c Release -- --conductor-worker explorer/bin/Release/net8.0-windows/EllipticCurves.Explorer.exe 12 3
```

The last arguments select workers and samples (defaults: 12 and 1). Each sample
uses a fresh process, checks the conductor and reports startup, calculation,
formatting/exit, total wall time and process CPU time in milliseconds. CPU time
sums work across threads; it is not elapsed time. Each sample has a five-minute limit.

## Reproduce the original comparison

From the repository root, with .NET 8 or newer installed:

```powershell
dotnet run --project tests/FactorizationBenchmark/FactorizationBenchmark.csproj -c Release
```

Run the GP benchmark separately, without simultaneous CPU-heavy work:

```powershell
& 'C:/path/to/gp.exe' -f -q tests/Fixtures/benchmark-factorization.gp
```

Both print CSV with the first run and median. Each C# run has a two-minute
cancellation limit. Tests check exact results and cancellation, not benchmark
times. See [algorithm details](native-arithmetic.md#integer-factorization).

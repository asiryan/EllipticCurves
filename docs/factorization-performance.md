# Certified factorization performance

Measured on 2026-09-16 on Windows x64, .NET 8.0.31, with 24 logical processors
available. SIQS uses at most **four workers** for inputs of at least 45 decimal
digits; the GP comparison explicitly uses **one thread**. These are wall-clock
measurements, not equal-CPU-work comparisons or guarantees for other inputs.

Both implementations prove primality of every returned factor. GP is PARI/GP
2.19.0 development 31236-59c418aa0c, with `factor_proven=1` and a stack growth
limit of 512 MB. Each measurement recomputes the result. GP conductor runs use
fresh `ellinit` objects to avoid the cached result of `ellglobalred`.

The reported median uses three runs after a first measured run. There is no
factorization cache, factor hint, database lookup or precomputed answer in the
native computation. Expected factors are used only to check the benchmark output.
Other benchmark/test processes were stopped during measurement.

## Reported curve

```
y^2 = x^3 + x^2 - 221556180740323405132844117936*x
      + 35386140191724122461245294467670188433973860
```

| Conductor computation | First run | Median of subsequent runs |
| --- | ---: | ---: |
| Native ECM + SIQS | 731.5 ms | 572.2 ms |
| PARI/GP, proven factors | 828 ms | 827 ms |

The earlier Pollard–Brent + basic MPQS implementation took approximately 5.3 s
on this curve. These figures exclude Explorer worker startup and UI updates.
The Explorer integration test also runs the same calculation in an actual worker
process and checks the complete conductor.

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

GP is faster on the smallest cases here. Parallel sieving gives the native
implementation lower wall-clock times on the larger cases. Factorization remains
input dependent; the same decimal size can have very different difficulty.

## Reproduce

From the repository root, with .NET 8 or newer installed:

```powershell
dotnet run --project tests/FactorizationBenchmark/FactorizationBenchmark.csproj -c Release
```

Run the GP benchmark separately, without simultaneous CPU-heavy work:

```powershell
& 'C:/path/to/gp.exe' -f -q tests/Fixtures/benchmark-factorization.gp
```

Both print CSV with the first run and median. Each C# run has a two-minute
cancellation limit. The timing is informational; unit tests check mathematical
results and use generous timeouts rather than asserting subsecond performance.
The factorization tests also exercise sequential SIQS, parallel SIQS, both ECM
stages, large perfect powers, repeated factors, pseudoprime rejection and
cancellation. See [algorithm details and references](native-arithmetic.md#integer-factorization).

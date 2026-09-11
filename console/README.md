# Elliptic Curves Console

A .NET 8 console example of the [EllipticCurves library](../README.md).
[Program.cs](Program.cs) studies the curve **y² = x³ − 17x² + 72x** from
[this paper](https://arxiv.org/abs/2510.11768).

The program prints the curve's invariants and torsion points, fetches LMFDB
metadata, and compares it with locally computed minimal models, rank bounds and
the conductor. It also estimates the analytic rank and attempts a rigorous
rank 0/1 certificate.

## Run

Install the .NET 8 SDK, then run this command from the repository root:

```sh
dotnet run --project console/EllipticCurves.Console.csproj
```

Alternatively, run `dotnet run` from the `console` directory. The project references
the library source directly; no separate NuGet installation is needed.

This example requires internet access for its LMFDB lookup. The native
computations themselves work offline after dependencies have been restored.

## Example output

The following excerpt shows the invariants, torsion points and comparison with
LMFDB. The program then prints the native analytic-rank result and its comparison
with the stored analytic rank.

```text
E: y^2 = x^3 - 17*x^2 + 72*x
Short Weierstrass: y^2 = x^3 - 73/3*x + 1190/27
Torsion: Z/2Z x Z/4Z
b2 = -68
b4 = 144
b6 = 0
b8 = -5184
D  = 82944
c4 = 1168
c6 = -38080
j  = 1556068/81
Torsion points:
O
(0, 0)
(8, 0)
(9, 0)
(6, 6)
(6, -6)
(12, 12)
(12, -12)
LMFDB: 48.a3
Url: https://www.lmfdb.org/EllipticCurve/Q/48.a3/
Minimal Weierstrass model: y^2 = x^3 + x^2 - 24*x + 36
Torsion: Z/2Z x Z/4Z
Rank(E) = 0
Analytic rank(E) = 0
Cond(E) = 48
Isomorphic to E: True
Native minimal Weierstrass model: y^2 = x^3 + x^2 - 24*x + 36
Native rank bounds(E) = 0
Exact native rank proved: True
Native Cond(E) = 48
Native minimal model matches LMFDB: True
Native conductor matches LMFDB: True
LMFDB rank is within native bounds: True
```

See the [library README](../README.md#rank-bounds-and-conductor) for the native
rank APIs and the meaning of their proof and certification statuses.

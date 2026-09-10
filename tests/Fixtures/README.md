# Independent arithmetic fixtures

Generated using PARI/GP 2.19.0 development 31236-59c418aa0c, compiled 2026-09-08,
from the [official Windows snapshot](https://pari.math.u-bordeaux.fr/download.html).
All inputs are explicit coefficients. The scripts use `ellglobalred`,
`ellminimalmodel`, `elllocalred`, `ellrank`, `ell2cover`, `ellrootno`, `ellan`, `ellanalyticrank`
and `lfun`; they do not load an elliptic-curve
database. See the [PARI elliptic-curve reference](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).

The committed CSV files contain unique rows encoded as UTF-8. Tests only read
these files; PARI/GP is not required to build or run the tests.

- `conductors.csv`: input a1,a2,a3,a4,a6; conductor; reduced minimal a1,a2,a3,a4,a6;
  PARI Kodaira codes at 2 and 3. The last two columns record coverage metadata.
- `ranks.csv`: a,b,proved rank for y²=x³+a*x²+b*x. Only cases where PARI's lower
  and upper rank bounds agree are included.
- `general-descent.csv`: a1,a2,a3,a4,a6; proved lower and upper algebraic rank
  bounds; dimension of Sel_2(E/Q), from `ellrank` and `ell2cover`. The first eight
  rows are fixed examples (ranks 0 through 4 and nontrivial Sha); the remainder
  contains small curves generated with a fixed seed, plus families with full or
  partial rational 2-torsion to exercise reducible resolvents. Tests compare the
  actual Selmer dimension on these small curves. Large fixed examples also document
  cases where quartic enumeration can exceed practical work limits.
- `root-numbers.csv`: input a1,a2,a3,a4,a6; global sign; local signs at 2 and 3,
  for the 5,975 conductor inputs.
- `local-roots.csv`: minimal a1,a2,a3,a4,a6; local signs at 2 and 3. The 1,472
  distinct models come from weighted short equations, including c4=0 and c6=0.
- `analytic.csv`: a1,a2,a3,a4,a6; conductor; sign; numerical analytic rank;
  leading derivative (not divided by rank factorial), evaluated with 60-digit
  precision. Includes ranks 0 through 4. PARI's `ellanalyticrank` output is a
  numerical reference, not itself a certificate of vanishing.
- `coefficients.csv`: one-based curve index into `analytic.csv`; n; a_n, for n<=100.
- `derivatives.csv`: curve index; derivative order 0 through 8; L^(k)(1), computed
  independently with `lfun` for the first two curves.
- `special-functions.csv`: numerator and denominator of x; exp(-x); E1(x),
  with x from 1/1000 through 200, including the interval integrator's cutoff.

Regenerate with a local GP executable (PowerShell, from the repository root):

```powershell
$gpPath = 'C:/path/to/gp.exe'
& $gpPath -q -f tests/Fixtures/generate-conductors.gp |
    Select-Object -Unique | Set-Content -Encoding utf8 tests/Fixtures/conductors.csv
& $gpPath -q -f tests/Fixtures/generate-ranks.gp |
    Select-Object -Unique | Set-Content -Encoding utf8 tests/Fixtures/ranks.csv
& $gpPath -q -f tests/Fixtures/generate-general-descent.gp |
    Where-Object { $_.Trim() } | Select-Object -Unique |
    Set-Content -Encoding utf8 tests/Fixtures/general-descent.csv
& $gpPath -q -f tests/Fixtures/generate-local-roots.gp |
    Select-Object -Unique | Set-Content -Encoding utf8 tests/Fixtures/local-roots.csv
New-Item -ItemType Directory -Force artifacts/native-validation | Out-Null
& $gpPath -q -f tests/Fixtures/generate-analytic.gp
foreach ($name in 'root-numbers', 'analytic', 'coefficients', 'derivatives', 'special-functions') {
    Copy-Item -LiteralPath "artifacts/native-validation/$name.csv" -Destination "tests/Fixtures/$name.csv"
}
dotnet test EllipticCurves.sln
```

The scripts fix the PARI random seed. A different PARI version may change its
random stream or certified rank bounds, so review regenerated fixture changes.

## Local invariants, heights and periods

`generate-extended.gp` uses seed 20260910, 70-digit real precision, random general
equations, weighted families at 2, 3, 5 and 7, short equations and fixed examples.
All three files are deduplicated:

- `local-data.csv`: input a1,a2,a3,a4,a6; prime; minimal discriminant valuation;
  conductor valuation; PARI Kodaira code; Tamagawa number; local root number.
  The wild cases include II, III, IV, I0*, II*, III*, IV* and long I_n* refinements.
- `periods.csv`: input a1,a2,a3,a4,a6; primitive positive real period; absolute
  imaginary part of the second period; BSD real period; period-lattice area.
  Values refer to `ellminimalmodel` and cover both signs of the discriminant.
- `heights.csv`: minimal a1,a2,a3,a4,a6; x; y; canonical height from `ellheight`.
  Points include torsion, integral coordinates and rational denominators.

```powershell
$rows = & $gpPath -q -f tests/Fixtures/generate-extended.gp
foreach ($entry in @(@('L,','local-data'), @('P,','periods'), @('H,','heights'))) {
    $rows | Where-Object { $_.StartsWith($entry[0]) } |
        ForEach-Object { $_.Substring(2) } | Select-Object -Unique |
        Set-Content -Encoding utf8 (Join-Path tests/Fixtures ($entry[1] + '.csv'))
}
```

`lmfdb-*.json` are unmodified aggregate JSON responses downloaded on 2026-09-10
from `https://www.lmfdb.org/EllipticCurve/Q/data/{label}?_format=json`, for labels
11.a1, 37.a1, 48.a3, 389.a1 and 5077.a1. Their original table names, counts and
decimal precision metadata are preserved. Sources:
[11.a1](https://www.lmfdb.org/EllipticCurve/Q/data/11.a1),
[37.a1](https://www.lmfdb.org/EllipticCurve/Q/data/37.a1),
[48.a3](https://www.lmfdb.org/EllipticCurve/Q/data/48.a3),
[389.a1](https://www.lmfdb.org/EllipticCurve/Q/data/389.a1),
[5077.a1](https://www.lmfdb.org/EllipticCurve/Q/data/5077.a1).
LMFDB numerical values are approximations, not certified error intervals; tests
compare them with an explicit tolerance. HTTP behavior tests use an in-memory
handler, including pagination, failure and cancellation, and never call the site.

## Rational-kernel isogenies and numerical real logarithms

`generate-basic-extensions.gp` uses seed 20260911 and 70-digit real precision.
It calls `ellisogeny`, `ellorder` and `ellpointtoz` on explicit equations, without
a database. These are reference outputs only, not translated PARI implementation code.

- `isogenies.csv`: source a1,a2,a3,a4,a6; rational kernel generator x,y; degree;
  target c4,c6. Degrees 2,3,4,5,6 and general source equations are represented.
- `real-logarithms.csv`: minimal a1,a2,a3,a4,a6; point x,y; real logarithm modulo
  the least positive real period; absolute imaginary part; least positive real
  period. Covers both real components, 2-torsion, both ordinate signs and multiples
  with rational denominators. The imaginary sign is normalized to the library's
  positive-imaginary period basis.

```powershell
$rows = & $gpPath -q -f tests/Fixtures/generate-basic-extensions.gp
$rows | Where-Object { $_.StartsWith('I,') } |
    ForEach-Object { $_.Substring(2) } |
    Set-Content -Encoding utf8 tests/Fixtures/isogenies.csv
$rows | Where-Object { $_.StartsWith('L,') } |
    ForEach-Object { $_.Substring(2) } | Select-Object -Unique |
    Set-Content -Encoding utf8 tests/Fixtures/real-logarithms.csv
```

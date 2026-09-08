# Independent arithmetic fixtures

Generated using PARI/GP 2.19.0 development 31236-59c418aa0c, compiled 2026-09-08,
from the [official Windows snapshot](https://pari.math.u-bordeaux.fr/download.html).
All inputs are explicit coefficients. The scripts use `ellglobalred`,
`ellminimalmodel`, `elllocalred`, `ellrank`, `ellrootno`, `ellan`, `ellanalyticrank`
and `lfun`; they do not load an elliptic-curve
database. See the [PARI elliptic-curve reference](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).

The committed CSV files contain unique rows encoded as UTF-8. Tests only read
these files; PARI/GP is not required to build or run the tests.

- `conductors.csv`: input a1,a2,a3,a4,a6; conductor; reduced minimal a1,a2,a3,a4,a6;
  PARI Kodaira codes at 2 and 3. The last two columns record coverage metadata.
- `ranks.csv`: a,b,proved rank for y²=x³+a*x²+b*x. Only cases where PARI's lower
  and upper rank bounds agree are included.
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

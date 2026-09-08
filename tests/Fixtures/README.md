# Independent arithmetic fixtures

Generated using PARI/GP 2.19.0 development 31236-59c418aa0c, compiled 2026-09-08,
from the [official Windows snapshot](https://pari.math.u-bordeaux.fr/download.html).
All inputs are explicit coefficients. The scripts use `ellglobalred`,
`ellminimalmodel`, `elllocalred` and `ellrank`; they do not load an elliptic-curve
database. See the [PARI elliptic-curve reference](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).

The committed CSV files contain unique rows encoded as UTF-8. Tests only read
these files; PARI/GP is not required to build or run the tests.

- `conductors.csv`: input a1,a2,a3,a4,a6; conductor; reduced minimal a1,a2,a3,a4,a6;
  PARI Kodaira codes at 2 and 3. The last two columns record coverage metadata.
- `ranks.csv`: a,b,proved rank for y²=x³+a*x²+b*x. Only cases where PARI's lower
  and upper rank bounds agree are included.

Regenerate with a local GP executable (PowerShell, from the repository root):

```powershell
$gpPath = 'C:/path/to/gp.exe'
& $gpPath -q -f tests/Fixtures/generate-conductors.gp |
    Select-Object -Unique | Set-Content -Encoding utf8 tests/Fixtures/conductors.csv
& $gpPath -q -f tests/Fixtures/generate-ranks.gp |
    Select-Object -Unique | Set-Content -Encoding utf8 tests/Fixtures/ranks.csv
dotnet test sources/EllipticCurves.sln
```

The scripts fix the PARI random seed. A different PARI version may change its
random stream or certified rank bounds, so review regenerated fixture changes.

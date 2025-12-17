<p align="center"><img width="25%" src="docs/ec_logo.png" /></p>
<p align="center"> .NET library for studying and working with elliptic curves </p>    

# About
**EllipticCurves** is a small C# library for studying and working with elliptic curves. It provides functionality to compute and explore:
* coefficients and group structure,  
* discriminant and j-invariant,  
* torsion/rational/integral points,  
* short and minimal Weierstrass model,  
* algebraic/analytic ranks,  
* LMFDB label/url,  
* conductor, etc.  

# Version
You can build **EllipticCurves** from sources or install to your own project using nuget package manager.
| Assembly | Specification | OS | Download | Package |
|-------------|:-------------:|:-------------:|:--------------:|:--------------:|
| [EllipticCurves](sources) | .NET Standard 2.0 | Cross-platform | [Release](https://github.com/asiryan/EllipticCurves/releases/) | [NuGet](https://www.nuget.org/packages/EllipticCurves/) | 

# Installation
C# interface  
```c#
using EllipticCurves;
```
To get started with **EllipticCurves** it is recommended to take a look at the [example project](examples).  
Here are some results for the [elliptic curve](https://arxiv.org/abs/2510.11768): **Y^2 = X^3 - 17X^2 + 72X**.
```
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
Minimal Weirstrass model: y^2 = x^3 + x^2 - 24*x + 36
Torsion: Z/2Z x Z/4Z
Rank(E) = 0
Analytic rank(E) = 0
Cond(E) = 48
Isomorphic to E: True
```

# License
**MIT**  

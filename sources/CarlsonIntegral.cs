using System;
using System.Threading;

namespace EllipticCurves
{
    // Real Carlson R_F. Duplication plus the degree-seven expansion, DLMF 19.26.18 and 19.36.1.
    // Numerical double arithmetic only; this routine supplies no interval certificate.
    internal static class CarlsonIntegral
    {
        internal static double Rf(double x, double y, double z, int maxIterations, CancellationToken token)
        {
            double scale = Math.Max(x, Math.Max(y, z));
            if (!(x >= 0 && y >= 0 && z >= 0 && scale > 0) || double.IsInfinity(scale) ||
                (x == 0 && y == 0) || (x == 0 && z == 0) || (y == 0 && z == 0))
                throw new ArithmeticException("Elliptic-integral arguments are outside the supported numerical range.");
            x /= scale; y /= scale; z /= scale;
            for (int i = 0; i < maxIterations; i++)
            {
                token.ThrowIfCancellationRequested();
                double mean = (x + y + z) / 3;
                double dx = (mean - x) / mean, dy = (mean - y) / mean, dz = (mean - z) / mean;
                if (Math.Max(Math.Abs(dx), Math.Max(Math.Abs(dy), Math.Abs(dz))) < 0.001)
                {
                    double e2 = dx * dy + dy * dz + dz * dx, e3 = dx * dy * dz;
                    double polynomial = 1 - e2 / 10 + e3 / 14 + e2 * e2 / 24 - 3 * e2 * e3 / 44
                        - 5 * e2 * e2 * e2 / 208 + 3 * e3 * e3 / 104 + e2 * e2 * e3 / 16;
                    return polynomial / Math.Sqrt(mean) / Math.Sqrt(scale);
                }
                double sx = Math.Sqrt(x), sy = Math.Sqrt(y), sz = Math.Sqrt(z);
                double lambda = sx * sy + sy * sz + sz * sx;
                x = (x + lambda) / 4; y = (y + lambda) / 4; z = (z + lambda) / 4;
            }
            throw new ArithmeticException("Elliptic-integral iteration limit reached.");
        }
    }
}

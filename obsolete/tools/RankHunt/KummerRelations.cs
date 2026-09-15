using System.Numerics;
using EllipticCurves;

// A zero character vector suggests a combination worth testing for divisibility.
// It never, by itself, asserts divisibility or rational dependence.
static class KummerRelations
{
    public static int[][] Find(EllipticCurveQ curve, EllipticCurvePoint[] points, int primeBound)
    {
        if (new[] { curve.A1, curve.A2, curve.A3, curve.A4, curve.A6 }.Any(a => a.Den != 1))
            throw new ArgumentException("An integral Weierstrass model is required for this experimental relation search.");
        if (points.Any(p => p.IsInfinity || !curve.IsOnCurve(p))) throw new ArgumentException("Invalid affine points");
        var images = new BigInteger[points.Length];
        int bit = 0;
        foreach (int p in Local302.Primes(primeBound))
        {
            if (curve.Discriminant.Num % p == 0) continue;
            int a = Mod(curve.B2.Num, p), b = Mod(8 * curve.B4.Num, p), c = Mod(16 * curve.B6.Num, p);
            for (int root = 0; root < p; root++)
            {
                if ((((long)root + a) * root * root + (long)b * root + c) % p != 0) continue;
                int derivative = (int)((3L * root * root + 2L * a * root + b) % p);
                for (int j = 0; j < points.Length; j++)
                {
                    var x = points[j].X;
                    if (x.Den % p == 0) continue; // reduction to infinity
                    var residue = Mod(4 * x.Num * BigInteger.ModPow(x.Den % p, p - 2, p) - root, p);
                    if (residue == 0) residue = derivative;
                    if (BigInteger.ModPow(residue, (p - 1) / 2, p) == p - 1) images[j] |= BigInteger.One << bit;
                }
                bit++;
            }
        }
        var basis = new Dictionary<int, (BigInteger Image, BigInteger Combination)>();
        var relations = new List<int[]>();
        for (int j = 0; j < images.Length; j++)
        {
            var image = images[j]; var combination = BigInteger.One << j;
            while (!image.IsZero)
            {
                int pivot = checked((int)image.GetBitLength() - 1);
                if (!basis.TryGetValue(pivot, out var row)) { basis[pivot] = (image, combination); break; }
                image ^= row.Image; combination ^= row.Combination;
            }
            if (image.IsZero)
                relations.Add(Enumerable.Range(0, points.Length).Where(i => !(combination & (BigInteger.One << i)).IsZero).ToArray());
        }
        return relations.ToArray();
    }
    static int Mod(BigInteger n, int p) => (int)((n % p + p) % p);
}

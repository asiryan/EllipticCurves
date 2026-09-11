using System.Windows.Media;
using System.Windows.Media.Media3D;

namespace EllipticCurves.Explorer.Controls;

/// <summary>A fixed topological embedding of R²/Z²; its radii do not encode the curve's complex structure.</summary>
internal static class TorusMesh
{
    private const double MajorRadius = 2.1;
    private const double MinorRadius = 0.75;
    public static Point3D Position(double u, double v, double lift = 0)
    {
        var a = 2 * Math.PI * u;
        var b = 2 * Math.PI * v;
        var radius = MinorRadius + lift;
        return new Point3D((MajorRadius + radius * Math.Cos(b)) * Math.Cos(a),
            radius * Math.Sin(b), (MajorRadius + radius * Math.Cos(b)) * Math.Sin(a));
    }

    private static Vector3D Normal(double u, double v) => new(Math.Cos(2 * Math.PI * v) * Math.Cos(2 * Math.PI * u),
        Math.Sin(2 * Math.PI * v), Math.Cos(2 * Math.PI * v) * Math.Sin(2 * Math.PI * u));

    public static MeshGeometry3D Surface() => Patch(96, 48, (u, v) => (Position(u, v), Normal(u, v)));

    public static MeshGeometry3D Cycle(bool firstPeriod, double coordinate, double halfWidth) =>
        Patch(firstPeriod ? 128 : 2, firstPeriod ? 2 : 96, (u, v) =>
        {
            if (firstPeriod) v = coordinate + (v * 2 - 1) * halfWidth;
            else u = coordinate + (u * 2 - 1) * halfWidth;
            return (Position(u, v, 0.012), Normal(u, v));
        });

    public static MeshGeometry3D Sphere() => Patch(16, 12, (u, v) =>
    {
        var a = 2 * Math.PI * u;
        var b = Math.PI * (v - 0.5);
        var normal = new Vector3D(Math.Cos(b) * Math.Cos(a), Math.Sin(b), Math.Cos(b) * Math.Sin(a));
        return (new Point3D(normal.X, normal.Y, normal.Z), normal);
    });

    private static MeshGeometry3D Patch(int columns, int rows, Func<double, double, (Point3D Position, Vector3D Normal)> sample)
    {
        var mesh = new MeshGeometry3D();
        for (var i = 0; i <= columns; i++)
            for (var j = 0; j <= rows; j++)
            {
                var vertex = sample((double)i / columns, (double)j / rows);
                mesh.Positions.Add(vertex.Position);
                mesh.Normals.Add(vertex.Normal);
            }
        for (var i = 0; i < columns; i++)
            for (var j = 0; j < rows; j++)
            {
                var a = i * (rows + 1) + j;
                var b = a + rows + 1;
                foreach (var index in new[] { a, a + 1, b, a + 1, b + 1, b }) mesh.TriangleIndices.Add(index);
            }
        mesh.Freeze();
        return mesh;
    }

    public static Material Material(string hex, bool emissive = false)
    {
        var brush = new SolidColorBrush((Color)ColorConverter.ConvertFromString(hex));
        Material material;
        if (emissive)
        {
            // Emission alone adds to the already rendered surface color. An opaque
            // black base makes cycles and markers match their legend colors.
            var unlit = new MaterialGroup();
            unlit.Children.Add(new DiffuseMaterial(Brushes.Black));
            unlit.Children.Add(new EmissiveMaterial(brush));
            material = unlit;
        }
        else material = new DiffuseMaterial(brush);
        material.Freeze();
        return material;
    }
}

using System.Runtime.InteropServices;

namespace EllipticCurves.Explorer.Windowing;

[StructLayout(LayoutKind.Sequential)]
internal struct NativeRect
{
    public int Left;
    public int Top;
    public int Right;
    public int Bottom;
}

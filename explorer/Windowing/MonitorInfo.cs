using System.Runtime.InteropServices;

namespace EllipticCurves.Explorer.Windowing;

[StructLayout(LayoutKind.Sequential)]
internal struct MonitorInfo
{
    public int Size;
    public NativeRect Monitor;
    public NativeRect WorkArea;
    public uint Flags;
}

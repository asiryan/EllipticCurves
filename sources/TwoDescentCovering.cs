namespace EllipticCurves
{
    internal sealed class TwoDescentCovering
    {
        internal readonly BinaryQuartic Form;
        internal readonly int[] Signature;
        internal bool HasPoint;
        internal TwoDescentCovering(BinaryQuartic form, int[] signature) { Form = form; Signature = signature; }
    }
}

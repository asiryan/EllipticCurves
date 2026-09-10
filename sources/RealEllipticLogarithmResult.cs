namespace EllipticCurves
{
    /// <summary>A numerical elliptic logarithm of a real point, modulo the periods of the minimal model.
    /// Double values are approximations, not certified intervals.</summary>
    public sealed class RealEllipticLogarithmResult
    {
        /// <summary>Minimal model and differential used for the logarithm.</summary>
        public EllipticCurveQ MinimalModel { get; }
        /// <summary>Real part in [0, PrimitiveRealPeriod), subject to floating-point rounding.</summary>
        public double RealPart { get; }
        /// <summary>Zero on the identity component, half the imaginary period on the other component.</summary>
        public double ImaginaryPart { get; }
        /// <summary>Zero for the identity component, one for the bounded real component.</summary>
        public int ComponentIndex { get; }
        /// <summary>Least positive real period of the minimal differential.</summary>
        public double PrimitiveRealPeriod { get; }
        internal RealEllipticLogarithmResult(EllipticCurveQ model, double real, double imaginary, int component, double period)
        { MinimalModel = model; RealPart = real; ImaginaryPart = imaginary; ComponentIndex = component; PrimitiveRealPeriod = period; }
    }
}

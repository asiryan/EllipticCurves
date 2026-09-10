namespace EllipticCurves
{
    /// <summary>Certified periods of the Neron differential on a reduced global minimal model.</summary>
    public sealed class PeriodResult
    {
        /// <summary>Model to which the differential and all returned periods refer.</summary>
        public EllipticCurveQ MinimalModel { get; }
        /// <summary>Smallest positive real period omega1.</summary>
        public RealEnclosure PrimitiveRealPeriod { get; }
        /// <summary>Real part of omega2; zero for a rectangular lattice, omega1/2 otherwise.</summary>
        public RealEnclosure SecondPeriodRealPart { get; }
        /// <summary>Positive imaginary part of omega2. The basis has Im(omega2/omega1) > 0.</summary>
        public RealEnclosure SecondPeriodImaginaryPart { get; }
        /// <summary>Integral over E(R): number of real components times omega1, as in LMFDB and BSD.</summary>
        public RealEnclosure RealPeriod { get; }
        /// <summary>Area of a fundamental parallelogram of the Neron lattice.</summary>
        public RealEnclosure Area { get; }
        internal PeriodResult(EllipticCurveQ e, RealEnclosure w, RealEnclosure re, RealEnclosure im, RealEnclosure real, RealEnclosure area)
        { MinimalModel = e; PrimitiveRealPeriod = w; SecondPeriodRealPart = re; SecondPeriodImaginaryPart = im; RealPeriod = real; Area = area; }
    }
}

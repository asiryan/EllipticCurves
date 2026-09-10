using System;
using System.Numerics;

namespace EllipticCurves
{
    public sealed partial class EllipticCurveQ
    {
        /// <summary>Whether the curve has complex multiplication over the algebraic closure of Q.</summary>
        public bool HasComplexMultiplication => CmDiscriminant != 0;

        /// <summary>Discriminant of the geometric CM endomorphism order; zero for a non-CM curve.
        /// Uses the complete list of thirteen rational CM j-invariants. It does not construct endomorphism maps.</summary>
        public int CmDiscriminant
        {
            get
            {
                if (IsSingular) throw new InvalidOperationException("CM requires a nonsingular curve.");
                var j = JInvariant;
                if (!j.Den.IsOne) return 0;
                var invariants = new long[] { 0, 1728, -3375, 8000, -32768, 54000, 287496,
                    -884736, -12288000, 16581375, -884736000, -147197952000, -262537412640768000 };
                var discriminants = new[] { -3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163 };
                for (int i = 0; i < invariants.Length; i++) if (j.Num == new BigInteger(invariants[i])) return discriminants[i];
                return 0;
            }
        }
    }
}

using System.Globalization;
using System.Numerics;

namespace EllipticCurves.Explorer.Models;

/// <summary>Published polynomial data, in descending powers, from Elkies, arXiv:2608.25406v1, Theorem 4.</summary>
public static class ElkiesSearchFamily
{
    public const string Id = "elkies-17-v1";
    public const string Name = "Elkies family · 17 known sections";
    public const string Source = "https://arxiv.org/html/2608.25406v1#S2.SS1";
    internal static readonly BigInteger[] S = Polynomial("307516108335972163537936 10476571005172375234427296 234256046667228607566274912 2020678721371875903158954848 8789387383568632081365832240 21430123310022469548285709072 28607402618712438778345257832 17860826619093915900857289304 4201305425690184127251888481");
    internal static readonly BigInteger[] T = Polynomial("1050290276365892761266194577222156800 67802587761728815952013525763236564480 2392486076703808362288120169049836903680 38126035250980128714796491999580538771200 372202978476351718721663756748866085220800 2373760737463050257069464720014664373086080 9904246958414858348647761354992989326760320 26905633537996991160744810870319331164617600 47243082583908684509409509915652973906060800 52862444598312784274784438443066814490530880 36435013603665838306995466090052055171475872 13865015501478235534002649882546248548532768 2193201312876924214657300134273061462776968");
    private static readonly BigInteger[][] X = new[]
    {
        "419884536396 -6900780974412 84146613883956 448019664127620 304456582100883",
        "960407733324 32982109084392 228946163155128 451794461129916 136981770876723",
        "191092188813708 1333867432517100 1886874275645632 -1605097821400112 -329325794912045",
        "881331598668 20106018946320 191145949680312 489728900722308 205182206512275",
        "7087168886668 60084894852268 119727776998960 70169296825056 124968924115923",
        "-2654105330292 -75993994150932 -436645278959760 -794565531069024 -379495133581677",
        "-2842297467828 -76265967628812 -435348555495516 -839785632779964 -362275017421677",
        "-2947994548863 -25598671575906 63026792718435 464786649518334 320073994727283",
        "-2414976971316 54450640822344 2834155382615496 -12715093268802228 9732633560757363",
        "13412195434209 -250429886278338 -242466751598877 666777676835166 521473683384723",
        "353434406988 62514191744628 271192708423620 -55497536934924 -125053466701677",
        "-1586365228500 -83483171473260 -486949775116428 -653155402766412 345607319019603",
        "5724934740993 -29562185743194 -764745203764737 -1492345666793982 2685202943830203",
        "2426651051916 18957579322812 164233784236041 694007861500356 964931580370398",
        "3734241561804 27061905787332 153585168336648 413697107315976 260813404752123",
        "-798764561556 -57890934328188 -354534107851872 -678508004328024 -287286844790877",
        "2607059076492 25332675721548 206429765168484 494998206601236 408936541867923"
    }.Select(Polynomial).ToArray();

    public static (EllipticCurveQ Curve, EllipticCurvePoint[] Points) Create(int numerator, int denominator, CancellationToken token = default)
    {
        if (denominator <= 0) throw new ArgumentOutOfRangeException(nameof(denominator));
        var parameter = new BigRational(numerator, denominator);
        var a = parameter.Num; var b = parameter.Den;
        // x'=4*b^4*x, y'=8*b^6*y gives an integral short model without factorization.
        var curve = new EllipticCurveQ(0, 0, 0, new BigRational(-432 * Homogeneous(S, a, b)), new BigRational(432 * Homogeneous(T, a, b)));
        if (curve.IsSingular) return (curve, Array.Empty<EllipticCurvePoint>());
        var points = new EllipticCurvePoint[X.Length];
        for (int i = 0; i < X.Length; i++)
        {
            token.ThrowIfCancellationRequested();
            var x = new BigRational(4 * Homogeneous(X[i], a, b));
            if (!BigRational.IsSquare(x * x * x + curve.A4 * x + curve.A6, out var y))
                throw new ArithmeticException("A published family section did not specialize to a rational point.");
            points[i] = new(x, y);
        }
        return (curve, points);
    }

    private static BigInteger Homogeneous(BigInteger[] coefficients, BigInteger a, BigInteger b)
    {
        BigInteger result = coefficients[0], power = 1;
        for (int i = 1; i < coefficients.Length; i++) { power *= b; result = result * a + coefficients[i] * power; }
        return result;
    }
    private static BigInteger[] Polynomial(string text) => text.Split(' ').Select(value => BigInteger.Parse(value, CultureInfo.InvariantCulture)).ToArray();
}

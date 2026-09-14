using System.Numerics;
using EllipticCurves;

// Independent C# transcription of ICARM #302's published homogeneous family.
static class Family302
{
    public static BigInteger[] Parameters(BigInteger u, BigInteger v)
    {
        var u2=u*u; var uv=u*v; var v2=v*v;
        var u4=u2*u2; var u3v=u2*uv; var u2v2=u2*v2; var uv3=uv*v2; var v4=v2*v2;
        var l=446667*u2+471466*uv+239031*v2;
        var p=318552*u2+368554*uv-72570*v2;
        var q=733413*u2-45082*uv-14960*v2;
        var d=5*(7174492962*u4-7114589515*u3v-22069002960*u2v2+3909144679*uv3-205134150*v4);
        var e=882769396002*u4+811447034567*u3v-1174040743*u2v2-32493137198*uv3-2386325360*v4;
        var b=p*q*(l+p+q)-p*e-q*d;
        return new[]{l,p,q,d,e,b};
    }
    public static EllipticCurveQ Curve(BigInteger u,BigInteger v)
    {
        var p=Parameters(u,v);
        return new(new(-p[0]),new(p[3]+p[4]),new(-p[5]),new(p[3]*p[4]),0);
    }
    public static EllipticCurvePoint[] Points(BigInteger u,BigInteger v)
    {
        var p=Parameters(u,v); var d=p[3];var e=p[4];var pq=p[1]*p[2];
        return new[]{new EllipticCurvePoint(0,0),new(new BigRational(-d),0),
            new(new BigRational(-e),0),new(new BigRational(-pq),new BigRational(-p[1]*(e-pq)))};
    }
}

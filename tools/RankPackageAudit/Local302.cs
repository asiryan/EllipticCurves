// Modular scorer ported for reproducibility on .NET; compared to independent
// project point counts and the archive's C++ output. Scores are NOT rank bounds.
sealed class Local302
{
    public int P {get;}
    readonly int[] chi;
    readonly double[] scores;
    public Local302(int p)
    {
        P=p;chi=Enumerable.Repeat(-1,p).ToArray();chi[0]=0;
        for(int x=1;x<p;x++)chi[Mul(x,x,p)]=1;
        scores=Enumerable.Repeat(double.NaN,p+1).ToArray();
    }
    public static int Mod(long n,int p)=>(int)((n%p+p)%p);
    static int Mul(long x,long y,int p)=>(int)(x*y%p);
    static int Add(int x,int y,int p){int z=x+y;return z>=p?z-p:z;}
    static int Pow(int x,int n,int p){int y=1;for(;n>0;n>>=1,x=Mul(x,x,p))if((n&1)!=0)y=Mul(x,y,p);return y;}
    public int Slot(int u,int v)=>v%P==0?P:Mul(Mod(u,P),Pow(Mod(v,P),P-2,P),P);
    public int Trace(int t)
    {
        int p=P,u=t==p?1:t,v=t==p?0:1;
        int u2=Mul(u,u,p),uv=Mul(u,v,p),v2=Mul(v,v,p);
        int u4=Mul(u2,u2,p),u3v=Mul(u2,uv,p),u2v2=Mul(u2,v2,p),uv3=Mul(uv,v2,p),v4=Mul(v2,v2,p);
        int l=Mod(446667L*u2+471466L*uv+239031L*v2,p);
        int q1=Mod(318552L*u2+368554L*uv-72570L*v2,p);
        int q2=Mod(733413L*u2-45082L*uv-14960L*v2,p);
        int d=Mod(5*(7174492962L*u4-7114589515L*u3v-22069002960L*u2v2+3909144679L*uv3-205134150L*v4),p);
        int e=Mod(882769396002L*u4+811447034567L*u3v-1174040743L*u2v2-32493137198L*uv3-2386325360L*v4,p);
        int b=Mod((long)Mul(q1,q2,p)*Mod(l+q1+q2,p)-(long)q1*e-(long)q2*d,p);
        int b2=Mod((long)l*l+4L*(d+e),p),b4=Mod((long)l*b+2L*d*e,p),b6=Mul(b,b,p);
        int c4=Mod((long)b2*b2-24L*b4,p),c6=Mod(-(long)Mul(b2,b2,p)*b2+36L*b2*b4-216L*b6,p);
        if(Mod((long)Mul(c4,c4,p)*c4-(long)c6*c6,p)==0)return p+1;
        int f=b6,d1=Mod(4L+b2+2L*b4,p),d2=Mod(24L+2L*b2,p),d3=24%p,sum=0;
        for(int x=0;x<p;x++){sum+=chi[f];f=Add(f,d1,p);d1=Add(d1,d2,p);d2=Add(d2,d3,p);}
        if((long)sum*sum>4L*p)throw new Exception("Hasse bound");
        return -sum;
    }
    public double Score(int t)
    {
        if(double.IsNaN(scores[t])){int ap=Trace(t);scores[t]=ap==P+1?0:Math.Log((double)(P+1-ap)/P);}
        return scores[t];
    }
    public static int[] Primes(int limit)
    {
        var composite=new bool[limit+1];var result=new List<int>();
        for(int p=2;p<=limit;p++)if(!composite[p])
        {if(p>=5)result.Add(p);for(long j=(long)p*p;j<=limit;j+=p)composite[j]=true;}
        return result.ToArray();
    }
}

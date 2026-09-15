using System.Diagnostics;
using System.Globalization;
using System.Numerics;

sealed class Scored(int u,int v,double first=0)
{
    public int U=u,V=v;public double First=first,Second,Final;
}
static class SieveReplay
{
    const int Scale=1048576;
    public static (List<Scored> Candidates,long Primitive,double Seconds) Grid(string output)
    {
        const int h=2000,b0=1021,b1=4093,b2=16381,keep=2048,finalKeep=32;
        var watch=Stopwatch.StartNew();var ps=Local302.Primes(b2);
        var tables=ps.Where(p=>p<=b0).Select(p=>new Local302(p)).ToArray();
        foreach(var t in tables)for(int j=0;j<=t.P;j++)t.Score(j);
        Console.WriteLine($"C# initial local tables: {watch.Elapsed.TotalSeconds:F3} s");
        var queue=new PriorityQueue<Scored,double>();var scores=new int[2*h+1];long primitive=0;
        for(int v=1;v<=h;v++)
        {
            Array.Clear(scores);
            foreach(var t in tables)
            {
                int p=t.P;
                if(v%p==0)
                {
                    int value=(int)Math.Round(Scale*t.Score(p),MidpointRounding.AwayFromZero);
                    for(int j=0;j<scores.Length;j++)scores[j]+=value;
                }
                else
                {
                    var pattern=new int[p];
                    int residue=t.Slot(-h,v),inverse=t.Slot(1,v);
                    for(int j=0;j<p;j++)
                    {
                        pattern[j]=(int)Math.Round(Scale*t.Score(residue),MidpointRounding.AwayFromZero);
                        residue+=inverse;if(residue>=p)residue-=p;
                    }
                    for(int off=0;off<scores.Length;off+=p)
                        for(int j=0,n=Math.Min(p,scores.Length-off);j<n;j++)scores[off+j]+=pattern[j];
                }
            }
            for(int j=0;j<scores.Length;j++)
            {
                int u=j-h;if(BigInteger.GreatestCommonDivisor(u,v)!=1)continue;primitive++;
                double s=(double)scores[j]/Scale;
                if(queue.Count<keep)queue.Enqueue(new(u,v,s),s);
                else if(s>queue.Peek().First){queue.Dequeue();queue.Enqueue(new(u,v,s),s);}
            }
        }
        Console.WriteLine($"C# grid: {primitive} primitive parameters, {watch.Elapsed.TotalSeconds:F3} s including tables");
        var candidates=queue.UnorderedItems.Select(x=>x.Element).ToList();var record=new Scored(164518,924945);
        foreach(var c in candidates.Append(record))c.Second=c.First=tables.Sum(t=>t.Score(t.Slot(c.U,c.V)));
        foreach(int p in ps.Where(p=>p>b0&&p<=b1))
        {var t=new Local302(p);foreach(var c in candidates.Append(record))c.Second+=t.Score(t.Slot(c.U,c.V));}
        candidates=candidates.OrderByDescending(c=>c.Second).Take(finalKeep).ToList();
        foreach(var c in candidates.Append(record))c.Final=c.Second;
        Console.WriteLine($"C# refined candidates: {watch.Elapsed.TotalSeconds:F3} s");
        foreach(int p in ps.Where(p=>p>b1))
        {var t=new Local302(p);foreach(var c in candidates.Append(record))c.Final+=t.Score(t.Slot(c.U,c.V));}
        candidates=candidates.OrderByDescending(c=>c.Final).ToList();
        using(var writer=new StreamWriter(output))
        {
            writer.WriteLine("u\tv\tscore_1021\tscore_4093\tscore_16381\trole");
            foreach(var c in candidates.Append(record))writer.WriteLine(FormattableString.Invariant($"{c.U}\t{c.V}\t{c.First:G17}\t{c.Second:G17}\t{c.Final:G17}\t{(ReferenceEquals(c,record)?"published_record_control":"candidate")}"));
        }
        Console.WriteLine($"C# full sieve: {watch.Elapsed.TotalSeconds:F3} s; best={candidates[0].U}/{candidates[0].V}");
        return(candidates,primitive,watch.Elapsed.TotalSeconds);
    }

    public static List<Scored> Rescore(IEnumerable<(int U,int V)> input,int selectionBound,bool holdout,string output)
    {
        var watch=Stopwatch.StartNew();var all=input.Select(x=>new Scored(x.U,x.V)).ToList();
        var record=new Scored(164518,924945);
        foreach(int p in Local302.Primes(65521))
        {
            var t=new Local302(p);
            foreach(var c in all.Append(record))
            {
                double score=t.Score(t.Slot(c.U,c.V));
                if(p<=selectionBound)c.First+=score;
                if(p<=4093)c.Second+=score;
                c.Final+=score;
            }
            if(p==4093&&all.Count>32)all=all.OrderByDescending(c=>c.Second-(holdout?c.First:0)).Take(32).ToList();
        }
        all=all.OrderByDescending(c=>c.Final-(holdout?c.First:0)).ToList();
        using(var writer=new StreamWriter(output))
        {
            writer.WriteLine("u\tv\tselection_score\tfull_score\tholdout_score\trole");
            foreach(var c in all.Append(record))writer.WriteLine(FormattableString.Invariant($"{c.U}\t{c.V}\t{c.First:G17}\t{c.Final:G17}\t{c.Final-c.First:G17}\t{(ReferenceEquals(c,record)?"published_record_control":"candidate")}"));
        }
        Console.WriteLine($"C# rescore: {watch.Elapsed.TotalSeconds:F3} s; best={all[0].U}/{all[0].V}; full={all[0].Final:F9}; holdout={all[0].Final-all[0].First:F9}; record full={record.Final:F9}; record holdout={record.Final-record.First:F9}");
        return all;
    }
}

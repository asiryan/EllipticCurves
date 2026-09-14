using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.Json;
using EllipticCurves;

CultureInfo.CurrentCulture=CultureInfo.InvariantCulture;
var root=Path.GetFullPath("artifacts/rank-package-audit");
var original=Path.Combine(root,"original/elliptic_rank_search");
var output=Path.Combine(root,"reproduced");
Directory.CreateDirectory(output);
var summary=new Dictionary<string,object>();

var record=Load(Path.Combine(original,"record302.json"));
var family=Family302.Curve(164518,924945);
var u=new BigRational(3);
var r=(u*u*record.Curve.B2-family.B2)/12;
var s=(u*record.Curve.A1-family.A1)/2;
var t=(BigRational.Pow(u,3)*record.Curve.A3-family.A3-r*family.A1)/2;
var map=family.ChangeModel(u,r,s,t);
Require(Ainvs(map.Target).SequenceEqual(Ainvs(record.Curve)),"Family isomorphism coefficients");
var mapped=record.Points.Select(map.MapBack).ToArray();
var timer=Stopwatch.StartNew();
var proof=family.GetRankLowerBound(mapped,503);
Require(proof.LowerBound==31,"Project rank 31 certificate");
summary["project_record_certificate"]=new{proof.LowerBound,proof.ImageDimension,proof.NoTwoTorsionPrime,elapsed_ms=timer.Elapsed.TotalMilliseconds};
Console.WriteLine($"Project: family-to-record isomorphism checked, all 31 mapped points certified in {timer.Elapsed.TotalMilliseconds:F3} ms");
foreach(var parameter in new[]{(0,1),(-58,237),(164518,924945)})
{
    var curve=Family302.Curve(parameter.Item1,parameter.Item2);
    var seeds=Family302.Points(parameter.Item1,parameter.Item2);
    Require(seeds.All(curve.IsOnCurve),"Elementary points on curve");
    Require(curve.Add(curve.Add(seeds[0],seeds[1]),seeds[2]).IsInfinity,"Three collinear points sum to zero");
    var cert=curve.GetRankLowerBound(seeds,1009);
    Require(cert.LowerBound==3,"Three elementary independent points");
    Console.WriteLine($"Project elementary sections at {parameter.Item1}/{parameter.Item2}: 4 listed points, lower bound {cert.LowerBound}");
}
var t0=Load(Path.Combine(original,"gp_driver_t0.json"));
var t0Fresh=Load(Path.Combine(output,"gp_t0.json"));
Require(t0.Curve.GetRankLowerBound(t0.Points).LowerBound==12,"Archived T0 certificate");
Require(t0Fresh.Curve.GetRankLowerBound(t0Fresh.Points).LowerBound==12,"Fresh PARI T0 certificate");
summary["project_t0_archived_and_fresh_lower_bound"]=12;
var direct=t0.Curve.SearchRankLowerBound(new()
{NumeratorRadius=10000000000,DenominatorRootBound=1,TimeLimit=TimeSpan.FromSeconds(10),TargetLowerBound=12});
summary["project_t0_equation_only"]=new{direct.LowerBound,point_count=direct.Points.Count,seconds=direct.Elapsed.TotalSeconds,stop=direct.StopReason.ToString()};
Console.WriteLine($"Project T0 equation-only search: lower={direct.LowerBound}, points={direct.Points.Count}, seconds={direct.Elapsed.TotalSeconds:F3}, stop={direct.StopReason}");
File.WriteAllText(Path.Combine(output,"csharp_t0_equation_only.json"),JsonSerializer.Serialize(new
{
    ainvs=Ainvs(t0.Curve).Select(a=>a.ToString()).ToArray(),
    points=direct.Points.Select(p=>new[]{p.X.ToString(),p.Y.ToString()}).ToArray(),
    source="C# SearchRankLowerBound; equation only; no supplied points",
    lower_bound=direct.LowerBound
},new JsonSerializerOptions{WriteIndented=true}));

int checkedRows=0;
foreach(var line in File.ReadLines(Path.Combine(original,"local_check.tsv")).Skip(1))
{
    var v=line.Split('\t').Select(int.Parse).ToArray();int p=v[0],slot=v[1],expected=v[2];
    var curve=slot==p?Family302.Curve(1,0):Family302.Curve(slot,1);
    int fast=new Local302(p).Trace(slot);
    int independent=curve.Discriminant.Num%p==0?p+1:checked((int)(p+1-new EllipticCurveFp(
        p,curve.A1.Num,curve.A2.Num,curve.A3.Num,curve.A4.Num,curve.A6.Num).CountPoints()));
    Require(fast==expected&&independent==expected,$"Local table {p},{slot}");checkedRows++;
}
summary["project_local_tables_checked"]=checkedRows;
Console.WriteLine($"Project and C# scorer agree with all {checkedRows} archived local traces");
if(args.Contains("--sieve"))
{
    var grid=SieveReplay.Grid(Path.Combine(output,"csharp_candidates.tsv"));
    Require(grid.Primitive==4866351,"Primitive parameter count");
    var expected=ReadCandidates(Path.Combine(original,"candidates.tsv"));
    Require(expected.Select(c=>(c.U,c.V)).ToHashSet().SetEquals(grid.Candidates.Select(c=>(c.U,c.V))),"Sieve finalist set");
    double maxError=0;
    foreach(var c in grid.Candidates)
    {
        var other=expected.Single(x=>x.U==c.U&&x.V==c.V);
        maxError=Math.Max(maxError,Math.Abs(c.Final-other.Final));
    }
    Require(maxError<1e-10,"Final score agreement");
    summary["csharp_sieve"]=new{grid.Primitive,grid.Seconds,finalist_count=grid.Candidates.Count,max_score_error=maxError};
    var deep=SieveReplay.Rescore(grid.Candidates.Select(c=>(c.U,c.V)),1021,false,Path.Combine(output,"csharp_deep.tsv"));
    Require(deep[0].U==-58&&deep[0].V==237,"Deep score winner");
    var crt=File.ReadLines(Path.Combine(output,"crt_candidates.tsv")).Skip(1).Select(l=>l.Split('\t'))
        .Select(v=>(int.Parse(v[0]),int.Parse(v[1]))).ToArray();
    var held=SieveReplay.Rescore(crt,37,true,Path.Combine(output,"csharp_crt_rescored.tsv"));
    Require(held[0].U==109321&&held[0].V==104861,"CRT holdout winner");
    summary["deep_best"]=new{deep[0].U,deep[0].V,deep[0].Final};
    summary["crt_holdout_best"]=new{held[0].U,held[0].V,score=held[0].Final-held[0].First};
}
var summaryName=args.Contains("--sieve")?"project_audit.json":"project_core_audit.json";
File.WriteAllText(Path.Combine(output,summaryName),JsonSerializer.Serialize(summary,new JsonSerializerOptions{WriteIndented=true}));
Console.WriteLine("All requested C# audit checks passed.");

static void Require(bool value,string message){if(!value)throw new Exception(message);}
static BigRational Q(string text){var v=text.Split('/');return new(BigInteger.Parse(v[0]),v.Length==1?BigInteger.One:BigInteger.Parse(v[1]));}
static BigRational[] Ainvs(EllipticCurveQ e)=>new[]{e.A1,e.A2,e.A3,e.A4,e.A6};
static (EllipticCurveQ Curve,EllipticCurvePoint[] Points) Load(string path)
{
    using var doc=JsonDocument.Parse(File.ReadAllText(path));
    var a=doc.RootElement.GetProperty("ainvs").EnumerateArray().Select(v=>Q(v.GetString()!)).ToArray();
    var points=doc.RootElement.GetProperty("points").EnumerateArray().Select(v=>new EllipticCurvePoint(Q(v[0].GetString()!),Q(v[1].GetString()!))).ToArray();
    return(new(a[0],a[1],a[2],a[3],a[4]),points);
}
static List<Scored> ReadCandidates(string path)=>File.ReadLines(path).Skip(1)
    .Where(l=>!l.Contains("published_record_control")).Select(l=>l.Split('\t'))
    .Select(v=>new Scored(int.Parse(v[0]),int.Parse(v[1]),double.Parse(v[2])){Second=double.Parse(v[3]),Final=double.Parse(v[4])}).ToList();

"""Equation-only first point via p-adic jets, weighted LLL and degree <= 6 intersections.

Every used auxiliary polynomial satisfies exact jet congruences and S < p^5.
This is a bounded search, not a guarantee of finding a point or computing rank.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'tools/RankHunt'))
from bootstrap import GP, equation, checked_script, save, vec
from point_arithmetic import on_curve, multiply


def model(a):
    a1, a2, a3, a4, a6 = map(int, a)
    if a1 == a3 == 0:
        return (a2, a4, a6), 1
    # X=4*x, Y=8*y+4*a1*x+4*a3. No point or reduction hints needed.
    return (a1*a1+4*a2, 8*(a1*a3+2*a4), 16*(a3*a3+4*a6)), 4


def disc(f):
    a, b, c = f
    return a*a*b*b - 4*b**3 - 4*a**3*c - 27*c*c + 18*a*b*c


def iroot(n, k):
    lo, hi = 0, 1 << ((n.bit_length()+k-1)//k)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**k <= n: lo = mid
        else: hi = mid
    return lo


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s + (s*s != n)


def weights(U, W, V):
    return [V**6, U*V**4, W*V**3, U*U*V*V, U*W*V, U**3]


def boxes(f, numerator=None, denominator=1):
    a, b, c = map(abs, f)
    natural = max(16, a, ceil_sqrt(b), iroot(c, 3)+1)
    base = 1 << (natural-1).bit_length()
    schedule = [(numerator, denominator)] if numerator else [
        (base, 1), (4*base, 1), (base, 2), (16*base, 1), (4*base, 2), (base, 4)]
    result = []
    for U, V in schedule:
        W = ceil_sqrt(U**3+a*U*U*V*V+b*U*V**4+c*V**6)
        # Strict inequality for the sufficient bound; actual S is always checked.
        start = max(3, V+1, iroot(6**6*math.prod(weights(U, W, V)), 15)+1)
        result.append([U, W, V, start])
    return result


GP_CODE = r'''
jetmatrix(ff,pp,xb,yb)={
  my(mm=pp^5,xx,yy,z,inv,vv);
  if(yb,
    z=yb;
    for(k=1,3,z=lift(Mod(z+(subst(ff,x,xb)-z^2)/Mod(2*z,mm),mm)));
    xx=Mod(1,mm)*(xb+t+O(t^5)); yy=Mod(z,mm)+O(t^5); inv=1/Mod(2*z,mm);
    for(k=1,4,yy+=polcoef(subst(ff,x,xx)-yy^2,k)*inv*t^k),
    z=xb;
    for(k=1,3,z=lift(Mod(z-subst(ff,x,z)/Mod(subst(deriv(ff),x,z),mm),mm)));
    xx=Mod(z,mm)+O(t^5); yy=Mod(1,mm)*t+O(t^5); inv=1/Mod(subst(deriv(ff),x,z),mm);
    for(k=1,4,xx+=polcoef(yy^2-subst(ff,x,xx),k)*inv*t^k)
  );
  if(lift(yy^2-subst(ff,x,xx))!=O(t^5),error("Jet identity"));
  vv=[1+O(t^5),xx,yy,xx^2,xx*yy,xx^3];
  matrix(5,6,i,j,lift(polcoef(vv[j],i-1)))
};
auxiliary(ff,pp,xb,yb,ww)={
  my(H=jetmatrix(ff,pp,xb,yb),D=vector(5,i,pp^(6-i))~,B,C,best=0,S,win=0);
  B=matsolvemod(H,D,[0,0,0,0,0]~,1)[2];
  C=B*qflll(matdiagonal(ww)*B);
  for(j=1,6,
    S=sum(i=1,6,ww[i]*abs(C[i,j]));
    if(!win || S<best,win=j;best=S)
  );
  if(best>=pp^5,return([]));
  C=C[,win];
  if(!C,error("Zero auxiliary polynomial"));
  for(i=1,5,if((H*C)[i]%D[i],error("Jet congruence")));
  [Vec(C),best]
};
intersections(ff,cc)={
  my(aa=cc[1]+cc[2]*x+cc[4]*x^2+cc[6]*x^3,bb=cc[3]+cc[5]*x,
     R=aa^2-ff*bb^2,F,xx,yy,out=List());
  if(!R,error("Identically zero intersection"));
  F=factor(R);
  for(i=1,matsize(F)[1],if(poldegree(F[i,1])==1,
    xx=-polcoef(F[i,1],0)/polcoef(F[i,1],1);
    if(subst(bb,x,xx),
      yy=-subst(aa,x,xx)/subst(bb,x,xx);listput(out,[xx,yy]),
      if(!subst(aa,x,xx) && issquare(subst(ff,x,xx),&yy),
        listput(out,[xx,yy]);if(yy,listput(out,[xx,-yy]))
      )
    )
  ));
  Vec(out)
};
inbox(P,U,W,V)={
  my(v);
  if(!issquare(denominator(P[1]),&v) || v>V,return(0));
  denominator(P[2]*v^3)==1 && abs(P[1]*v^2)<=U && abs(P[2]*v^3)<=W
};
runsearch()={
  my(ff=x^3+ffv[1]*x^2+ffv[2]*x+ffv[3],E=ellinit(ainvs),D=poldisc(ff),
     N=0,red=List(),q=3,pp,U,W,V,ww,yy,rr,zz,aux,P,Q,bx,
     branches=0,certified=0,localcount=0,localcert=0,wall=getwalltime());
  while(#red<3,
    if(E.disc%q,N=gcd(N,ellcard(ellinit(ainvs,q)));listput(red,[q,ellcard(ellinit(ainvs,q))]));
    q=nextprime(q+1)
  );
  print("REDUCTIONS ",Vec(red));
  for(job=1,#jobs,
    bx=jobs[job];U=bx[1];W=bx[2];V=bx[3];pp=nextprime(bx[4]);
    while(D%pp==0 || !isprime(pp),pp=nextprime(pp+1));
    ww=[V^6,U*V^4,W*V^3,U^2*V^2,U*W*V,U^3];
    print("BOX ",[job,vector(4,i,Str([U,W,V,pp][i]))]);
    localcount=0;localcert=0;
    for(xb=0,pp-1,
      if(getwalltime()-wall>=budgetms,print("BUDGET ",[branches,certified]);return(0));
      rr=Mod(subst(ff,x,xb),pp);
      if(kronecker(lift(rr),pp)>=0,
        yy=lift(sqrt(rr));zz=if(yy,Set([yy,pp-yy]),[0]);
        for(sign=1,#zz,
          if(getwalltime()-wall>=budgetms,print("BUDGET ",[branches,certified]);return(0));
          branches++;localcount++;
          aux=auxiliary(ff,pp,xb,zz[sign],ww);
          if(#aux,
            certified++;localcert++;
            P=intersections(ff,aux[1]);
            for(k=1,#P,
              if(P[k][2]^2!=subst(ff,x,P[k][1]),error("Intersection recovery"));
              if(inbox(P[k],U,W,V) && lift(Mod(P[k][1],pp))==xb && lift(Mod(P[k][2],pp))==zz[sign],
                Q=if(scale==1,P[k],[P[k][1]/4,(P[k][2]-ainvs[1]*P[k][1]-4*ainvs[3])/8]);
                if(!ellisoncurve(E,Q),error("Inverse model map"));
                if(ellmul(E,Q,N)!=[0],
                  print("HIT ",[vector(2,i,Str(Q[i])),
                    vector(3,i,Str([pp,xb,zz[sign]][i])),vector(3,i,Str([U,W,V][i])),
                    vector(6,i,Str(aux[1][i])),Str(aux[2]),[branches,certified]]);
                  return(1)
                )
              )
            )
          );
          if(branches%32==0,print("PROGRESS ",[branches,certified]))
        )
      )
    );
    print("EXHAUSTED ",[job,localcount,localcert]);
  );
  print("DONE ",[branches,certified]);0
};
'''


def gp_run(script, seconds):
    checked_script(script)
    try:
        p = subprocess.run([str(GP), '-fq', '-s', '64M'], input=script, text=True,
                           capture_output=True, timeout=max(.01, seconds))
        out, err, timed_out = p.stdout, p.stderr, False
        if p.returncode or '***' in err:
            raise ArithmeticError(err[-2000:])
    except subprocess.TimeoutExpired as e:
        out, err, timed_out = e.stdout or b'', e.stderr or b'', True
        if isinstance(out, bytes): out = out.decode(errors='replace')
        if isinstance(err, bytes): err = err.decode(errors='replace')
        if '***' in err: raise ArithmeticError(err[-2000:])
    return out, timed_out


def convolution(a, b, m):
    return [sum(a[j]*b[i-j] for j in range(i+1)) % m for i in range(5)]


def jets(f, p, xb, yb):
    """Independent integer recomputation of all five Taylor congruences."""
    m = p**5
    def value(x): return ((x+f[0])*x+f[1])*x+f[2]
    def derivative(x): return 3*x*x+2*f[0]*x+f[1]
    def series_value(x):
        x2 = convolution(x, x, m); x3 = convolution(x2, x, m)
        return [(x3[i]+f[0]*x2[i]+f[1]*x[i]+(f[2] if i==0 else 0)) % m for i in range(5)]
    if yb:
        z = yb
        for _ in range(3): z = (z+(value(xb)-z*z)*pow(2*z, -1, m)) % m
        x, y = [xb, 1, 0, 0, 0], [z, 0, 0, 0, 0]
        inv = pow(2*z, -1, m)
        for i in range(1, 5): y[i] = (series_value(x)[i]-convolution(y, y, m)[i])*inv % m
    else:
        z = xb
        for _ in range(3): z = (z-value(z)*pow(derivative(z), -1, m)) % m
        x, y = [z, 0, 0, 0, 0], [0, 1, 0, 0, 0]
        inv = pow(derivative(z), -1, m)
        for i in range(1, 5): x[i] = (convolution(y, y, m)[i]-series_value(x)[i])*inv % m
    if convolution(y, y, m) != series_value(x): raise ArithmeticError('Taylor identity')
    x2 = convolution(x, x, m)
    return list(zip([1, 0, 0, 0, 0], x, y, x2, convolution(x, y, m), convolution(x2, x, m)))


def prime(n):
    if n < 2: return False
    if n < 1000000:
        return all(n % d for d in range(2, math.isqrt(n)+1))
    # Exact primality, not a probable-prime label, for large branch moduli.
    out, timeout = gp_run(f'print(isprime({n}));quit;\n', 2)
    return not timeout and out.strip() == '1'


def reduction_certificate(ainvs, point, reductions):
    a = list(map(Q, ainvs)); P = tuple(map(Q, point))
    if not on_curve(a, P): raise ArithmeticError('Point off original equation')
    f, _ = model(a); N = 0; used = set()
    for p, count in reductions:
        if p in used or p < 3 or not prime(p) or disc(f) % p == 0:
            raise ArithmeticError('Invalid good reduction')
        used.add(p)
        actual = 1
        for x in range(p):
            d = ((int(a[0])*x+int(a[2]))**2 + 4*(x**3+int(a[1])*x*x+int(a[3])*x+int(a[4]))) % p
            actual += 1 if d == 0 else 2 if pow(d, (p-1)//2, p) == 1 else 0
        if count != actual: raise ArithmeticError('Incorrect finite group order')
        N = math.gcd(N, actual)
    if len(used) < 2: raise ArithmeticError('Need two distinct good odd primes')
    multiple = multiply(a, P, N)
    if multiple is None: raise ArithmeticError('Infinite order not established')
    return {'method': 'good-reduction torsion bound and exact nonzero multiple',
            'reductions': reductions, 'torsion_annihilator': N,
            'nonzero_multiple': list(map(str, multiple)), 'rank_lower_bound': 1}


def verify_hit(ainvs, hit, reductions):
    point, branch, bounds, coefficients, claimed_S, counts = hit
    p, xb, yb = map(int, branch); U, W, V = map(int, bounds); c = list(map(int, coefficients))
    if not (min(U,W,V)>0 and p>max(2,V) and prime(p) and 0<=xb<p and 0<=yb<p):
        raise ArithmeticError('Invalid prime, box or branch')
    f, scale = model(ainvs)
    if not disc(f) % p or (yb*yb-((xb+f[0])*xb+f[1])*xb-f[2]) % p:
        raise ArithmeticError('Bad reduction or invalid modular point')
    if len(c)!=6 or not any(c): raise ArithmeticError('Zero auxiliary polynomial')
    H = jets(f, p, xb, yb)
    if any(sum(h*z for h,z in zip(row,c)) % p**(5-j) for j,row in enumerate(H)):
        raise ArithmeticError('False jet congruence')
    S = sum(abs(z)*w for z,w in zip(c, weights(U,W,V)))
    if S != int(claimed_S) or S >= p**5: raise ArithmeticError('False smallness certificate')
    a = list(map(Q, ainvs)); x, y = map(Q, point)
    X,Y = (x,y) if scale==1 else (4*x, 8*y+4*a[0]*x+4*a[2])
    v = math.isqrt(X.denominator)
    if v*v!=X.denominator or v>V or (Y*v**3).denominator!=1 or abs(X*v*v)>U or abs(Y*v**3)>W:
        raise ArithmeticError('Point outside certified box')
    if X.numerator*pow(X.denominator,-1,p)%p!=xb or Y.numerator*pow(Y.denominator,-1,p)%p!=yb:
        raise ArithmeticError('Point outside certified branch')
    if c[0]+c[1]*X+c[2]*Y+c[3]*X*X+c[4]*X*Y+c[5]*X**3 != 0:
        raise ArithmeticError('Point not on auxiliary polynomial')
    cert = reduction_certificate(ainvs, point, reductions)
    return {'points': [point], 'certificate': cert, 'auxiliary': {
        'prime': str(p), 'branch': branch[1:], 'bounds_U_W_V': bounds,
        'coefficients_1_x_y_x2_xy_x3': coefficients, 'S': str(S), 'p5': str(p**5)},
        'rank_lower_bound': 1, 'verified': True, 'hit': hit}


def search(data, seconds=3, numerator=None, denominator=1):
    start = time.perf_counter(); data = equation(data); a = list(map(int, data['ainvs']))
    f, scale = model(a); jobs = boxes(f, numerator, denominator)
    script = GP_CODE + '\n' + 'ainvs='+vec(a)+';ffv='+vec(f)+';scale='+str(scale)+';\n'
    script += 'jobs=['+','.join(vec(j) for j in jobs)+'];budgetms='+str(max(1,int(seconds*1000)-30))+';\nrunsearch();quit;\n'
    remaining = seconds-(time.perf_counter()-start)
    out, timed_out = gp_run(script, remaining)
    result = dict(data, method='p-adic jets / certified weighted LLL / polynomial intersection',
                  points=[], status='budget_exhausted' if timed_out else 'no_point_in_tested_boxes',
                  boxes=[], exhausted_boxes=[], branches=0, certified_branches=0,
                  input_points=0, full_rank_calls=False, direct_search_comparison=False)
    reductions = []; hit = None
    for line in out.splitlines():
        if ' ' not in line: continue
        tag, raw = line.split(' ', 1)
        if tag not in ('REDUCTIONS','BOX','EXHAUSTED','PROGRESS','BUDGET','HIT','DONE'): continue
        value = json.loads(raw)
        if tag == 'REDUCTIONS': reductions = value
        elif tag == 'BOX': result['boxes'].append(value)
        elif tag == 'EXHAUSTED': result['exhausted_boxes'].append(value)
        elif tag == 'HIT': hit = value
        elif tag in ('PROGRESS','BUDGET','DONE'):
            result['branches'], result['certified_branches'] = value
            if tag == 'BUDGET': result['status'] = 'budget_exhausted'
    if hit:
        result.update(verify_hit(a, hit, reductions), status='point_found')
        result['branches'], result['certified_branches'] = hit[-1]
    result['total_seconds'] = time.perf_counter()-start
    result['budget_seconds'] = seconds
    result['limitation'] = 'Only tested boxes/branches; failure gives no rank upper bound. Bounds use the integral cubic model.'
    return result


def benchmark(args):
    # Equation fixtures only; this import is in the coordinator, never the search worker.
    from first_cover import CASES
    cases = [('attachment-example', None, [0,0,0,-37,-1686])] + CASES
    if args.cases:
        names = set(args.cases.split(','))
        if names-set(c[0] for c in cases): raise ValueError('Unknown benchmark case')
        cases = [c for c in cases if c[0] in names]
    output = Path(args.output).resolve(); output.parent.mkdir(parents=True, exist_ok=True)
    report = {'method':'p-adic auxiliary polynomial first point', 'direct_search_comparison':False,
              'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'selection':'attachment example plus previously used coefficient-only curve fixtures; not a random sample',
              'trials':[]}
    start = time.perf_counter()
    for repeat in range(args.repeats):
        for name, published, a in cases:
            with tempfile.TemporaryDirectory(prefix='first-jet-', dir=output.parent) as temporary:
                source, dest = Path(temporary)/'input.json', Path(temporary)/'result.json'
                save(source, {'ainvs':list(map(str,a))}); t = time.perf_counter()
                command = [sys.executable, str(Path(__file__)), '--input', str(source), '--output',str(dest),
                           '--seconds',str(args.seconds)]
                p = subprocess.run(command, capture_output=True, text=True, timeout=args.seconds+5)
                trial = {'case':name, 'published_lower_bound':published, 'repeat':repeat+1,
                         'process_seconds':time.perf_counter()-t}
                if p.returncode or not dest.exists(): trial.update(status='worker_error', error=p.stderr[-2000:])
                else:
                    value = json.loads(dest.read_text()); trial.update(value)
                    if value.get('points'): verify_hit(a,value['hit'],value['certificate']['reductions'])
                report['trials'].append(trial); save(output,report)
                print(json.dumps({k:trial.get(k) for k in ('case','repeat','status','process_seconds','points')}),flush=True)
    report.update(complete=True,wall_seconds=time.perf_counter()-start,
                  successes=sum(t.get('status')=='point_found' for t in report['trials']))
    save(output,report)


def self_test():
    f=(0,-37,-1686);c=[16512658,-149240,-6959,-111,107,-17]
    hit=[['70','582'],['53','17','52'],['100','1003','1'],list(map(str,c)),'67258635',[1,1]]
    verify_hit([0,0,0,-37,-1686],hit,[[3,4],[11,14]])
    # Differential jet checks, including y=0, across nonsingular small cubics.
    scripts=[]; expected=[]
    for f in [(0,-37,-1686),(0,-1,0),(2,-3,1)]:
        for p in [3,5,7,11]:
            if disc(f)%p==0: continue
            for xb in range(p):
                for yb in range(p):
                    if (yb*yb-((xb+f[0])*xb+f[1])*xb-f[2])%p: continue
                    expected.append(jets(f,p,xb,yb))
                    scripts.append(f'H=jetmatrix(x^3+({f[0]})*x^2+({f[1]})*x+({f[2]}),{p},{xb},{yb});print(vector(5,i,Vec(H[i,])));')
    out, timeout=gp_run(GP_CODE+'\n'+'\n'.join(scripts)+'\nquit;\n',3)
    got=[json.loads(s) for s in out.splitlines()]
    assert not timeout and len(got)==len(expected)
    assert all([list(row) for row in H]==J for H,J in zip(expected,got))
    # Recover a rational point through its denominator-aware auxiliary polynomial.
    U,W,V=129,383,10; start=iroot(6**6*math.prod(weights(U,W,V)),15)+1
    script=GP_CODE+f'''\np=nextprime({start});X=129/100;Y=383/1000;
aa=auxiliary(x^3-2,p,lift(Mod(X,p)),lift(Mod(Y,p)),{vec(weights(U,W,V))});
print("A ",[Str(p),vector(6,i,Str(aa[1][i])),Str(aa[2])]);
print("P ",vector(#intersections(x^3-2,aa[1]),i,vector(2,j,Str(intersections(x^3-2,aa[1])[i][j]))));quit;\n'''
    out, timeout=gp_run(script,3)
    rows={s.split(' ',1)[0]:json.loads(s.split(' ',1)[1]) for s in out.splitlines()}
    assert not timeout and ['129/100','383/1000'] in rows['P']
    pp,cc,ss=rows['A'];p=int(pp)
    rat_hit=[['129/100','383/1000'],[pp,str(129*pow(100,-1,p)%p),str(383*pow(1000,-1,p)%p)],
             list(map(str,[U,W,V])),cc,ss,[1,1]]
    verify_hit([0,0,0,0,-2],rat_hit,[[5,6],[7,7]])
    # Reject a torsion point even when it lies on a nonsingular curve.
    try: reduction_certificate([0,0,0,-1,0],['0','0'],[[3,4],[5,8]])
    except ArithmeticError: pass
    else: raise AssertionError('Torsion accepted')
    for target in ['point','polynomial','bound']:
        bad=json.loads(json.dumps(hit))
        if target=='point': bad[0][0]='71'
        if target=='polynomial': bad[3][0]=str(int(bad[3][0])+1)
        if target=='bound': bad[4]=str(int(bad[4])+1)
        try: verify_hit([0,0,0,-37,-1686],bad,[[3,4],[11,14]])
        except ArithmeticError: pass
        else: raise AssertionError('Tampered certificate accepted')
    for a in [[0,0,0,-1,1],[0,0,0,-25,0],[1,-1,0,-106384,13075804]]:
        got=search({'ainvs':a},3)
        assert got['status']=='point_found',got
    print('PASS: attachment certificate; GP/integer jets including y=0; rational denominator recovery; torsion/tampering rejection; equation-only end-to-end searches')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input');parser.add_argument('--output')
    parser.add_argument('--seconds',type=float,default=3)
    parser.add_argument('--numerator',type=int);parser.add_argument('--denominator',type=int,default=1)
    parser.add_argument('--benchmark',action='store_true');parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--repeats',type=int,default=2);parser.add_argument('--cases')
    args=parser.parse_args()
    if not (0<args.seconds<=60 and args.denominator>=1 and (args.numerator is None or args.numerator>=1)
            and 1<=args.repeats<=10): parser.error('Invalid bounded settings')
    if args.self_test: self_test()
    elif args.benchmark:
        if not args.output: parser.error('--output required')
        if args.numerator or args.denominator!=1: parser.error('Benchmark uses automatic boxes')
        benchmark(args)
    else:
        if not args.input or not args.output: parser.error('--input and --output required')
        from run import data_boundary
        source,output=Path(args.input).resolve(),Path(args.output).resolve()
        reads=data_boundary(source,output.parent)
        result=search(json.loads(source.read_text(encoding='utf-8-sig')),args.seconds,args.numerator,args.denominator)
        result['unexpected_data_reads']=[p for p in reads if Path(p)!=source and Path(p).suffix not in ('.py','.pyc')]
        save(output,result)
        print(json.dumps({k:result[k] for k in ('status','points','total_seconds')}))

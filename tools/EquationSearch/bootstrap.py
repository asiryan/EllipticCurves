"""Equation-only point discovery. All models and trial points come from coefficients."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

import certificate
from point_arithmetic import on_curve, negate

ROOT = Path(__file__).resolve().parents[2]
GP = ROOT/'artifacts/native-validation/gp.exe'
FORBIDDEN = re.compile(r'\b(?:ellrank|ellrankinit|ell2cover|ellanalyticrank|ellgenerators|ellsearch|ellidentify|read|readstr|system|extern)\s*\(')


def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    temporary.replace(path)


def checked_script(script):
    if FORBIDDEN.search(script): raise ValueError('Forbidden rank, database or external-data operation')
    return script


def vec(items):
    return '['+','.join(str(Q(x)) for x in items)+']'


def poly(items):
    return 'Polrev('+vec(items)+')'


def equation(data):
    if not isinstance(data, dict) or set(data) != {'ainvs'}:
        raise ValueError('Input must contain only ainvs; points and family metadata are forbidden')
    a = list(map(Q, data['ainvs']))
    if len(a) != 5 or any(x.denominator != 1 for x in a):
        raise ValueError('Expected five integral Weierstrass coefficients')
    a1,a2,a3,a4,a6 = a
    b2=a1*a1+4*a2; b4=a1*a3+2*a4; b6=a3*a3+4*a6
    b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    if -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6 == 0:
        raise ValueError('Singular equation')
    return {'ainvs':list(map(str,a))}


def gp(script, seconds):
    checked_script(script); started=time.perf_counter(); timed_out=False
    try:
        p=subprocess.run([str(GP),'-fq','-s','64M'],input=script,text=True,
                         capture_output=True,timeout=max(.05,seconds))
        stdout,stderr,code=p.stdout,p.stderr,p.returncode
    except subprocess.TimeoutExpired as error:
        stdout,stderr,code=error.stdout or '',error.stderr or '',None; timed_out=True
        if isinstance(stdout,bytes): stdout=stdout.decode(errors='replace')
        if isinstance(stderr,bytes): stderr=stderr.decode(errors='replace')
    if code not in (0,None) or any('***' in line and 'Warning:' not in line for line in stderr.splitlines()):
        raise RuntimeError('Point-search arithmetic failed: '+stderr[-3000:])
    return stdout, {'seconds':time.perf_counter()-started,'timed_out':timed_out,
                    'script_sha256':hashlib.sha256(script.encode()).hexdigest()}


def prefix(a):
    digits=max(80, max(len(str(abs(int(x)))) for x in a)+40)
    return ('default(parisizemax,536870912);\n'+f'default(realprecision,{digits});\n'
            +'E=ellinit('+vec(a)+');\n')


def prepare(data, seconds):
    script=prefix(data['ainvs'])
    script+=('M=ellminimalmodel(E,&change);\n'
        'print("MINIMAL ",vector(5,i,Str(M[i])));print("CHANGE ",vector(4,i,Str(change[i])));\n'
        'D=4*x^3+M.b2*x^2+2*M.b4*x+M.b6;'
        'R=concat(polrootsreal(D),polrootsreal(deriv(D)));'
        'print("CENTERS ",vector(#R,i,Str(round(R[i]))));\n'
        'F=factor(M.b6,100000);print("FACTORS ",vector(matsize(F)[1],i,[Str(F[i,1]),F[i,2]]));\n'
        'print("PREPARED");quit;\n')
    out,stats=gp(script,seconds)
    if 'PREPARED' not in out: return None,stats
    values={line.split(' ',1)[0]:json.loads(line.split(' ',1)[1]) for line in out.splitlines() if ' ' in line}
    return {'ainvs':values['MINIMAL'],'change':values['CHANGE'],
            'centers':values['CENTERS'],'partial_b6_factors':values['FACTORS']},stats


def model_script(original, prepared, center, stride):
    a1,a2,a3,a4,a6=prepared['ainvs']
    base=prefix(original['ainvs'])+'change='+vec(prepared['change'])+';\n'
    base+='F='+poly([a6,a4,a2,1])+';Q='+poly([a3,a1])+';\n'
    base+=f'R=[subst(F,x,{center}+{stride}*x),subst(Q,x,{center}+{stride}*x)];\n'
    base+=('C=hyperellred(R,&m);d=m[2][2,1]*x+m[2][2,2];t=(m[2][1,1]*x+m[2][1,2])/d;\n'
        'if(d^4*subst(R[1],x,t)-m[3]^2-d^2*subst(R[2],x,t)*m[3]!=m[1]^2*C[1],error("Model polynomial identity"));\n'
        'if(2*m[1]*m[3]+m[1]*d^2*subst(R[2],x,t)!=m[1]^2*C[2],error("Model linear identity"));\n'
        'print("MODEL ",[vector(5,i,Str(polcoef(C[1],i-1))),vector(3,i,Str(polcoef(C[2],i-1))),'
        '[Str(m[1]),vector(4,i,Str(m[2][(i-1)\\2+1,(i-1)%2+1])),vector(3,i,Str(polcoef(m[3],i-1)))]]);\n'
        'print("MODEL_END");quit;\n')
    return base


def recipes(prepared):
    a=list(map(int,prepared['ainvs']))
    size=max(1,math.isqrt(abs(a[3])))
    # Cubic-root scale using integer arithmetic (no coordinate guess from a point).
    lo,hi=0,1 << ((abs(a[4]).bit_length()+2)//3)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**3<=abs(a[4]): lo=mid
        else: hi=mid
    size=max(size,lo); exponent=max(0,len(str(size))-1)
    scales=sorted({1,10**(exponent//2),10**exponent,10**(exponent+1)})
    centers=list(dict.fromkeys([0]+list(map(int,prepared['centers']))))
    return [(r,s) for s in scales for r in centers]


def canonical_point(a, point):
    p=tuple(map(Q,point))
    if len(p)!=2 or not on_curve(a,p): raise ValueError('Discovered point is off the input curve')
    return tuple(map(str,min(p,negate(a,p))))


def search_script(original, prepared, model, numerator, denominator, stop_first=None):
    f,q,m=model['reduced']; e,matrix,h=m; aa,b,c,d=matrix
    base=prefix(original['ainvs'])+'change='+vec(prepared['change'])+';\n'
    base+='C=['+poly(f)+','+poly(q)+'];\n'
    base+=('emit(t,z)={my(den='+str(c)+'*t+'+str(d)+',u,v,W);if(den==0,return);'
        f'u={model["center"]}+{model["stride"]}*({aa}*t+{b})/den;'
        'v=('+str(e)+'*z+subst('+poly(h)+',x,t))/den^2;'
        'W=ellchangepointinv([u,v],change);if(!ellisoncurve(E,W),error("Bootstrap inverse map"));'
        'print("POINT ",W);};\n')
    # Finite points hidden at infinity by reduction are part of the search too.
    base+=('infq=polcoef(C[2],2);infp=polcoef(C[1],4);'
        f'if({c}!=0 && issquare(infq^2+4*infp,&z),'
        'Z=Set([(-infq+z)/2,(-infq-z)/2]);for(i=1,#Z,'
        f'W=ellchangepointinv([{model["center"]}+{model["stride"]}*({aa})/({c}),'
        f'({e}*Z[i]+({h[2]}))/({c})^2],change);'
        'if(!ellisoncurve(E,W),error("Infinity inverse map"));print("POINT ",W)));\n')
    # If reduction keeps the original point at infinity at infinity, every
    # affine hit is useful: return it immediately, before a batch can time out.
    # Otherwise the first affine hit may map back to infinity, so retain the
    # complete-list mode and explicitly skip the pole of the inverse map.
    if stop_first is None: stop_first=Q(c)==0
    base+=f'H=hyperellratpoints(C,[{numerator},{denominator}],{int(stop_first)});\n'
    base+='for(i=1,#H,emit(H[i][1],H[i][2]));print("SEARCH_END");quit;\n'
    return base


def parse_points(out, a):
    result=[]
    for line in out.splitlines():
        if not line.startswith('POINT '): continue
        match=re.fullmatch(r'POINT \[(-?\d+(?:/\d+)?), (-?\d+(?:/\d+)?)\]',line)
        if match is None: raise ValueError('Invalid point output')
        result.append(canonical_point(a,match.groups()))
    return list(dict.fromkeys(result))


def select_basis(data):
    cert=certificate.build_certificate(data,max_prime=2000)
    if cert['torsion_2_dimension_upper_bound']: return None
    pivots={}; selected=[]
    for j in range(len(data['points'])):
        mask=sum((row['bits'][j]=='1') << i for i,row in enumerate(cert['independent_rows']))
        if certificate.add_row(pivots,mask): selected.append(j)
    basis={'ainvs':data['ainvs'],'points':[data['points'][j] for j in selected]}
    cert=certificate.build_certificate(basis,max_prime=2000)
    claim=certificate.verify_certificate(basis,cert)
    if not claim['all_points_independent_modulo_torsion']: raise ValueError('Invalid bootstrap basis')
    return {**basis,'rank_lower_bound':claim['rank_lower_bound'],'certificate':cert,'verification':claim}


def discover(data, output, seconds=120, workers=12, job_seconds=2, seed_limit=None):
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter(); deadline=start+seconds
    a=list(map(Q,data['ainvs'])); observed={}; jobs=[]; models=[]; errors=[]
    state={'equation_only':True,'status':'preparing','points':[], 'lower_bound':0,
           'workers':workers,'job_seconds':job_seconds,'budget_seconds':seconds,'seed_limit':seed_limit}
    def checkpoint(status):
        state.update(status=status,seconds=time.perf_counter()-start,points=list(observed),
                     jobs=len(jobs),models=len(models),errors=errors)
        save(output/'bootstrap.json',state)
    def keep(points,origin):
        for p in points:
            if p not in observed: observed[p]=origin
    def log(record):
        jobs.append(record)
        with (output/'bootstrap-jobs.jsonl').open('a',encoding='utf-8') as stream:
            stream.write(json.dumps(record)+'\n')
    checkpoint('preparing')
    prepared,stats=prepare(data,min(20,seconds))
    log({'phase':'minimal_model',**stats})
    if prepared is None:
        checkpoint('preparation_timeout'); return None,state
    save(output/'minimal.json',prepared)
    if seed_limit==1:
        from geometry import small_points
        u,r,s,t=map(Q,prepared['change'])
        for (x,y),origin in small_points(prepared['ainvs']):
            point=canonical_point(a,(u*u*x+r,u**3*y+s*u*u*x+t))
            keep([point],{'method':'square_denominator',**origin})
            basis=select_basis({'ainvs':data['ainvs'],'points':[point]})
            if basis and basis['rank_lower_bound']==1:
                state['lower_bound']=1
                save(output/'seed.json',basis)
                save(output/'origins.json',[{'point':point,**observed[point]}])
                checkpoint('seed_found')
                print(json.dumps({'phase':'bootstrap','lower_bound':1,'seconds':state['seconds'],
                                  'input_points':0,'method':'square_denominator'}),flush=True)
                return basis,state
            if time.perf_counter()>=deadline: break
    # An additional coefficient-derived search: divisors of b6, with no assumption
    # that the remaining cofactor is prime. Every generated divisor is exact.
    divisors=[1]
    for p,e in prepared['partial_b6_factors']:
        p=int(p)
        if p<=1: continue
        if len(divisors)*(e+1)>65536: break
        divisors=[d*p**k for d in divisors for k in range(e+1)]
    ma=list(map(Q,prepared['ainvs'])); u,r,s,t=map(Q,prepared['change'])
    divisor_candidates=0
    for d in divisors:
        for sign in (-1,1):
            divisor_candidates+=1
            x=Q(sign*d); b=ma[0]*x+ma[2]
            square=b*b+4*(x**3+ma[1]*x*x+ma[3]*x+ma[4]); z=math.isqrt(max(0,int(square)))
            if square>=0 and z*z==square:
                y=(-b-z)/2; point=(u*u*x+r,u**3*y+s*u*u*x+t)
                keep([canonical_point(a,point)],{'method':'b6_divisor','divisor':str(sign*d)})
        if time.perf_counter()>=deadline: break
    log({'phase':'b6_divisors','candidate_count':divisor_candidates,'points':len(observed)})
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures={}
        for center,stride in recipes(prepared):
            if time.perf_counter()>=deadline: break
            script=model_script(data,prepared,center,stride)
            futures[executor.submit(gp,script,min(job_seconds,deadline-time.perf_counter()))]=(center,stride)
        for future in as_completed(futures):
            center,stride=futures[future]
            try: out,stats=future.result()
            except RuntimeError as error:
                errors.append(str(error)); continue
            record={'phase':'model','center':str(center),'stride':str(stride),**stats}; log(record)
            if 'MODEL_END' not in out: continue
            reduced=json.loads(next(line[6:] for line in out.splitlines() if line.startswith('MODEL ')))
            model={'id':len(models),'center':str(center),'stride':str(stride),'reduced':reduced}
            if any(m['reduced'][:2]==reduced[:2] for m in models): continue
            models.append(model)
        models.sort(key=lambda m:(int(m['stride']),int(m['center'])))
        for i,m in enumerate(models): m['id']=i
        save(output/'models.json',models)
        for numerator,denominator in ((2**16,1),(2**20,16),(2**24,256),(2**28,4096),
                                      (2**32,65536),(2**36,1048576)):
            if time.perf_counter()>=deadline: break
            pending=list(models)
            for offset in range(0,len(pending),workers):
                remaining=deadline-time.perf_counter()
                if remaining<=.05: break
                group=pending[offset:offset+workers]
                futures={executor.submit(gp,search_script(data,prepared,m,numerator,denominator),
                                         min(job_seconds,remaining)):m for m in group}
                for future in as_completed(futures):
                    model=futures[future]
                    try: out,stats=future.result()
                    except RuntimeError as error:
                        errors.append(str(error)); continue
                    points=parse_points(out,a)
                    origin={'model':model['id'],'numerator':numerator,'denominator':denominator}
                    keep(points,origin)
                    log({'phase':'points',**origin,**stats,'finished':'SEARCH_END' in out,'point_count':len(points)})
                checkpoint('searching')
                if observed:
                    basis=select_basis({'ainvs':data['ainvs'],'points':list(observed)})
                    if basis and basis['rank_lower_bound']:
                        if seed_limit is not None and len(basis['points'])>seed_limit:
                            points=sorted(basis['points'],key=lambda p:tuple(
                                max(abs(Q(v).numerator),Q(v).denominator) for v in p))[:seed_limit]
                            basis=select_basis({'ainvs':data['ainvs'],'points':points})
                        state['lower_bound']=basis['rank_lower_bound']
                        save(output/'seed.json',basis)
                        save(output/'origins.json',[{'point':p,**origin} for p,origin in observed.items()])
                        checkpoint('seed_found')
                        print(json.dumps({'phase':'bootstrap','lower_bound':state['lower_bound'],
                                          'seconds':state['seconds'],'input_points':0}),flush=True)
                        return basis,state
            print(json.dumps({'phase':'bootstrap','seconds':round(time.perf_counter()-start,3),
                              'models':len(models),'numerator':numerator,'denominator':denominator,
                              'observed_points':len(observed),'lower_bound':0}),flush=True)
    checkpoint('budget_completed' if time.perf_counter()>=deadline else 'portfolio_completed')
    return None,state

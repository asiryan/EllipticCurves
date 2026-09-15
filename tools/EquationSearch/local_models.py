"""Exact, bounded integral neighbours of generalized quartic search models.

For C: y^2+q(x)y=f(x), a neighbour uses det(M)=e=p in
    old_x=(a*x+b)/(c*x+d), old_y=(e*y+H(x))/(c*x+d)^2.
At odd p, a repeated root r of F=4f+q^2 with p^2|F(r) gives M=[p,r;0,1].
The infinity branch is M=[1,0;0,p]. We retain integral models only, including
the separately enumerated p=2 branches. These changes preserve the quartic
invariants and rational points; finding new independent points still requires
search. The bounded beam is a height heuristic, not a completeness claim.
"""
from fractions import Fraction as Q
import hashlib
import json
import time

from bootstrap import gp, poly
from point_models import complete_records


def _compose(first, second):
    """Compose inverse coordinate maps, with exact rational coefficients."""
    e1,m1,h1=first;e2,m2,h2=second
    e1,e2=Q(e1),Q(e2);a,b,c,d=map(Q,m1);u,v,w,z=map(Q,m2)
    h0,hx,hxx=map(Q,h1);k0,kx,kxx=map(Q,h2)
    h=[e1*k0+h0*z*z+hx*v*z+hxx*v*v,
       e1*kx+2*h0*w*z+hx*(u*z+v*w)+2*hxx*u*v,
       e1*kxx+h0*w*w+hx*u*w+hxx*u*u]
    return [str(e1*e2),list(map(str,[a*u+b*w,a*v+b*z,c*u+d*w,c*v+d*z])),list(map(str,h))]


def _script(model, count, beam, depth, prime_bound):
    script='default(parisizemax,536870912);default(realprecision,100);\n'
    script+='C0=['+poly(model['f'])+','+poly(model['q'])+'];\n'
    script+=r'''
hom(f,m,k)={sum(j=0,k,polcoef(f,j)*(m[1,1]*x+m[1,2])^j*(m[2,1]*x+m[2,2])^(k-j))};
compose(m,n)={[m[1]*n[1],m[2]*n[2],m[1]*n[3]+hom(m[3],n[2],2)]};
size(C)={my(v=concat(Vec(C[1]),Vec(C[2])));vecmax(abs(v))};
check(C,m)={my(qh=hom(C0[2],m[2],2));if(hom(C0[1],m[2],4)-m[3]^2-qh*m[3]!=m[1]^2*C[1],error("Neighbour polynomial identity"));if(2*m[1]*m[3]+m[1]*qh!=m[1]^2*C[2],error("Neighbour linear identity"));1};
seen=Map();mapput(seen,Str(C0),1);front=[[C0,[1,matid(2),0],[]]];best=[];
offer(node,m,p)={my(C=hyperellchangecurve(node[1],m),red,mm,k,row);if(denominator(content(C[1]))!=1||denominator(content(C[2]))!=1,return());C=hyperellred(C,&red);k=Str(C);if(mapisdefined(seen,k),return());mapput(seen,k,1);mm=compose(node[2],compose(m,red));check(C,mm);row=[C,mm,concat(node[3],[p])];listput(nextfront,row);listput(allbest,row)};
'''
    script+=f'for(level=1,{depth},nextfront=List();allbest=List(best);'
    script+=(f'for(ni=1,#front,node=front[ni];C=node[1];F=4*C[1]+C[2]^2;'
        f'forprime(p=2,{prime_bound},'
        'if(p==2,R=[0,1],G=Mod(1,p)*F;R=if(G==0,vector(p,j,j-1),lift(polrootsmod(G))));'
        'for(ri=1,#R,r=R[ri];if(p!=2 && (subst(F,x,r)%(p*p)!=0 || subst(deriv(F),x,r)%p!=0),next);'
        'm=[p,[p,r;0,1],0];if(p==2,for(h=0,1,m[3]=h;offer(node,m,p)),'
        'm[3]=sum(j=0,2,centerlift(Mod(-1/2,p)*polcoef(hom(C[2],m[2],2),j))*x^j);offer(node,m,p)));'
        'if(p==2 || (polcoef(F,4)%(p*p)==0 && polcoef(F,3)%p==0),'
        'm=[p,[1,0;0,p],0];if(p==2,for(h=0,1,m[3]=h*x^2;offer(node,m,p)),'
        'm[3]=sum(j=0,2,centerlift(Mod(-1/2,p)*polcoef(hom(C[2],m[2],2),j))*x^j);offer(node,m,p)))));'
        f'best=vecsort(Vec(allbest),v->size(v[1]));best=best[1..min({count},#best)];'
        'print("NEIGHBOURS ",vector(#best,k,[vector(5,j,Str(polcoef(best[k][1][1],j-1))),'
        'vector(3,j,Str(polcoef(best[k][1][2],j-1))),[Str(best[k][2][1]),'
        'vector(4,j,Str(best[k][2][2][(j-1)\\2+1,(j-1)%2+1])),'
        'vector(3,j,Str(polcoef(best[k][2][3],j-1)))],best[k][3]]));'
        f'front=vecsort(Vec(nextfront),v->size(v[1]));front=front[1..min({beam},#front)];'
        'if(!#front,break));print("NEIGHBOURS_END");quit;\n')
    return script


def _spread_records(records, count):
    """Keep several search scales, not only the smallest final coefficients.

    A coordinate change can lower coefficients while moving a rational point
    from integral x to a large denominator. Intermediate charts therefore
    complement the final reduced charts even when their coefficients are larger.
    Selection uses only equations/maps and the fixed geometric depth schedule.
    """
    if not records:return []
    chosen=[];seen=set()
    def take(rows, level, limit):
        added=0
        for row in rows:
            key=json.dumps(row[:3],separators=(',',':'))
            if key in seen:continue
            seen.add(key);chosen.append((row,level));added+=1
            if added>=limit or len(chosen)>=count:break
    take(records[-1],len(records),min(2,count))
    level=1
    while level<=len(records) and len(chosen)<count:
        take(records[level-1],level,1);level*=2
    if len(chosen)<count:take(records[-1],len(records),count-len(chosen))
    for level,rows in enumerate(records,1):
        if len(chosen)>=count:break
        take(rows,level,count-len(chosen))
    return chosen


def expand_models(models, timeout, cache, count=6, roots=2, beam=6, depth=16, prime_bound=97):
    """Return compatible search models, keeping current anchor metadata on hits.

    The caller supplies its existing ordered models and owns scheduling. Cached
    maps are independent of basis indices. An interrupted expansion retains its
    last complete output and records its achieved depth and time allowance; a
    larger allowance retries it, rather than treating that partial beam as final.
    Neither a rank nor any external points are read.
    """
    if not (0<=count<=64 and 1<=roots<=8 and 1<=beam<=32 and 1<=depth<=32 and 2<=prime_bound<=1009):
        raise ValueError('Invalid neighbour search bounds')
    if not count or timeout<=0:return []
    selected=[];seen=set()
    for model in models:
        if model.get('neighbour_origin') or model['key'] in seen:continue
        seen.add(model['key']);selected.append(model)
        if len(selected)>=roots:break
    deadline=time.perf_counter()+timeout;groups=[]
    arithmetic_cache=cache.setdefault('local_models',{})
    per_root=max(1,(count+len(selected)-1)//max(1,len(selected)))
    for index,model in enumerate(selected):
        key=json.dumps([model['key'],per_root,beam,depth,prime_bound,'stratified-v1'])
        entry=arithmetic_cache.get(key)
        allowance=(deadline-time.perf_counter())/(len(selected)-index)
        # A small tolerance avoids restarting a partial beam merely because the
        # Python bookkeeping took a few microseconds less than on the last call.
        reuse=entry is not None and (entry['complete'] or allowance<=entry['requested_seconds']*1.05)
        if reuse:
            rows=entry['models']
        else:
            if allowance<.03:break
            output,stats=gp(_script(model,per_root,beam,depth,prime_bound),allowance)
            records=list(complete_records(output,'NEIGHBOURS ',stats['timed_out']))
            rows=[]
            if records:
                for (f,q,transform,primes),snapshot in _spread_records(records,per_root):
                    combined=_compose(model['transform'],transform)
                    digest=hashlib.sha256(json.dumps([model['key'],f,q,combined]).encode()).hexdigest()
                    rows.append({'key':digest,'f':f,'q':q,'transform':combined,
                        'coefficient_bits':max(abs(int(v)).bit_length() for v in f+q),
                        'neighbour_origin':model['key'],'neighbour_depth':len(primes),
                        'neighbour_snapshot_level':snapshot,'neighbour_primes':primes})
            completed='NEIGHBOURS_END' in output
            # An unexpectedly slower retry must not erase useful prior maps.
            if entry and not completed and len(records)<entry['completed_levels']:
                rows=entry['models'];levels=entry['completed_levels']
            else:levels=len(records)
            arithmetic_cache[key]={'models':rows,'complete':completed,
                'completed_levels':levels,'requested_seconds':allowance,
                'elapsed_seconds':stats.get('seconds')}
        groups.append([dict(model,**row) for row in rows])
    return [group[i] for i in range(max(map(len,groups),default=0)) for group in groups if i<len(group)][:count]

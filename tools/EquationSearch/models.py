"""Cached pointed models, parity classes and resumable denominator slices."""
import hashlib
import json
from fractions import Fraction as Q

from anchor_diversity import diverse_vectors
from point_search import quartic_reduction_code
from bootstrap import gp, prefix, vec, poly


def parity_vectors(gram, limit, count):
    """Balance H/2H, retaining several small representatives of each class.

    j_P(Q)=-P-Q, and j_(P+2R)(Q-R)=j_P(Q)-R. Equal parity
    therefore gives translation-equivalent projections, not equal search boxes.
    This is a priority heuristic, not a completeness or independence test.
    """
    # Preserve the original candidate pool. Selecting parity representatives
    # before reduction can discard an anchor with much smaller model coefficients.
    rows, stats = diverse_vectors(gram, limit, count)
    groups = {}
    for score, vector in rows:
        key = tuple(v % 2 for v in vector)
        groups.setdefault(key, []).append((score, vector))
    groups = sorted((sorted(group) for group in groups.values()), key=lambda g:g[0])
    selected = []
    for layer in range(max(map(len, groups), default=0)):
        for group in groups:
            if layer < len(group): selected.append(group[layer])
            if len(selected) == count: break
        if len(selected) == count: break
    stats.update(selection='balanced parity classes with short representatives',
                 available_parity_classes=len(groups),
                 selected_parity_classes=len({tuple(v % 2 for v in row) for _,row in selected}))
    return selected, stats


def prepare_models(pool, timeout, cache=None):
    """Reduce each new anchor once, retaining its exact inverse map."""
    cache = {} if cache is None else cache
    a = pool['ainvs']
    def cache_key(p): return json.dumps([a,p], separators=(',',':'))
    missing = [p for p in pool['points'] if cache_key(p) not in cache]
    if missing:
        script = prefix(a)+'E0=E;E=ellminimalmodel(E0,&change);\n'
        script += 'P=['+','.join(vec(p) for p in missing)+'];\n'
        script += ('for(k=1,#P,A=ellchangepoint(P[k],change);x0=A[1];v0=2*A[2]+E.a1*x0+E.a3;'
            'D=x^4-2*(12*x0+E.b2)*x^2+32*v0*x+E.b2^2-8*E.b2*x0-48*x0^2-32*E.b4;'
            'den=denominator(content(D));F=den^2*D;'+quartic_reduction_code(True)+
            'dd=m[2][2,1]*x+m[2][2,2];tt=(m[2][1,1]*x+m[2][1,2])/dd;'
            'if(dd^4*subst(F,x,tt)-m[3]^2!=m[1]^2*C[1],error("Cached polynomial identity"));'
            'if(2*m[1]*m[3]!=m[1]^2*C[2],error("Cached linear identity"));'
            'print("CACHED ",[k,vector(5,j,Str(polcoef(C[1],j-1))),vector(3,j,Str(polcoef(C[2],j-1))),'
            'Str(den),vector(5,j,Str(E[j])),vector(4,j,Str(change[j])),vector(2,j,Str(A[j])),'
            '[Str(m[1]),vector(4,j,Str(m[2][(j-1)\\2+1,(j-1)%2+1])),vector(3,j,Str(polcoef(m[3],j-1)))]]));\n'
            'print("CACHE_END");quit;\n')
        out, _ = gp(script, timeout)
        for line in out.splitlines():
            if not line.startswith('CACHED '): continue
            i,f,q,den,ma,change,anchor,transform = json.loads(line[7:])
            point = missing[i-1]
            key = hashlib.sha256(json.dumps([point,f,q]).encode()).hexdigest()
            cache[cache_key(point)] = {'key':key,'anchor':point,'f':f,'q':q,'den':den,
                'minimal_ainvs':ma,'change':change,'minimal_anchor':anchor,'transform':transform,
                'coefficient_bits':max(abs(int(v)).bit_length() for v in f+q)}
    models = []
    for i,p in enumerate(pool['points']):
        if cache_key(p) not in cache: continue
        model = dict(cache[cache_key(p)],pool_index=i,
                     anchor_height=pool['approximate_heights'][i])
        if 'vectors' in pool:
            model['parity'] = ''.join(str(int(v)%2) for v in pool['vectors'][i])
        models.append(model)
    return models


def parity_order(models):
    groups = {}
    for model in sorted(models,key=lambda m:(m['coefficient_bits'],m['anchor_height'],m['key'])):
        groups.setdefault(model.get('parity',model['key']),[]).append(model)
    ordered = sorted(groups.values(),key=lambda g:(g[0]['coefficient_bits'],g[0]['anchor_height'],g[0]['key']))
    return [group[i] for i in range(max(map(len,ordered),default=0)) for group in ordered if i<len(group)]


def missing_intervals(coverage, n, d):
    """Subtract rectangles already exhausted with a numerator bound >= n."""
    intervals = sorted((max(1,lo),min(d,hi)) for nn,lo,hi in coverage if nn>=n and lo<=d and hi>=1)
    result = []; cursor = 1
    for lo,hi in intervals:
        if lo>cursor: result.append((cursor,lo-1))
        cursor=max(cursor,hi+1)
    if cursor<=d: result.append((cursor,d))
    return result


def record_coverage(coverage, n, lo, hi):
    intervals=sorted((l,h) for nn,l,h in coverage if nn==n)+[(lo,hi)]
    merged=[]
    for l,h in sorted(intervals):
        if merged and l<=merged[-1][1]+1: merged[-1][1]=max(h,merged[-1][1])
        else: merged.append([l,h])
    return [row for row in coverage if row[0]!=n]+[[n,l,h] for l,h in merged]


def search_script(data, models, n, d, coverage=None, slice_width=256):
    """Only completed slices get a marker; an interrupted slice is retried."""
    coverage = coverage or {}
    script=prefix(data['ainvs'])+'E0=E;\n'
    for i,model in enumerate(models,1):
        e,matrix,h = model['transform'];aa,b,c,dd=matrix
        x0,y0 = model['minimal_anchor']
        script += ('E=ellinit('+vec(model['minimal_ainvs'])+');change='+vec(model['change'])+';'
                   'C=['+poly(model['f'])+','+poly(model['q'])+'];\n')
        script += (f'x0={x0};v0=2*({y0})+E.a1*x0+E.a3;\n'
                   'D=x^4-2*(12*x0+E.b2)*x^2+32*v0*x+E.b2^2-8*E.b2*x0-48*x0^2-32*E.b4;\n')
        script += ('emit(t,z)={my(xx,yy,W);if(z^2!=subst(D,x,t),error("Cached quartic inverse"));'
                   'xx=(t^2-E.b2-4*x0+z)/8;yy=(v0+t*(xx-x0)-E.a1*xx-E.a3)/2;'
                   'W=ellchangepointinv([xx,yy],change);if(!ellisoncurve(E0,W),error("Cached curve inverse"));'
                   'print("POINT ",W);};\n')
        script += (f'print("ANCHOR_BEGIN ",{i});\n'
            f'if({c}!=0 && issquare(polcoef(C[2],2)^2+4*polcoef(C[1],4),&zz),'
            'Z=Set([(-polcoef(C[2],2)+zz)/2,(-polcoef(C[2],2)-zz)/2]);'
            f'for(j=1,#Z,emit(({aa})/({c}),(({e})*Z[j]+({h[2]}))/({c})^2/({model["den"]}))));\n')
        for low,high in missing_intervals(coverage.get(model['key'],[]),n,d):
            script += (f'lo={low};while(lo<={high},hi=min({high},min(lo+{slice_width-1},max(lo,4*lo-1)));'
                f'H=hyperellratpoints(C,[{n},[lo,hi]]);'
                'for(j=1,#H,t=H[j][1];z=H[j][2];'
                f'dd=({c})*t+({dd});if(dd==0,next);'
                f'emit((({aa})*t+({b}))/dd,(({e})*z+subst({poly(h)},x,t))/dd^2/({model["den"]})));'
                f'print("SLICE_DONE ",[{i},{n},lo,hi]);lo=hi+1);\n')
        script+=f'print("ANCHOR_DONE ",{i});\n'
    return script+'print("SEARCH_END");quit;\n'

"""Bounded point division and coefficient-derived 2-isogeny quartics.

For a rational 2-torsion root alpha on V^2=f(X), translate U=X-alpha:
    V^2=U*(U^2+A*U+B),  X=4*x, V=4*(2*y+a1*x+a3).
For any nonzero d, the quartic z^2=d*t^4+A*t^2+B/d maps to E by
    U=d*t^2, V=d*t*z.
These are classical 2-isogeny coverings, not a new rank algorithm. We try a
bounded list of square classes obtained from B and already known points.
No Selmer group, upper rank bound, or completeness assertion is computed.
"""
from fractions import Fraction as Q
import hashlib
import json
import time

from bootstrap import gp,prefix,vec
from point_arithmetic import add,multiply,on_curve
from point_search import quartic_reduction_code


def isogeny_data(ainvs, alpha):
    """E -> E' after X=4x, V=4(2y+a1*x+a3), U=X-alpha.

    E: V^2=U(U^2+A*U+B), E': v^2=u(u^2-2*A*u+A^2-4*B).
    The two degree-two maps compose to [2]. They can change search heights,
    but their finite kernels cannot create independent directions by themselves.
    """
    a1,a2,a3,a4,a6=map(Q,ainvs);alpha=Q(alpha)
    b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
    if alpha**3+b2*alpha**2+8*b4*alpha+16*b6:
        raise ValueError('Isogeny kernel is not rational two-torsion')
    A=3*alpha+b2;B=3*alpha**2+2*b2*alpha+8*b4
    if not B or not A*A-4*B:raise ValueError('Singular isogeny model')
    return A,B,list(map(str,(0,-2*A,0,A*A-4*B,0)))


def isogeny_point(ainvs, alpha, point, inverse=False):
    """Exact degree-two map or its dual, including the finite kernels."""
    A,B,dual=isogeny_data(ainvs,alpha);a=list(map(Q,ainvs));alpha=Q(alpha)
    source=list(map(Q,dual)) if inverse else a
    p=None if point is None else tuple(map(Q,point))
    if not on_curve(source,p):raise ValueError('Point off isogeny source')
    if p is None:return None
    if inverse:
        u,v=p
        if not u:return None
        U=(u-2*A+(A*A-4*B)/u)/4
        V=v*(1-(A*A-4*B)/(u*u))/8
        x=(U+alpha)/4;y=(V-4*a[0]*x-4*a[2])/8
        result=(x,y);target=a
    else:
        x,y=p;U=4*x-alpha;V=4*(2*y+a[0]*x+a[2])
        if not U:return None
        result=(U+A+B/U,V*(1-B/(U*U)));target=list(map(Q,dual))
    if not on_curve(target,result):raise ValueError('Isogeny image off target')
    return result


def isogenous_sources(data, timeout, cache):
    """Coefficient-derived neighbours, supplied only images of known points."""
    key=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
    old=cache.get('isogenous_sources',{})
    if old.get('basis_key')==key:return old['sources']
    deadline=time.perf_counter()+timeout
    script=prefix(data['ainvs'])+('F=factor(x^3+E.b2*x^2+8*E.b4*x+16*E.b6);'
        'for(k=1,matsize(F)[1],if(poldegree(F[k,1])==1,'
        'print("ROOT ",Str(-polcoef(F[k,1],0)/polcoef(F[k,1],1)))));print("ROOTS_END");quit;')
    out,stats=gp(script,timeout);sources=[]
    for line in out.splitlines():
        if not line.startswith('ROOT '):continue
        alpha=str(Q(line[5:]));_,_,a=isogeny_data(data['ainvs'],alpha);points=[]
        for p in data['points']:
            image=isogeny_point(data['ainvs'],alpha,p)
            if image is None:continue
            if isogeny_point(data['ainvs'],alpha,image,True)!=multiply(list(map(Q,data['ainvs'])),tuple(map(Q,p)),2):
                raise ValueError('Dual composition is not doubling')
            points.append(list(map(str,image)))
        if not points:continue
        # The divisor beam on a nonminimal isogenous equation can be quite
        # different, even though each resulting quartic is later minimized.
        script=prefix(a)+'M=ellminimalmodel(E,&c);P=['+','.join(vec(p) for p in points)+'];'
        script+=('print("MINIMAL ",[vector(5,j,Str(M[j])),vector(4,j,Str(c[j])),'
                 'vector(#P,k,vector(2,j,Str(ellchangepoint(P[k],c)[j])))]);quit;')
        outmin,st=gp(script,max(.02,deadline-time.perf_counter()))
        records=list(complete_records(outmin,'MINIMAL ',st['timed_out']))
        if not records:continue
        ma,change,mp=records[0];u,r,s,t=map(Q,change)
        if not u or len(points)!=len(mp):raise ValueError('Invalid isogenous model change')
        for old,p in zip(points,mp):
            x,y=map(Q,p)
            if (u*u*x+r,u**3*y+s*u*u*x+t)!=tuple(map(Q,old)) or not on_curve(list(map(Q,ma)),(x,y)):
                raise ValueError('Isogenous minimal model identity')
        sources.append({'ainvs':ma,'points':mp,'transport_alpha':alpha,'transport_change':change})
    if 'ROOTS_END' in out:cache['isogenous_sources']={'basis_key':key,'sources':sources}
    return sources


def complete_records(output, label, timed_out):
    lines=output.splitlines()
    for i,line in enumerate(lines):
        if not line.startswith(label):continue
        try:yield json.loads(line[len(label):])
        except json.JSONDecodeError:
            if timed_out and i==len(lines)-1 and not output.endswith(('\n','\r')):continue
            raise


def divide_basis(data, timeout=.2, depth=3):
    """Replace P by R only after exactly checking m*R=P+T, T finite-order.

This preserves the rational span and does not discover a new free direction.
The bounded loop tests small divisors and torsion shifts without rank calls.
"""
    script=prefix(data['ainvs'])+'P=['+','.join(vec(p) for p in data['points'])+'];'
    script+=('T=elltors(E);G=List([[0]]);for(i=1,#T[2],B=Vec(G);'
             'for(j=1,T[2][i]-1,S=ellmul(E,T[3][i],j);for(k=1,#B,listput(G,elladd(E,B[k],S)))));'
             'N=[2,3,5];'
             f'for(k=1,#P,for(level=1,{depth},ok=0;'
             'for(i=1,#G,for(j=1,#N,S=elladd(E,P[k],G[i]);'
             'if(ellisdivisible(E,S,N[j],&R),m=N[j];'
             'if(R[2]>-R[2]-E.a1*R[1]-E.a3,R=ellneg(E,R);m=-m);'
             'if(ellmul(E,R,m)!=S,error("Point division identity"));'
             'print("DIVIDED ",[k,m,vector(2,l,Str(R[l])),if(#G[i]==1,[],vector(2,l,Str(G[i][l]))),T[1]]);'
             'P[k]=R;ok=1;break));if(ok,break));if(!ok,break)));print("DIVISION_END");quit;\n')
    out,stats=gp(script,timeout)
    a=list(map(Q,data['ainvs']));points=[tuple(map(Q,p)) for p in data['points']];steps=[]
    for k,m,raw,shift,order in complete_records(out,'DIVIDED ',stats['timed_out']):
        r=tuple(map(Q,raw));t=tuple(map(Q,shift)) if shift else None
        if not 1<=k<=len(points) or abs(m) not in (2,3,5) or not 1<=order<=16:
            raise ValueError('Malformed division witness')
        if not on_curve(a,r) or not on_curve(a,t) or multiply(a,t,order) is not None:
            raise ValueError('Invalid division or torsion point')
        if multiply(a,r,m)!=add(a,points[k-1],t):raise ValueError('Incorrect exact division identity')
        steps.append({'index':k-1,'original':list(map(str,points[k-1])),
                      'point':raw,'multiplier':m,'torsion_shift':shift,'torsion_order_multiple':order})
        points[k-1]=r
    return {'ainvs':data['ainvs'],'points':[list(map(str,p)) for p in points]}, {'steps':steps,**stats}


def cover_point(ainvs,alpha,d,t,z):
    """Exact forward map, also used independently by the algebraic tests."""
    a=list(map(Q,ainvs));alpha,d,t,z=map(Q,(alpha,d,t,z))
    b2=a[0]**2+4*a[1];b4=a[0]*a[2]+2*a[3]
    A=3*alpha+b2;B=3*alpha**2+2*b2*alpha+8*b4
    if not d or z*z!=d*t**4+A*t*t+B/d:raise ValueError('Point off isogeny quartic')
    x=(d*t*t+alpha)/4;y=(d*t*z-4*a[0]*x-4*a[2])/8
    if not on_curve(a,(x,y)):raise ValueError('Isogeny inverse missed original curve')
    return x,y


def prepare_covers(data, timeout, cache, count=8):
    """Bounded candidate divisors; local residue tests are necessary filters only."""
    basis_key=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
    old=cache.get('isogeny_covers',{})
    if old.get('basis_key')==basis_key:return old['models']
    candidate_limit=4*count
    script=prefix(data['ainvs'])+'P=['+','.join(vec(p) for p in data['points'])+'];'
    script+=('RF=factor(x^3+E.b2*x^2+8*E.b4*x+16*E.b6);RR=List();'
        'for(j=1,matsize(RF)[1],if(poldegree(RF[j,1])==1,listput(RR,-polcoef(RF[j,1],0)/polcoef(RF[j,1],1))));\n'
        'soluble(F)={my(pp=[3,5,7,11,13,17,19],p,ok=0,Frev=x^4*subst(F,x,1/x));'
        'for(t=0,15,if(setsearch([0,1,4,9],lift(Mod(subst(F,x,t),16))),ok=1;break));'
        'if(!ok,for(t=0,7,if(setsearch([0,1,4,9],lift(Mod(subst(Frev,x,2*t),16))),ok=1;break)));'
        'if(!ok,return(0));for(j=1,#pp,p=pp[j];'
        'ok=issquare(Mod(polcoef(F,4),p));if(!ok,for(t=0,p-1,if(issquare(Mod(subst(F,x,t),p)),ok=1;break)));'
        'if(!ok,return(0)));1};\n'
        'emitted=0;'
        'for(ri=1,#RR,root_emitted=0;alpha=RR[ri];A=3*alpha+E.b2;B=3*alpha^2+2*E.b2*alpha+8*E.b4;'
        'FB=factor(abs(B),10000);pr=Vec(FB[,1]);'
        'K=List([1]);for(k=1,#P,u=4*P[k][1]-alpha;if(u==0,next);'
        'ff=factor(abs(numerator(u)),10000);dd=sign(u)*prod(j=1,matsize(ff)[1],ff[j,1]^(ff[j,2]%2));'
        'if(B%dd==0,listput(K,dd));'
        # Known x-coordinates split unresolved cofactors by exact gcds.
        # This is useful structure from the anchor, not a full factorization.
        'pieces=List();for(j=1,#pr,g=gcd(pr[j],abs(numerator(u)));'
        'if(g>1 && g<pr[j],listput(pieces,g);listput(pieces,pr[j]/g),listput(pieces,pr[j])));pr=Vec(pieces));'
        'L=[1];for(j=1,#pr,L=Set(concat(L,L*pr[j]));L=vecsort(L,d->max(abs(d),abs(B/d)));L=L[1..min(256,#L)]);'
        'D=Set(concat(L,-L));D=select(d->B%d==0,D);D=vecsort(D,d->max(abs(d),abs(B/d)));'
        'D=concat(Vec(K),D);seen=List();'
        'for(di=1,#D,d=D[di];if(setsearch(Set(seen),d)||setsearch(Set(seen),B/d),next);listput(seen,d);'
        'if(d<0 && B/d<0 && (A<=0 || A^2<4*B),next);'
        'F=d*x^4+A*x^2+B/d;if(!soluble(F),next);'
        +quartic_reduction_code(True)+
        'dd=m[2][2,1]*x+m[2][2,2];tt=(m[2][1,1]*x+m[2][1,2])/dd;'
        'if(dd^4*subst(F,x,tt)-m[3]^2!=m[1]^2*C[1],error("Cover polynomial identity"));'
        'if(2*m[1]*m[3]!=m[1]^2*C[2],error("Cover linear identity"));'
        'print("COVER ",[Str(alpha),Str(d),setsearch(Set(Vec(K)),d)!=0,vector(5,j,Str(polcoef(C[1],j-1))),'
        'vector(3,j,Str(polcoef(C[2],j-1))),'
        '[Str(m[1]),vector(4,j,Str(m[2][(j-1)\\2+1,(j-1)%2+1])),vector(3,j,Str(polcoef(m[3],j-1)))]]);'
        f'emitted++;root_emitted++;if(root_emitted>=ceil({candidate_limit}/#RR),break));if(emitted>={candidate_limit},break));'
        'print("COVERS_END");quit;\n')
    # GP uses a single backslash for integer division (kept explicit above).
    script=script.replace('\\\\2','\\2')
    out,stats=gp(script,timeout);models=[]
    for alpha,d,known_class,f,q,transform in complete_records(out,'COVER ',stats['timed_out']):
        a1,_,a3,_,_=map(Q,data['ainvs']);tx=Q(alpha)/4;ty=-(a1*tx+a3)/2
        model={'kind':'isogeny_cover','alpha':alpha,'d':d,'f':f,'q':q,'transform':transform,
               'known_class':bool(known_class),
               'anchor':[str(tx),str(ty)],'anchor_height':0,'pool_index':None,
               'coefficient_bits':max(abs(int(v)).bit_length() for v in f+q)}
        # Coverage belongs to the equation and exact map, not its changing
        # priority or whether its square class has become known meanwhile.
        model['key']=hashlib.sha256(json.dumps([data['ainvs'],alpha,d,f,q,transform],sort_keys=True).encode()).hexdigest()
        models.append(model)
    groups={}
    for model in models:groups.setdefault(model['alpha'],[]).append(model)
    ordered=[]
    for group in groups.values():
        group.sort(key=lambda m:m['coefficient_bits'])
        known=next((m for m in group if m['known_class']),None)
        ordered.append(([known] if known else [])+[m for m in group if m is not known])
    models=[group[i] for i in range(max(map(len,ordered),default=0)) for group in ordered if i<len(group)][:count]
    # A completed or partially successful preparation is reused. An interrupted
    # empty attempt may be retried after a different certified basis is found.
    if models or 'COVERS_END' in out:cache['isogeny_covers']={'basis_key':basis_key,'models':models}
    return models

"""Exact independence certificates with witnessed rational torsion.

Append finite-order points T_1,...,T_t before selecting free columns. If their
character columns have rank t and t equals an upper bound for dim E(Q)[2],
they account for the entire torsion image. Independent remaining columns in
the quotient then prove independence modulo torsion. This includes generators
of order 4 or 8: their doubles alone need not have nonzero Kummer images.

The standalone verifier replays the unchanged v1 checker on the augmented
list and checks m*T=O by exact group arithmetic. It trusts neither PARI's
torsion claim nor any published torsion metadata. No completeness claim about
the supplied finite subgroup is necessary.
"""
from functools import lru_cache
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'RankHunt'))
import certificate as local
from point_arithmetic import add, multiply

FORMAT = 'torsion-augmented-local-2-kummer-v1'
SINGLETON = 'rational-nontorsion-mazur-v1'


@lru_cache(maxsize=64)
def torsion_points(ainvs):
    # Usually a tiny rootless reduction rules out even torsion immediately.
    _, _, cubic, disc = local.validate_curve_and_points({'ainvs':ainvs,'points':[]})
    for p in local.primes_up_to(43):
        if p != 2 and disc % p and not local.roots_mod(cubic,p): return ()
    from bootstrap import gp, prefix
    script = prefix(ainvs) + ('T=elltors(E);for(i=1,#T[2],if(T[2][i]%2==0,'
        'print("TORSION ",[T[2][i],vector(2,j,Str(T[3][i][j]))])));'
        'print("TORSION_END");quit;\n')
    stdout, _ = gp(script, 1)
    if 'TORSION_END' not in stdout: raise ValueError('Torsion witness preparation timed out')
    result = tuple((m,tuple(p)) for m,p in
                   (json.loads(l[8:]) for l in stdout.splitlines() if l.startswith('TORSION ')))
    check_torsion(ainvs,result)
    return result


def check_torsion(ainvs, witnesses):
    data={'ainvs':ainvs,'points':[p for _,p in witnesses]}
    a, points, _, _ = local.validate_curve_and_points(data)
    for (m,_),p in zip(witnesses,points):
        if type(m) is not int or not 2<=m<=12 or multiply(a,p,m) is not None:
            raise ValueError('Invalid finite-order witness')


def singleton_certificate(data):
    """Mazur: every rational torsion point has order at most twelve."""
    a,points,_,_=local.validate_curve_and_points(data)
    if len(points)!=1: raise ValueError('A singleton witness requires exactly one point')
    current=None
    for _ in range(12):
        current=add(a,current,points[0])
        if current is None: raise ValueError('The point has finite order')
    return {'format':SINGLETON,'input_sha256':local.fingerprint(a,points),
            'multiples_checked':12,'rank_lower_bound':1,
            'all_points_independent_modulo_torsion':True,'conditional_assumptions':[]}


def translate_pool(pool, limit):
    """Use torsion cosets, and torsion itself, as additional exact projections."""
    ts=torsion_points(tuple(pool['ainvs']))
    if not ts: return pool
    a=list(map(Q,pool['ainvs']));group={None}
    for m,raw in ts:
        step=tuple(map(Q,raw));old=list(group);shift=None
        for _ in range(1,m):
            shift=add(a,shift,step)
            group.update(add(a,p,shift) for p in old)
    shifts=[None]+sorted(p for p in group if p is not None)
    points=[];heights=[];vectors=[];translations=[];seen=set()
    def insert(p,height,vector,shift):
        if p is None or p in seen or len(points)>=limit:return
        seen.add(p);points.append(list(map(str,p)));heights.append(height)
        vectors.append(vector);translations.append(None if shift is None else list(map(str,shift)))
    for i,raw in enumerate(pool['points']):
        p=tuple(map(Q,raw));v=pool['vectors'][i]
        for shift in shifts: insert(add(a,p,shift),pool['approximate_heights'][i],v,shift)
        if i==0:
            for shift in shifts[1:]: insert(shift,0,[0]*len(v),shift)
    return {**pool,'points':points,'approximate_heights':heights,'vectors':vectors,
            'torsion_translations':translations}


def restrict_certificate(cert, ainvs, points, indices):
    """Restrict an already computed character matrix without rescanning primes."""
    rows=[];basis={};n=len(indices)
    for row in cert['independent_rows']:
        bits=''.join(row['bits'][j] for j in indices)
        mask=sum((bit=='1')<<j for j,bit in enumerate(bits))
        if local.add_row(basis,mask): rows.append({**row,'bits':bits})
    selected=[points[j] for j in indices]
    a,pp=local.normalized_input({'ainvs':ainvs,'points':selected})
    t=cert['torsion_2_dimension_upper_bound']
    return {**cert,'input_sha256':local.fingerprint(a,pp),'points_checked':n,
            'independent_rows':rows,'matrix_rank':len(basis),
            'rank_lower_bound':max(0,len(basis)-t),
            'all_points_independent_modulo_torsion':len(basis)==n and t==0}


def select_basis(data, witnesses=None, max_prime=1009):
    witnesses=torsion_points(tuple(data['ainvs'])) if witnesses is None else witnesses
    check_torsion(data['ainvs'],witnesses)
    tpoints=[list(p) for _,p in witnesses];k=len(tpoints)
    augmented={'ainvs':data['ainvs'],'points':tpoints+data['points']}
    cert=local.build_certificate(augmented,max_prime=max_prime)
    rows=cert['independent_rows'];basis={};torsion_indices=[];selected=[]
    for j in range(len(augmented['points'])):
        column=sum((r['bits'][j]=='1')<<i for i,r in enumerate(rows))
        if j==k and len(basis)!=cert['torsion_2_dimension_upper_bound']: break
        if local.add_row(basis,column):
            (torsion_indices if j<k else selected).append(j)
    free={'ainvs':data['ainvs'],'points':[augmented['points'][j] for j in selected]}
    if not selected: return None
    reduced=restrict_certificate(cert,data['ainvs'],augmented['points'],torsion_indices+selected)
    if not torsion_indices:
        proof=reduced
    else:
        a,pp=local.normalized_input(free)
        proof={'format':FORMAT,'input_sha256':local.fingerprint(a,pp),
               'torsion_witnesses':[[witnesses[j][0],list(witnesses[j][1])] for j in torsion_indices],
               'augmented_certificate':reduced,'rank_lower_bound':len(selected),
               'all_points_independent_modulo_torsion':True,'conditional_assumptions':[]}
    return {**free,'rank_lower_bound':len(selected),'certificate':proof}


def verify_certificate(data, cert):
    if cert.get('format')==SINGLETON:
        if cert!=singleton_certificate(data):raise ValueError('Incorrect singleton certificate')
        return {'points_checked':1,'rank_lower_bound':1,
                'all_points_independent_modulo_torsion':True,'conditional_assumptions':[]}
    if cert.get('format')!=FORMAT: return local.verify_certificate(data,cert)
    a,pp,_,_=local.validate_curve_and_points(data)
    if cert['input_sha256']!=local.fingerprint(a,pp): raise ValueError('Certificate belongs to different inputs')
    witnesses=cert['torsion_witnesses'];check_torsion(a,witnesses)
    augmented={'ainvs':a,'points':[p for _,p in witnesses]+data['points']}
    claim=local.verify_certificate(augmented,cert['augmented_certificate'])
    n=len(pp);t=len(witnesses)
    if claim['torsion_2_dimension_upper_bound']!=t or claim['matrix_rank']!=n+t:
        raise ValueError('Torsion-augmented matrix does not prove independence')
    expected={'rank_lower_bound':n,'all_points_independent_modulo_torsion':True,'conditional_assumptions':[]}
    for key,value in expected.items():
        if cert.get(key)!=value: raise ValueError('Incorrect certificate claim: '+key)
    return {'points_checked':n,'torsion_columns_checked':t,**expected}


def certify(path):
    """Use the fast original selector when torsion correction is unnecessary."""
    from blind_search import certify as original
    data=json.loads(Path(path).read_text())
    torsion=torsion_points(tuple(data['ainvs']))
    if not torsion: return original(path)
    result=select_basis(data,torsion)
    if result is None:
        for p in data['points']:
            try:singleton_certificate({'ainvs':data['ainvs'],'points':[p]})
            except ValueError:continue
            return {'LowerBound':1,'all_selected_independent':True,'points':[p]}
        return {'LowerBound':0,'all_selected_independent':False,'points':[]}
    return {'LowerBound':result['rank_lower_bound'],'all_selected_independent':True,
            'points':result['points']}

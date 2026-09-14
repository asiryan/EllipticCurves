"""Construct curves and points from split polynomials, without record-family data.

P=product(x-r_i), Q=the polynomial part of sqrt(P), F=Q^2-P.
Eight roots give a cubic F; ten give a quartic. This classical identity
construction supplies points, not an independence claim or a route to rank 32.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import time

from point_arithmetic import on_curve

ROOT = Path(__file__).resolve().parents[2]


def verifier():
    path = ROOT/'results/record-hunt-20260914/independent_certificate_verifier.py'
    if hashlib.sha256(path.read_bytes()).hexdigest() != '9e0d0d2562fc53705e92a2eaa9a3f6e7c923f1cd3fd82b14df68b60268f4ad54':
        raise ValueError('Independent verifier changed')
    spec = importlib.util.spec_from_file_location('synthesis_verifier', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def save(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    temporary.replace(path)


def product(a, b):
    result = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): result[i+j] += x*y
    return result


def evaluate(a, x):
    result = 0
    for value in reversed(a): result = result*x+value
    return result


def canonical_roots(values):
    values = sorted(set(map(int, values)))
    if len(values) not in (8, 10): raise ValueError('Need 8 or 10 distinct integer roots')
    values = [r-values[0] for r in values]
    divisor = math.gcd(*values); values = tuple(r//divisor for r in values)
    return min(values, tuple(values[-1]-r for r in reversed(values)))


def identity(values):
    roots = canonical_roots(values); p = [Q(1)]
    for r in roots: p = product(p, [-r, 1])
    degree = len(roots)//2; q = [Q(0)]*degree+[Q(1)]
    for j in range(degree-1, -1, -1):
        q[j] = (p[degree+j]-sum(q[k]*q[degree+j-k]
                 for k in range(j+1, degree)))/2
    scale = math.lcm(*(v.denominator for v in q))
    g = [int(v*scale) for v in q]; square = product(g, g)
    residual = [int(square[i]-scale*scale*p[i]) for i in range(len(p))]
    if any(residual[degree:]): raise ArithmeticError('Polynomial cancellation failed')
    f = residual[:degree]
    while f and f[-1] == 0: f.pop()
    return roots, g, f, scale


def model(f, g):
    """Exact maps from cubic/quartic y^2=F(x); no factorization or rank call."""
    if len(f) == 4:
        d, c, b, a = f
        return [0, b, 0, a*c, a*a*d], lambda x, y: (Q(a*x), Q(a*y))
    if len(f) != 5 or g[0] == 0: raise ValueError('Unusable degree or quartic origin')
    e, d, c, b, a = f; q = g[0]
    if e != q*q: raise ArithmeticError('Quartic origin is not a rational point')
    curve = [0, c, 0, b*d-4*a*e, a*d*d+e*b*b-4*a*e*c]
    def point(x, y):
        x, y = Q(x), Q(y)
        if x == 0:
            if y == q: return None  # the chosen group origin
            if y != -q: raise ValueError('Not a point at the quartic origin')
            X = Q(d*d, 4*e)-c
            return X, Q(-d*X-2*e*b, 2*q)
        X = (2*q*(y+q)+d*x)/(x*x)
        Y = ((X*X-4*a*e)*x-d*X-2*e*b)/(2*q)
        return X, Y
    return curve, point


def normalize(a, points):
    """Remove small common scaling factors; not a claim of global minimality."""
    a = list(a); points = list(points); scale = 1
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        while a[1] % p**2 == a[3] % p**4 == a[4] % p**6 == 0:
            a[1] //= p**2; a[3] //= p**4; a[4] //= p**6; scale *= p
    points = [(x/scale**2, y/scale**3) for x, y in points]
    points = sorted(set((x, -abs(y)) for x, y in points))
    if any(not on_curve(a, p) for p in points): raise ArithmeticError('Point transport failed')
    return {'ainvs':list(map(str, a)), 'points':[[str(x), str(y)] for x, y in points]}, scale


def construct(roots, scan):
    roots, g, f, identity_scale = identity(roots)
    a, mapping = model(f, g)
    # Cubic discriminant: reject singular outputs before any normalization loop.
    _, b, _, c, d = a
    disc = b*b*c*c-4*c**3-4*b**3*d-27*d*d+18*b*c*d
    if disc == 0: raise ValueError('Singular construction')
    seeds = []
    for x in roots:
        y = evaluate(g, x)
        if y*y != evaluate(f, x): raise ArithmeticError('Root witness failed')
        point = mapping(x, y)
        if point is not None: seeds.append(point)
    observed = list(seeds); extra = []
    for x in range(-scan, roots[-1]+scan+1):
        value = evaluate(f, x)
        if value < 0: continue
        y = math.isqrt(value)
        if y*y != value: continue
        for sign in (1, -1):
            point = mapping(x, sign*y)
            if point is not None: observed.append(point)
        if x not in roots: extra.append([str(x), str(y)])
    data, scale = normalize(a, observed)
    data['j_invariant'] = str(Q((16*b*b-48*c)**3, 16*disc))
    seed_data, seed_scale = normalize(a, seeds)
    if scale != seed_scale or data['ainvs'] != seed_data['ainvs']: raise ArithmeticError('Model mismatch')
    data['construction'] = {'method':'split-polynomial identity', 'roots':roots, 'g':list(map(str,g)),
        'f':list(map(str,f)), 'identity_scale':identity_scale, 'model_scale':scale,
        'integer_scan_radius':scan, 'extra_quartic_points':extra,
        'seed_points':seed_data['points'], 'record_curve_or_family_used':False}
    return data


def select_basis(data, prime_bound=251):
    v = verifier(); cert = v.build_certificate(data, max_prime=prime_bound)
    # Conservative: only export a basis when rational 2-torsion is excluded.
    if cert['torsion_2_dimension_upper_bound']:
        cert = v.build_certificate(data, max_prime=max(1009, prime_bound))
        if cert['torsion_2_dimension_upper_bound']: return None
    columns = {}; selected = []
    for j in range(len(data['points'])):
        mask = sum((row['bits'][j] == '1') << i for i, row in enumerate(cert['independent_rows']))
        if v.add_row(columns, mask): selected.append(j)
    basis = {**data, 'points':[data['points'][j] for j in selected]}
    proof = v.build_certificate(basis, max_prime=max(prime_bound, cert['largest_prime_examined']))
    claim = v.verify_certificate(basis, proof)
    if not claim['all_points_independent_modulo_torsion']: raise ArithmeticError('Invalid selected basis')
    return {**basis, 'rank_lower_bound':len(selected), 'certificate':proof, 'verification':claim}


def trial(task):
    index, seed, height, root_count, scan, *parent = task
    rng = random.Random((seed << 32)+index)
    roots = [0]+rng.sample(range(1, height+1), root_count-1)
    if parent:
        roots = list(parent[0]); height = max(height, roots[-1])
        for j in rng.sample(range(1,len(roots)),rng.choice([1,1,2])):
            available = [r for r in range(max(1,roots[j]-8),min(height,roots[j]+8)+1) if r not in roots]
            if not available: available = [r for r in range(1,height+1) if r not in roots]
            if available: roots[j] = rng.choice(available)
    try: data = construct(roots, scan)
    except ValueError as error: return {'status':'unusable', 'reason':str(error)}
    result = select_basis(data)
    if result is None: return {'status':'torsion_not_excluded'}
    v = verifier()
    seed_data = {'ainvs':data['ainvs'], 'points':data['construction']['seed_points']}
    result['construction']['seed_lower_bound'] = v.build_certificate(seed_data, max_prime=251)['rank_lower_bound']
    if parent: result['construction']['mutation_parent_roots'] = parent[0]
    result.update(status='certified', trial=index, sampling_height=height,
        observed_point_count=len(data['points']), coefficient_bits=max(abs(int(a)).bit_length() for a in data['ainvs']))
    return result


def ordering(row):
    return (-row['rank_lower_bound'], row['coefficient_bits'], row['trial'])


def run(args):
    output = Path(args.output).resolve()
    if output == ROOT or not output.is_relative_to(ROOT): raise ValueError('Use a dedicated workspace directory')
    result_path = Path(args.result).resolve()
    if not result_path.is_relative_to(ROOT) or result_path == output/'pool.json': raise ValueError('Invalid result path')
    config = {k:getattr(args,k) for k in ('samples','heights','roots','scan','seed','keep','workers','adaptive')}
    config['code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    pool_path = output/'pool.json'
    if pool_path.exists():
        state = json.loads(pool_path.read_text())
        if state['config'] != config: raise ValueError('Synthesis settings changed; use a new output directory')
    else: state = {'config':config, 'next_trial':0, 'rows':[], 'counts':{}, 'seconds':0.0}
    started = time.perf_counter(); before = state['seconds']
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for start in range(state['next_trial'], args.samples, 128):
            end = min(start+128, args.samples)
            tasks = [(i,args.seed,args.heights[(i//len(args.roots)) % len(args.heights)],
                      args.roots[i % len(args.roots)],args.scan) for i in range(start,end)]
            if args.adaptive:
                for j, task in enumerate(tasks):
                    # Half the draws remain fresh. Mutations use only our own
                    # earlier constructions, and recompute every point from scratch.
                    parents = [r for r in state['rows'] if len(r['construction']['roots']) == task[3]][:8]
                    if parents and task[0] % 4 < 2:
                        ancestor = parents[(task[0]//4) % len(parents)]
                        tasks[j] = (*task, ancestor['construction']['roots'])
            for row in executor.map(trial, tasks, chunksize=4):
                key = row['status']; state['counts'][key] = state['counts'].get(key,0)+1
                if key != 'certified': continue
                if 'mutation_parent_roots' in row['construction']:
                    state['counts']['certified_mutations'] = state['counts'].get('certified_mutations',0)+1
                key = 'bound_'+str(row['rank_lower_bound']); state['counts'][key] = state['counts'].get(key,0)+1
                state['rows'].append(row)
            unique = {}
            for row in sorted(state['rows'], key=ordering): unique.setdefault(row['j_invariant'], row)
            state['rows'] = list(unique.values())[:args.keep]
            state.update(next_trial=end, seconds=before+time.perf_counter()-started)
            save(pool_path,state)
            if end % 1024 == 0 or end == args.samples:
                print(json.dumps({'generated':end,'best_lower_bound':state['rows'][0]['rank_lower_bound'] if state['rows'] else 0,
                                  'seconds':round(state['seconds'],2)}),flush=True)
    if not state['rows']: raise ValueError('No certified nonsingular candidates')
    # Exact certification precedes expensive seeded point search.
    candidates = state['rows'][:args.search_count]
    prepared = {}
    if args.seconds:
        # Resolve and create all paths serially. Concurrent creation of the
        # shared parent caused a Windows resolve() workspace-check failure.
        for row in candidates:
            folder = (output/'search'/f'c{row["trial"]:06d}').resolve()
            source = (output/'inputs'/f'c{row["trial"]:06d}.json').resolve()
            if not folder.is_relative_to(output) or not source.is_relative_to(output):
                raise ValueError('Candidate path escapes synthesis output')
            folder.mkdir(parents=True,exist_ok=True)
            if source.exists():
                if json.loads(source.read_text()) != row: raise ValueError('Candidate input differs from generated curve')
            else: save(source,row)
            prepared[row['trial']] = source,folder
    def deepen(row):
        source,folder = prepared[row['trial']]
        checkpoint = folder/'checkpoint.json'
        used = json.loads(checkpoint.read_text())['wall_seconds'] if checkpoint.exists() else 0
        if used < args.seconds:
            found = search(source,folder,max(1,args.seconds-used),args.point_workers,
                           args.anchors,32,batch_size=4)
        else: found = json.loads((folder/'points.json').read_text())
        return {**row, **independent_result(found,found['rank_lower_bound']),
                'search_checkpoint':str(checkpoint.relative_to(ROOT)),
                'search_seconds':json.loads(checkpoint.read_text())['wall_seconds']}
    best = state['rows'][0]
    if candidates and args.seconds:
        from fast_search import search, independent_result
        with ThreadPoolExecutor(max_workers=args.parallel_curves) as executor:
            futures = {executor.submit(deepen,row):row for row in candidates}
            for future in as_completed(futures):
                try: row = future.result()
                except Exception as error:
                    print(json.dumps({'failed_trial':futures[future]['trial'],'error':str(error)}),flush=True)
                    raise
                print(json.dumps({'finished_trial':row['trial'],'seed_lower_bound':row['construction']['seed_lower_bound'],
                                  'lower_bound':row['rank_lower_bound']}),flush=True)
                if ordering(row) < ordering(best): best = row
    best = {**best, 'target':32, 'target_reached':best['rank_lower_bound']>=32,
            'global_novelty_not_established':True}
    if not result_path.exists() or json.loads(result_path.read_text())['rank_lower_bound'] < best['rank_lower_bound']:
        save(result_path,best)
    print(json.dumps({'status':'completed','best_this_run':best['rank_lower_bound'],
                      'saved_best':json.loads(result_path.read_text())['rank_lower_bound'],
                      'generation_seconds':round(state['seconds'],2),'counts':state['counts']}),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples',type=int,default=2048); p.add_argument('--heights',type=int,nargs='+',default=[24,48,96])
    p.add_argument('--roots',type=int,nargs='+',default=[8,10]); p.add_argument('--scan',type=int,default=256)
    p.add_argument('--seed',type=int,default=20260917); p.add_argument('--keep',type=int,default=24)
    p.add_argument('--adaptive',action='store_true',help='Mutate successful generated root sets on half the draws')
    p.add_argument('--workers',type=int,default=20); p.add_argument('--point-workers',type=int,default=6)
    p.add_argument('--parallel-curves',type=int,default=4); p.add_argument('--search-count',type=int,default=8)
    p.add_argument('--seconds',type=float,default=60); p.add_argument('--anchors',type=int,default=1024)
    p.add_argument('--output',default='artifacts/independent-synthesis/first')
    p.add_argument('--result',default='results/independent-synthesis-best.json')
    args = p.parse_args()
    if not (1<=args.samples<=1000000 and set(args.roots)<= {8,10} and 1<=len(args.roots)<=2
        and all(max(args.roots)<=h<=10000 for h in args.heights) and 1<=len(args.heights)<=10
        and 0<=args.scan<=100000 and 0<=args.seed<2**31 and 1<=args.keep<=1000
        and 1<=args.workers<=32 and 1<=args.point_workers<=8 and 1<=args.parallel_curves<=4
        and 0<=args.search_count<=args.keep and 0<=args.seconds<=3600 and 1<=args.anchors<=4096):
        p.error('Invalid bounded synthesis settings')
    run(args)

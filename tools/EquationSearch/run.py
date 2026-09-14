"""Find and certify points starting only from an integral Weierstrass equation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools/RankHunt'))

from bootstrap import discover, equation, save


def data_boundary(source, output):
    """Audit Python data reads; GP receives arithmetic generated from the equation.

    This is a reproducibility guard, not an operating-system sandbox.
    """
    reads=set(); source=source.resolve(); output=output.resolve()
    def guard(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)): return
        path=Path(args[0]).resolve()
        if path.is_relative_to(ROOT):
            code=path.suffix in {'.py','.pyc','.dll','.deps.json','.runtimeconfig.json'}
            if path!=source and not path.is_relative_to(output) and not code:
                raise PermissionError('Equation-only data boundary: '+str(path))
            mode=args[1]
            if mode is None or isinstance(mode,str) and 'r' in mode: reads.add(str(path))
    sys.addaudithook(guard)
    return reads


def run(args):
    started=time.perf_counter()
    source=Path(args.input).resolve(); output=Path(args.output).resolve()
    if output==ROOT or not output.is_relative_to(ROOT): raise ValueError('Use a dedicated workspace output')
    if output.exists() and any(output.iterdir()): raise ValueError('Use an empty output for a fresh speed measurement')
    output.mkdir(parents=True,exist_ok=True)
    reads=data_boundary(source,output)
    raw=source.read_bytes(); data=equation(json.loads(raw.decode('utf-8-sig')))
    save(output/'equation.json',data)
    modules=[HERE/name for name in ('run.py','bootstrap.py','seeded.py','certificate.py','geometry.py','models.py','lattice.py')]
    modules += [ROOT/'tools/RankHunt'/name for name in ('blind_search.py','point_search.py',
                 'bounded_anchor_pool.py','anchor_diversity.py','point_arithmetic.py')]
    config={'input_sha256':hashlib.sha256(raw).hexdigest(),
            'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in modules},
            'bootstrap_seconds':args.bootstrap_seconds,'search_seconds':args.search_seconds,
            'workers':args.workers,'anchors':args.anchors,'target':args.target,'job_seconds':args.job_seconds,
            'seed_limit':args.seed_limit,
            'anchor_mode':args.anchor_mode,'lattice_seconds':args.lattice_seconds,
            'equation_only':True,'initial_points':0,'full_rank_calls':False}
    save(output/'config.json',config)
    print(json.dumps({'phase':'start','input_points':0,'target':args.target,
                      'bootstrap_budget':args.bootstrap_seconds,'search_budget':args.search_seconds}),flush=True)
    basis,bootstrap=discover(data,output,args.bootstrap_seconds,args.workers,args.job_seconds,args.seed_limit,args.lattice_seconds)
    search_seconds=0.0; verification_seconds=0.0
    if basis:
        from seeded import search, independent_result
        seed=output/'seed.json'
        t=time.perf_counter()
        if basis['rank_lower_bound']<args.target and args.search_seconds:
            found=search(seed,output/'expansion',args.search_seconds,args.workers,args.anchors,args.target,
                         batch_size=4,anchor_mode=args.anchor_mode)
        else: found=basis
        search_seconds=time.perf_counter()-t
        t=time.perf_counter(); result=independent_result(found,found['rank_lower_bound'])
        verification_seconds=time.perf_counter()-t
        status='target_reached' if result['rank_lower_bound']>=args.target else 'search_incomplete'
    else:
        result={**data,'points':[],'rank_lower_bound':0}
        status='no_certified_seed_found' if bootstrap['points'] else 'no_initial_point_found'
    result.update(equation_only=True,initial_points=0,target=args.target,
                  target_reached=result['rank_lower_bound']>=args.target,status=status,
                  bootstrap_seconds=bootstrap['seconds'],search_seconds=search_seconds,
                  final_verification_seconds=verification_seconds,total_seconds=time.perf_counter()-started,
                  limitation='Bounded heuristic; failure to find a point gives no upper bound on rank.',
                  bootstrap_status=bootstrap['status'],bootstrap_errors=bootstrap['errors'])
    save(output/'result.json',result)
    save(output/'data-reads.json',sorted(reads))
    print(json.dumps({k:v for k,v in result.items() if k not in {'ainvs','points','certificate','verification'}}),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default=str(HERE/'curve.json'))
    parser.add_argument('--output',required=True)
    parser.add_argument('--bootstrap-seconds',type=float,default=120)
    parser.add_argument('--search-seconds',type=float,default=600)
    parser.add_argument('--workers',type=int,default=12)
    parser.add_argument('--anchors',type=int,default=2048)
    parser.add_argument('--target',type=int,default=31)
    parser.add_argument('--job-seconds',type=float,default=3)
    parser.add_argument('--seed-limit',type=int,default=None,help='Pass at most this many independently found points to expansion')
    parser.add_argument('--anchor-mode',choices=('adaptive','fixed','frozen','parity'),default='adaptive')
    parser.add_argument('--lattice-seconds',type=float,default=0,help='Optional tangent-lattice stage within the bootstrap budget')
    args=parser.parse_args()
    if not (0<args.bootstrap_seconds<=3600 and 0<=args.search_seconds<=7200 and
            1<=args.workers<=24 and 1<=args.anchors<=4096 and 1<=args.target<=100 and
            .05<=args.job_seconds<=60 and 0<=args.lattice_seconds<=args.bootstrap_seconds and
            (args.seed_limit is None or 1<=args.seed_limit<=100)):
        parser.error('Invalid bounded search settings')
    run(args)

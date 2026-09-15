"""Bounded equation-only first point search using partial cubic 2-descent.

Collect projected linear elements, cancel good-prime ideal parities, certify
nontrivial square classes, solve a conic, and search its binary quartic.
No complete class/unit group, Selmer group, rank call, or reference points.
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
from point_arithmetic import on_curve
from point_search import quartic_reduction_code
import certificate


def primes(bound):
    sieve = bytearray(b'\1') * (bound + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, math.isqrt(bound) + 1):
        if sieve[p]:
            sieve[p*p:bound+1:p] = b'\0' * ((bound-p*p)//p+1)
    return [p for p in range(2, bound+1) if sieve[p]]


def cubic_model(ainvs):
    a1, a2, a3, a4, a6 = map(int, ainvs)
    if a1 == a3 == 0:
        return (a2, a4, a6), 1
    return (a1*a1+4*a2, 8*(a1*a3+2*a4), 16*(a3*a3+4*a6)), 4


def discriminant(f):
    a, b, c = f
    return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c


def mul(f, v, w):
    a, b, c = f
    z = [sum(v[i]*w[k-i] for i in range(3) if 0 <= k-i < 3) for k in range(5)]
    for k in (4, 3):
        z[k-1] -= a*z[k]
        z[k-2] -= b*z[k]
        z[k-3] -= c*z[k]
    return z[:3]


def determinant(m):
    a, b, c = m
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def norm(f, v):
    cols = [mul(f, v, [int(i == j) for i in range(3)]) for j in range(3)]
    return determinant(list(map(list, zip(*cols))))


def matrices(f, delta):
    out = [[[0]*3 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            e = [0, 0, 0]; e[i] = 1
            g = [0, 0, 0]; g[j] = 1
            value = mul(f, delta, mul(f, e, g))
            for k in range(3): out[k][i][j] = value[k]
    if determinant(out[2]) != -norm(f, delta):
        raise ArithmeticError('Conic determinant identity')
    return out


def gp_matrix(m):
    return '[' + ';'.join(','.join(str(x) for x in row) for row in m) + ']'


def run_gp(script, seconds):
    checked_script(script)
    started = time.perf_counter()
    try:
        p = subprocess.run([str(GP), '-fq', '-s', '64M'], input=script, text=True,
                           capture_output=True, timeout=max(.01, seconds))
        out, err, status = p.stdout, p.stderr, 'finished'
        if p.returncode or '***' in err: status = 'arithmetic_error'
    except subprocess.TimeoutExpired as e:
        out, err, status = e.stdout or b'', e.stderr or b'', 'timeout'
        if isinstance(out, bytes): out = out.decode(errors='replace')
        if isinstance(err, bytes): err = err.decode(errors='replace')
    return out, {'seconds': time.perf_counter()-started, 'status': status,
                 'error': err[-1200:] if err else None}


class Collector:
    def __init__(self, f, bound):
        self.f = f; self.disc = discriminant(f)
        self.primes = primes(max(bound, 1000))
        self.base = [p for p in self.primes if p <= bound and (2*self.disc) % p]
        self.root_cache = {}; self.labels = {}; self.pivots = {}; self.rows = []
        self.tested = 0; self.norm_bits = []; self.square_cofactors = 0

    def roots(self, p):
        if p not in self.root_cache:
            a, b, c = (x % p for x in self.f)
            self.root_cache[p] = [r for r in range(p) if ((r+a)*r*r+b*r+c) % p == 0]
        return self.root_cache[p]

    def accept(self, a, b, divisors=None):
        if math.gcd(a, b) != 1: return None
        A, B, C = self.f
        n = a**3+A*a*a*b+B*a*b*b+C*b**3
        if not n: return None
        self.tested += 1; self.norm_bits.append(abs(n).bit_length())
        u = abs(n)
        while True:
            g = math.gcd(u, 2*abs(self.disc))
            if g == 1: break
            u //= g
        odd = []
        for p in self.base if divisors is None else divisors:
            e = 0
            while u % p == 0: u //= p; e ^= 1
            if e: odd.append(p)
            if u == 1: break
        root = math.isqrt(u)
        if root*root != u: return None
        self.square_cofactors += int(u != 1)
        labels = []
        for p in odd:
            roots = self.roots(p)
            r = a*pow(b, -1, p) % p
            if r not in roots or len(roots) not in (1, 3):
                raise ArithmeticError('Linear element splitting identity')
            labels += [(p, -1)] if len(roots) == 1 else [(p, t) for t in roots if t != r]
        column = 0
        for label in labels:
            if label not in self.labels: self.labels[label] = len(self.labels)
            column ^= 1 << self.labels[label]
        row = {'a': a, 'b': b, 'norm': n, 'odd_good_components': labels,
               'remaining_square': root}
        index = len(self.rows); self.rows.append(row)
        support = 1 << index
        while column:
            bit = column.bit_length()-1
            if bit not in self.pivots:
                self.pivots[bit] = column, support
                return None
            other, other_support = self.pivots[bit]
            column ^= other; support ^= other_support
        return support

    def candidate(self, support):
        rows = [r for i, r in enumerate(self.rows) if support >> i & 1]
        if len(rows) > 8: return None
        gamma = [1, 0, 0]; product = 1
        for row in rows:
            gamma = mul(self.f, gamma, [row['a'], -row['b'], 0])
            product *= row['norm']
        # n = d*q^2, without factoring the remaining large cofactor.
        rest = abs(product); d = -1 if product < 0 else 1; q = 1
        for p in self.primes:
            e = 0
            while rest % p == 0: rest //= p; e += 1
            q *= p**(e//2)
            if e & 1: d *= p
            if rest == 1: break
        sq = math.isqrt(rest)
        if sq*sq == rest: q *= sq
        else: d *= rest
        if d*q*q != product: raise ArithmeticError('Rational square extraction')
        delta = [d*x for x in gamma]
        content = math.gcd(*delta); scale = 1
        for p in self.primes:
            while content % (p*p) == 0:
                content //= p*p; scale *= p
        sq = math.isqrt(content)
        if sq*sq == content: scale *= sq
        delta = [x//(scale*scale) for x in delta]
        nu = Q(d*d*q, scale**3)
        if norm(self.f, delta) != nu*nu: raise ArithmeticError('Class norm identity')
        if max(abs(x).bit_length() for x in delta) > 2048: return None
        witness = None
        for p in self.primes:
            if p == 2 or (self.disc*nu.numerator*nu.denominator) % p == 0: continue
            for r in self.roots(p):
                value = sum(delta[i]*pow(r, i, p) for i in range(3)) % p
                if pow(value, (p-1)//2, p) == p-1:
                    witness = {'prime': p, 'root': r, 'residue': value}; break
            if witness: break
        if witness is None: return None  # Inconclusive, not a proof of squarehood.
        return {'delta': list(map(str, delta)), 'norm_sqrt': str(nu), 'rows': rows,
                'product_to_delta_square': str(q*scale), 'nontrivial_witness': witness}


def cover_script(ainvs, f, scale, candidate, height):
    delta = list(map(int, candidate['delta'])); q0, q1, q2 = matrices(f, delta)
    A, B, C = f
    script = ('default(parisizemax,134217728);default(realprecision,100);setrand(20260915);\n'
              'E=ellinit('+vec(ainvs)+');\n'
              'Q0='+gp_matrix(q0)+';Q1='+gp_matrix(q1)+';Q2='+gp_matrix(q2)+';\n'
              'nu='+candidate['norm_sqrt']+';\n'
              f'nm(v)={{my(u=v[1],b=v[2],w=v[3]);matdet([u,-({C})*w,-({C})*b+({A*C})*w;'
              f'b,u-({B})*w,-({B})*b+({A*B-C})*w;w,b-({A})*w,u-({A})*b+({A*A-B})*w])}};\n'
              'print("STAGE conic ",getwalltime());sol=qfsolve(Q2);'
              'if(type(sol)!="t_COL",print("NO_CONIC ",sol);quit);'
              'pm=qfparam(Q2,sol);bp=pm*[x^2,x,1]~;'
              'if(bp~*Q2*bp!=0,error("Conic parametrization"));'
              'rawF=-bp~*Q1*bp;G=bp~*Q0*bp;'
              f'if(G^3+({A})*G^2*rawF+({B})*G*rawF^2+({C})*rawF^3-nu^2*nm(bp)^2!=0,error("Cover identity"));'
              'den=denominator(content(rawF));F=den^2*rawF;'
              'print("STAGE reduction ",getwalltime());'+quartic_reduction_code(True)+
              'md=m[2][2,1]*x+m[2][2,2];mt=(m[2][1,1]*x+m[2][1,2])/md;'
              'if(md^4*subst(F,x,mt)-m[3]^2!=m[1]^2*C[1],error("Reduction identity"));'
              'if(2*m[1]*m[3]!=m[1]^2*C[2],error("Reduction linear identity"));'
              'print("MODEL ",[vector(5,i,Str(polcoef(C[1],i-1))),vector(3,i,Str(polcoef(C[2],i-1)))]);\n')
    inverse = ('P=[XX,YY];' if scale == 1 else
               f'P=[XX/4,(YY-({ainvs[0]})*XX-4*({ainvs[2]}))/8];')
    script += ('emit(s,t,z)={my(v,XX,YY,P);if(z==0,return(0));v=pm*[s^2,s*t,t^2]~;'
               'if(v~*Q2*v!=0 || v~*Q1*v!=-z^2,error("Point on cover"));'
               'XX=(v~*Q0*v)/z^2;YY=nu*nm(v)/z^3;'+inverse+
               'if(!ellisoncurve(E,P),error("Point recovery"));'
               'print("POINT ",vector(2,i,Str(P[i])));'
               'print("LIFT ",[vector(3,i,Str(v[i])),Str(z)]);return(1);};\n'
               'print("STAGE search ",getwalltime());'
               'if(issquare(polcoef(C[2],2)^2+4*polcoef(C[1],4),&zz),'
               'yy=(-polcoef(C[2],2)+zz)/2;'
               'if(emit(m[2][1,1],m[2][2,1],(m[1]*yy+polcoef(m[3],2))/den),print("DONE");quit));'
               f'H=hyperellratpoints(C,{int(height)},1);'
               'for(j=1,#H,h=H[j][1];yy=H[j][2];'
               'if(emit(m[2][1,1]*h+m[2][1,2],m[2][2,1]*h+m[2][2,2],'
               '(m[1]*yy+subst(m[3],x,h))/den),break));'
               'print("STAGE end ",getwalltime());print("DONE");quit;\n')
    return script


def verify_seed(ainvs, f, scale, candidate, lift, point):
    """Independent exact replay of relation, class character and covering map."""
    col = Collector(f, 2)  # Only field arithmetic; no relation recollection.
    product = [1, 0, 0]; labels = set()
    for row in candidate['rows']:
        a, b, n = row['a'], row['b'], row['norm']
        if math.gcd(a, b) != 1 or norm(f, [a, -b, 0]) != n:
            raise ArithmeticError('Relation input')
        product = mul(f, product, [n*a, -n*b, 0])
        labels.symmetric_difference_update(map(tuple, row['odd_good_components']))
    if labels: raise ArithmeticError('Uncancelled good valuations')
    delta = list(map(Q, candidate['delta'])); ratio = Q(candidate['product_to_delta_square'])
    if product != [ratio*ratio*x for x in delta]: raise ArithmeticError('Class representative')
    witness = candidate['nontrivial_witness']; p, r = witness['prime'], witness['root']
    if not certificate.is_prime(p) or p == 2 or col.disc % p == 0:
        raise ArithmeticError('Character prime')
    A, B, C = f
    if (r**3+A*r*r+B*r+C) % p: raise ArithmeticError('Character root')
    value = sum(int(delta[i])*pow(r, i, p) for i in range(3)) % p
    if value != witness['residue'] or pow(value, (p-1)//2, p) != p-1:
        raise ArithmeticError('Nontrivial character')
    beta, z = list(map(Q, lift[0])), Q(lift[1])
    x, y = map(Q, point)
    X, Y = (x, y) if scale == 1 else (4*x, 8*y+4*Q(ainvs[0])*x+4*Q(ainvs[2]))
    if mul(f, delta, mul(f, beta, beta)) != [X*z*z, -z*z, 0]:
        raise ArithmeticError('Kummer lift')
    if Y*z**3 != Q(candidate['norm_sqrt'])*norm(f, beta): raise ArithmeticError('Norm lift')
    if not on_curve(list(map(Q, ainvs)), (x, y)): raise ArithmeticError('Original curve')
    data = {'ainvs': list(map(str, ainvs)), 'points': [list(map(str, (x, y)))]}
    proof = certificate.build_certificate(data, max_prime=2000)
    checked = certificate.verify_certificate(data, proof)
    if checked['rank_lower_bound'] != 1: raise ArithmeticError('Independent point certificate')
    return {**data, 'rank_lower_bound': 1, 'certificate': proof, 'verification': checked}


def pairs(centers, collector, deadline):
    """Sieve relation strips by a = b*r (mod p); never scan curve x for squares."""
    seen = set()
    previous = -1
    root_primes = [(p, collector.roots(p)) for p in collector.base if collector.roots(p)]
    for width in (8, 32, 128, 512, 2048, 8192, 32768):
        for b in range(1, 33):
            for c in centers:
                center = round(c*b)
                intervals = [(-width, width)] if previous < 0 else [(-width, -previous-1), (previous+1, width)]
                for lo, hi in intervals:
                    if time.perf_counter() >= deadline: return
                    start = center+lo; count = hi-lo+1
                    divisors = [[] for _ in range(count)]
                    for p, roots in root_primes:
                        if b % p == 0: continue
                        for r in roots:
                            for j in range((b*r-start) % p, count, p): divisors[j].append(p)
                    for j in range(count):
                        a = start+j
                        if (a, b) not in seen:
                            seen.add((a, b)); yield a, b, divisors[j]
        previous = width


def search(data, seconds=3, factor_bound=20000, height=10000, max_covers=12):
    started = time.perf_counter(); deadline = started+seconds
    data = equation(data); ainvs = data['ainvs']; f, scale = cubic_model(ainvs)
    A, B, C = f
    result = {'ainvs': ainvs, 'input_points': 0, 'status': 'no_point_within_budget',
              'config': {'seconds': seconds, 'factor_bound': factor_bound, 'height': height,
                         'max_covers': max_covers}, 'stages': {}, 'cover_trials': [],
              'points': [], 'rank_lower_bound': 0}
    prep = ('default(realprecision,160);f=x^3+('+str(A)+')*x^2+('+str(B)+')*x+('+str(C)+');'
            'print("IRREDUCIBLE ",polisirreducible(f));'
            'rr=polrootsreal(f);print("CENTERS ",vector(#rr,i,Str(round(rr[i]*10^40))));'
            f'D=2*poldisc(f);forprime(pp=3,{factor_bound},if(D%pp,rr=polrootsmod(f,pp);'
            'if(#rr,print("ROOTS ",[pp,Vec(lift(rr))]))));print("ROOTS_END");quit;')
    out, stats = run_gp(prep, min(.7, max(.01, deadline-time.perf_counter())))
    result['stages']['preparation_seconds'] = stats['seconds']
    if 'IRREDUCIBLE 1' not in out:
        result.update(status='unsupported_reducible_cubic' if 'IRREDUCIBLE 0' in out else 'preparation_incomplete',
                      total_seconds=time.perf_counter()-started, preparation=stats)
        return result
    centers = [Q(0)]
    for line in out.splitlines():
        if line.startswith('CENTERS '): centers = [Q(int(x), 10**40) for x in json.loads(line[8:])] + centers
    collector = Collector(f, factor_bound)
    if 'ROOTS_END' not in out:
        result.update(status='preparation_incomplete', total_seconds=time.perf_counter()-started, preparation=stats)
        return result
    collector.root_cache.update({p: [] for p in collector.base})
    for line in out.splitlines():
        if line.startswith('ROOTS '):
            p, roots = json.loads(line[6:]); collector.root_cache[p] = roots
    candidates = []; seen = set(); deps = 0; supports = []
    collect_started = time.perf_counter()
    collect_deadline = deadline-.1
    def retain(support):
        candidate = collector.candidate(support)
        if candidate is None: return
        key = tuple(candidate['delta'])
        if key not in seen:
            seen.add(key); candidates.append(candidate)
    for a, b, divisors in pairs(centers, collector, collect_deadline):
        if time.perf_counter() >= collect_deadline: break
        if candidates and time.perf_counter()-collect_started >= seconds*.6: break
        support = collector.accept(a, b, divisors)
        if support is None: continue
        deps += 1
        retain(support)
        # Retain products even when individual generators fail local tests later.
        # This is a bounded partial span, not a claim to enumerate all of D.
        for prior in supports[:4]:
            if len(candidates) >= max_covers: break
            retain(support ^ prior)
        supports.append(support)
        if len(candidates) >= max_covers: break
    result['stages']['relations_seconds'] = time.perf_counter()-collect_started
    result['relations'] = {'tested_pairs': collector.tested, 'accepted_elements': len(collector.rows),
                           'ideal_matrix_rows': len(collector.labels), 'ideal_matrix_rank': len(collector.pivots),
                           'dependencies': deps, 'certified_nontrivial_candidates': len(candidates),
                           'norm_bits': [min(collector.norm_bits), max(collector.norm_bits)] if collector.norm_bits else [],
                           'square_cofactors': collector.square_cofactors}
    candidates.sort(key=lambda c: (max(abs(int(x)).bit_length() for x in c['delta']), len(c['rows'])))
    for candidate in candidates:
        remaining = deadline-time.perf_counter()
        if remaining <= .03: break
        t = time.perf_counter()
        script = cover_script(ainvs, f, scale, candidate, height)
        out, stats = run_gp(script, min(.6, max(.01, deadline-time.perf_counter())))
        trial = {'candidate': candidate, 'status': stats['status'], 'seconds': time.perf_counter()-t,
                 'stages_ms': {}, 'coefficient_bits': None}
        stage_events = []; point = lift = None
        for line in out.splitlines():
            if line.startswith('STAGE '):
                _, name, stamp = line.split(); stage_events.append((name, int(stamp)))
            elif line.startswith('NO_CONIC '): trial['status'] = 'conic_obstruction'
            elif line.startswith('MODEL '):
                model = json.loads(line[6:]); trial['coefficient_bits'] = max(abs(int(x)).bit_length() for p in model for x in p)
            elif line.startswith('POINT '): point = json.loads(line[6:])
            elif line.startswith('LIFT '): lift = json.loads(line[5:])
        for (name, stamp), (_, end) in zip(stage_events, stage_events[1:]): trial['stages_ms'][name] = end-stamp
        trial['last_stage'] = stage_events[-1][0] if stage_events else None
        if stats['error']: trial['error'] = stats['error']
        if point and lift:
            t = time.perf_counter()
            checked = verify_seed(ainvs, f, scale, candidate, lift, point)
            result['stages']['verification_seconds'] = time.perf_counter()-t
            trial.update(status='point_found', lift=lift)
            result.update(checked, status='point_found')
        result['cover_trials'].append(trial)
        if result['points']: break
    result['stages']['covers_seconds'] = sum(c['seconds'] for c in result['cover_trials'])
    if not result['points']:
        if not candidates: result['status'] = 'no_certified_class'
        elif all(c['status'] == 'conic_obstruction' for c in result['cover_trials']) and result['cover_trials']:
            result['status'] = 'all_tested_conics_obstructed'
        elif any(c['status'] == 'arithmetic_error' for c in result['cover_trials']): result['status'] = 'arithmetic_error'
    result['total_seconds'] = time.perf_counter()-started
    result['limitation'] = 'Bounded partial descent; no point found gives no rank upper bound.'
    return result


CASES = [
    ('small', 1, [0, 0, 0, -1, 1]),
    ('icarm-47', 8, [1, -1, 0, -106384, 13075804]),
    ('icarm-157', 12, [1, -1, 0, -227292004, 882331831684]),
    ('icarm-644', 18, [1, 0, 1, -174327769786445569143685, 26865014546188354102268547317071580]),
    ('icarm-718', 19, [1, 0, 0, -416891163445623604369297428825, 103908982619412549877181724313142587147100601]),
    ('icarm-631', 20, [1, 0, 0, -428377857530303119736330313, 2949746475279594219167674546858042856217]),
    ('record-31', 31, [1, 1, 1, -1284727764113567728281797636015784768866707681415849262157224232063,
                     560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181]),
]


def benchmark(args):
    output = Path(args.output).resolve(); output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    report = {'method': 'partial cubic descent, first point only', 'direct_search_comparison': False,
              'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'selection': 'small check plus previously used ICARM leaders and the user equation; not a random sample',
              'trials': []}
    cases = CASES if not args.cases else [r for r in CASES if r[0] in args.cases.split(',')]
    for repeat in range(args.repeats):
        for name, target, ainvs in cases:
            with tempfile.TemporaryDirectory(prefix='first-cover-', dir=output.parent) as temporary:
                directory = Path(temporary); source = directory/'input.json'; dest = directory/'result.json'
                save(source, {'ainvs': list(map(str, ainvs))})
                t = time.perf_counter()
                command = [sys.executable, str(Path(__file__)), '--input', str(source), '--output', str(dest),
                           '--seconds', str(args.seconds), '--factor-bound', str(args.factor_bound),
                           '--height', str(args.height), '--max-covers', str(args.max_covers)]
                p = subprocess.run(command, capture_output=True, text=True, timeout=args.seconds+10)
                elapsed = time.perf_counter()-t
                trial = {'case': name, 'published_lower_bound': target, 'repeat': repeat+1, 'process_seconds': elapsed}
                if p.returncode or not dest.exists(): trial.update(status='worker_error', error=p.stderr[-2000:])
                else:
                    value = json.loads(dest.read_text()); trial.update(value)
                    if value.get('points'):
                        certificate.verify_certificate(value, value['certificate'])
                        f, scale = cubic_model(ainvs)
                        won = next(c for c in value['cover_trials'] if c['status'] == 'point_found')
                        verify_seed(ainvs, f, scale, won['candidate'], won['lift'], value['points'][0])
                report['trials'].append(trial)
                save(output, report)
                print(json.dumps({k: trial.get(k) for k in ('case', 'repeat', 'status', 'process_seconds')}), flush=True)
    report['wall_seconds'] = time.perf_counter()-started
    report['complete'] = True
    report['successes'] = sum(t.get('status') == 'point_found' for t in report['trials'])
    save(output, report)


def self_test():
    f = (0, -1, 1)
    col = Collector(f, 31)
    rows = [(-3, 1), (3, 2)]
    for a, b in rows: col.accept(a, b)
    candidate = col.candidate(3)
    assert candidate and list(map(int, candidate['delta'])) == [9, -3, -2]
    assert norm(f, [9, -3, -2]) == 529
    got = search({'ainvs': ['0', '0', '0', '-1', '1']}, seconds=3)
    assert got['status'] == 'point_found', got
    forged = json.loads(json.dumps(got))
    forged['points'][0][0] = str(Q(forged['points'][0][0])+1)
    try:
        won = next(c for c in forged['cover_trials'] if c['status'] == 'point_found')
        verify_seed(got['ainvs'], f, 1, won['candidate'], won['lift'], forged['points'][0])
    except ArithmeticError: pass
    else: raise AssertionError('Modified point accepted')
    print('PASS: relation/class identities, complete covering recovery, independent certificate, tampered point rejected')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input'); parser.add_argument('--output')
    parser.add_argument('--seconds', type=float, default=3)
    parser.add_argument('--factor-bound', type=int, default=20000)
    parser.add_argument('--height', type=int, default=10000)
    parser.add_argument('--max-covers', type=int, default=12)
    parser.add_argument('--benchmark', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--cases')
    args = parser.parse_args()
    if not (0 < args.seconds <= 60 and 3 <= args.factor_bound <= 100000 and 1 <= args.max_covers <= 64
            and 1 <= args.height <= 1000000 and 1 <= args.repeats <= 10): parser.error('Invalid bounded settings')
    if args.self_test: self_test()
    elif args.benchmark:
        if not args.output: parser.error('--output required')
        benchmark(args)
    else:
        if not args.input or not args.output: parser.error('--input and --output required')
        from run import data_boundary
        source, output = Path(args.input).resolve(), Path(args.output).resolve()
        reads = data_boundary(source, output.parent)
        result = search(json.loads(source.read_text(encoding='utf-8-sig')), args.seconds,
                        args.factor_bound, args.height, args.max_covers)
        result['unexpected_data_reads'] = [p for p in reads if Path(p) != source and Path(p).suffix not in ('.py', '.pyc')]
        save(output, result)
        print(json.dumps({'status': result['status'], 'seconds': result['total_seconds']}))

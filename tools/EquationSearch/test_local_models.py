"""Independent rational identities for local quartic coordinate changes."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from models import prepare_models
from local_models import expand_models, _script
from bootstrap import gp


def add(a,b):
    return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]


def mul(a,b):
    c=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def power(a,k):
    p=[Q(1)]
    for _ in range(k):p=mul(p,a)
    return p


def homogeneous(f,m,k):
    a,b,c,d=map(Q,m);out=[Q(0)]*(k+1)
    for j,v in enumerate(f):out=add(out,[v*t for t in mul(power([b,a],j),power([d,c],k-j))])
    return out


def value(f,x):return sum(v*x**j for j,v in enumerate(f))


def invariants(model):
    f,q=map(lambda z:list(map(Q,z)),[model['f'],model['q']])
    e,d,c,b,a=add([4*v for v in f],mul(q,q))
    return 12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3


class LocalModelTests(unittest.TestCase):
    def test_exact_polynomial_and_independently_supplied_point_maps(self):
        examples=[([0,0,0,-25,4],[0,2],[5,2]),
                  ([0,0,0,-400,256],[0,16],[20,16]),
                  ([2,-1,2,-38,-1],[12,23],[Q(25,4),Q(-93,8)])]
        checked=0;used_two=False
        for a,p,other in examples:
            pool={'ainvs':list(map(str,a)),'points':[list(map(str,p))],'approximate_heights':[1]}
            base=prepare_models(pool,1,{})[0]
            charts=expand_models([base],1,{},count=6,roots=1,beam=4,depth=6,prime_bound=19)
            self.assertTrue(charts)
            ma=list(map(Q,base['minimal_ainvs']));x0,y0=map(Q,base['minimal_anchor'])
            b2=ma[0]**2+4*ma[1];b4=ma[0]*ma[2]+2*ma[3];v0=2*y0+ma[0]*x0+ma[2]
            D=[b2*b2-8*b2*x0-48*x0*x0-32*b4,32*v0,-2*(12*x0+b2),Q(0),Q(1)]
            den=Q(base['den']);F=[den*den*v for v in D]
            u,r,s,t=map(Q,base['change']);ox,oy=map(Q,other)
            xm=(ox-r)/(u*u);ym=(oy-s*(ox-r)-t)/u**3
            slope=(2*ym+ma[0]*xm+ma[2]-v0)/(xm-x0)
            zz=8*xm-slope*slope+b2+4*x0
            self.assertEqual(zz*zz,value(D,slope))
            for model in charts:
                used_two |= 2 in model['neighbour_primes']
                self.assertEqual(invariants(base),invariants(model))
                e,m,h=model['transform'];e=Q(e);h=list(map(Q,h));f=list(map(Q,model['f']));q=list(map(Q,model['q']))
                self.assertEqual(add(homogeneous(F,m,4),[-v for v in mul(h,h)]),[e*e*v for v in f])
                self.assertEqual([2*e*v for v in h],[e*e*v for v in q])
                aa,b,c,d=map(Q,m)
                if c*slope==aa:continue
                xx=(b-d*slope)/(c*slope-aa);yy=(den*zz*(c*xx+d)**2-value(h,xx))/e
                self.assertEqual(yy*yy+value(q,xx)*yy,value(f,xx))
                checked+=1
        self.assertTrue(used_two)
        self.assertGreaterEqual(checked,4)

    def test_cache_preserves_current_basis_provenance(self):
        pool={'ainvs':['0','0','0','-25','4'],'points':[['0','2']],'approximate_heights':[1]}
        base=prepare_models(pool,1,{})[0];cache={}
        first=expand_models([base],1,cache,count=4,roots=1,beam=4,depth=4,prime_bound=19)
        self.assertTrue(first)
        current=dict(base,pool_index=17,parity='101',anchor_height=Q(7,3))
        with patch('local_models.gp',side_effect=AssertionError('Cache hit must not run GP')):
            second=expand_models([current],.1,cache,count=4,roots=1,beam=4,depth=4,prime_bound=19)
        self.assertEqual([m['key'] for m in first],[m['key'] for m in second])
        self.assertTrue(all(m['pool_index']==17 and m['parity']=='101' and m['anchor_height']==Q(7,3) for m in second))

    def test_partial_cache_retries_with_a_larger_budget(self):
        pool={'ainvs':['0','0','0','-25','4'],'points':[['0','2']],'approximate_heights':[1]}
        base=prepare_models(pool,1,{})[0];cache={}
        full,stats=gp(_script(base,4,4,4,19),1)
        partial=next(line for line in full.splitlines() if line.startswith('NEIGHBOURS '))+'\n'
        with patch('local_models.gp',return_value=(partial,{'timed_out':True,'seconds':.1})):
            first=expand_models([base],.2,cache,count=4,roots=1,beam=4,depth=4,prime_bound=19)
        entry=next(iter(cache['local_models'].values()))
        self.assertFalse(entry['complete']);self.assertEqual(entry['completed_levels'],1)
        self.assertTrue(first)
        with patch('local_models.gp',side_effect=AssertionError('Same budget reuses partial work')):
            same=expand_models([base],.2,cache,count=4,roots=1,beam=4,depth=4,prime_bound=19)
        self.assertEqual(first,same)
        with patch('local_models.gp',return_value=(full,stats)) as process:
            expanded=expand_models([base],.5,cache,count=4,roots=1,beam=4,depth=4,prime_bound=19)
        process.assert_called_once();entry=next(iter(cache['local_models'].values()))
        self.assertTrue(entry['complete']);self.assertGreater(entry['completed_levels'],1)
        self.assertTrue(expanded)

    def test_intermediate_charts_find_an_independent_direction_from_one_seed(self):
        # Frozen ICARM #724 equation and first point only. Neither a target rank
        # nor another published point enters the construction or point search.
        from seeded import batch_search, clean
        from torsion_certificate import select_basis, verify_certificate
        data={'ainvs':['1','0','0',
            '-500155818938812672789174015511502418542914025471083342999327356190',
            '135726649501359436415916160439700483690089900096812327554427167585750287637491352376425892437884100'],
            'points':[['-175168091909712984914985049975244485862218184197220/271194694797617761',
            '-1943143601681747607708811842618040351218746916988056096393304886241390544910/141228317156673157916149009']]}
        base=prepare_models(dict(data,approximate_heights=[1]),1,{})
        charts=expand_models(base,.5,{},count=6,roots=1,beam=8,depth=16)
        self.assertTrue(any(m['neighbour_snapshot_level']==1 for m in charts))
        observations=batch_search(data,charts,1048576,1,.5)['observations']
        points=clean(dict(data,points=data['points']+[o['point'] for o in observations]))
        proof=select_basis(points,max_prime=2000)
        self.assertIsNotNone(proof)
        self.assertGreaterEqual(verify_certificate(proof,proof['certificate'])['rank_lower_bound'],2)


if __name__=='__main__':unittest.main()

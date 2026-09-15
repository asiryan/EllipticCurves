"""Deterministic discoveries through the real exact unified feedback loop."""
from contextlib import ExitStack
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from bootstrap import ROOT,save
from point_arithmetic import add,multiply
from point_models import divide_basis
from torsion_certificate import select_basis,verify_certificate
import unified


def exact_certificate(path):
    result=select_basis(json.loads(Path(path).read_text()),max_prime=2000)
    if result is None:return {'LowerBound':0,'all_selected_independent':False,'points':[]}
    verify_certificate(result,result['certificate'])
    return {'LowerBound':result['rank_lower_bound'],'all_selected_independent':True,
            'points':result['points']}


def exact_fixture_pool(source,directory,count,candidates,timeout,selector=None):
    data=json.loads(Path(source).read_text());n=len(data['points'])
    return {**data,'vectors':[[int(i==j) for j in range(n)] for i in range(n)],
            'approximate_heights':[1]*n}


def fixture_models(pool,timeout,cache,minimal=True):
    return [{'key':'basis-'+str(len(pool['points']))+'-'+str(i),'anchor':p,
             'anchor_height':1,'coefficient_bits':8,'pool_index':i}
            for i,p in enumerate(pool['points'])]


class UnifiedFeedbackTests(unittest.TestCase):
    def fixture(self, discoveries, calls):
        def search(data,models,n,d,timeout,coverage=None):
            calls.append({'models':[m['key'] for m in models],'n':n,'d':d,
                          'coverage':json.loads(json.dumps(coverage or {}))})
            observations=[]
            if len(calls)==1:
                observations=[{'point':p,'method':'pointed','model':models[0]['key'],
                               'anchor':models[0]['anchor'],'n':n,'d':d} for p in discoveries]
            return {'complete':[f'{m["key"]}:{n}:{d}' for m in models],
                    'observations':observations,'slices':[(m['key'],n,1,d) for m in models],
                    'timed_out':False,'seconds':0}
        stack=ExitStack()
        for name,value in [('unified.certify',exact_certificate),
                           ('seeded.generate',exact_fixture_pool),
                           ('unified.prepare_models',fixture_models),
                           ('seeded.batch_search',search),
                           ('seeded.boxes',lambda:iter([(4,1),(4,2),(4,4)])),
                           ('unified.divide_basis',lambda data,*args:(data,{'steps':[]})),
                           ('unified.prepare_covers',lambda *args,**kw:[]),
                           ('unified.expand_models',lambda *args,**kw:[])]:
            stack.enter_context(patch(name,value))
        return stack

    def test_hidden_direction_feedback_preserves_observation_and_coverage_on_resume(self):
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2']]}
        hidden=['4/25','-8/125'];calls=[]
        with tempfile.TemporaryDirectory(prefix='unified-regression-',dir=ROOT/'artifacts/equation-search') as raw:
            root=Path(raw);source=root/'input.json';output=root/'search';save(source,data)
            with self.fixture([hidden],calls):
                result=unified.search(source,output,seconds=1,workers=1,anchors=2,target=3,batch_size=1)
                state=json.loads((output/'checkpoint.json').read_text())
                self.assertEqual(result['rank_lower_bound'],2)
                self.assertEqual(state['lower_bound'],2)
                self.assertIn(hidden,state['found']['points'])
                self.assertGreater(len(state['found']['points']),len(state['certified_basis']))
                self.assertTrue(any(r['steps'] for r in state['refinements']))
                self.assertTrue(any(o['method']=='relation_division' for e in state['events']
                                    for o in e['selected_new_points']))
                self.assertEqual({m['key'] for m in state['models']},{'basis-1-0','basis-2-0','basis-2-1'})
                self.assertTrue(any(any(m.startswith('basis-2-') for m in c['models']) and c['coverage'].get('basis-1-0')
                                    for c in calls))
                self.assertEqual(sum('basis-1-0' in c['models'] and c['n']==4 and c['d']==1 for c in calls),1)
                old_calls=len(calls);old_coverage=state['coverage'];old_generation=state['generation']
                resumed=unified.search(source,output,seconds=.2,workers=1,anchors=2,target=3,batch_size=1)
                after=json.loads((output/'checkpoint.json').read_text())
            self.assertEqual(resumed['rank_lower_bound'],2)
            self.assertEqual(len(calls),old_calls)
            self.assertEqual(after['coverage'],old_coverage)
            self.assertEqual(after['generation'],old_generation)
            self.assertEqual(after['status'],'frontier_exhausted')
            proof=select_basis(resumed,max_prime=2000)
            self.assertEqual(verify_certificate(proof,proof['certificate'])['rank_lower_bound'],2)

    def test_torsion_and_torsion_translate_do_not_inflate_claim(self):
        a=list(map(Q,[0,0,0,-36,0]));p=(Q(12),Q(36));t=(Q(0),Q(0))
        data={'ainvs':list(map(str,a)),'points':[list(map(str,p))]}
        discoveries=[list(map(str,t)),list(map(str,add(a,p,t)))];calls=[]
        with tempfile.TemporaryDirectory(prefix='unified-torsion-',dir=ROOT/'artifacts/equation-search') as raw:
            root=Path(raw);source=root/'input.json';output=root/'search';save(source,data)
            with self.fixture(discoveries,calls):
                result=unified.search(source,output,seconds=1,workers=1,anchors=2,target=2,batch_size=1)
            state=json.loads((output/'checkpoint.json').read_text())
        self.assertEqual(result['rank_lower_bound'],1)
        self.assertEqual(state['events'],[])
        self.assertGreater(len(state['found']['points']),1)
        self.assertEqual(len(state['certified_basis']),1)

    def test_torsion_pool_fallback_keeps_exact_vector_provenance(self):
        data={'ainvs':['0','0','0','-36','0'],'points':[['12','36']]};calls=[]
        with tempfile.TemporaryDirectory(prefix='unified-fallback-',dir=ROOT/'artifacts/equation-search') as raw:
            root=Path(raw);save(root/'input.json',data)
            with self.fixture([],calls),patch('seeded.generate',side_effect=subprocess.TimeoutExpired('fixture',.01)):
                result=unified.search(root/'input.json',root/'search',seconds=1,workers=1,anchors=2,target=2,batch_size=1)
            state=json.loads((root/'search'/'checkpoint.json').read_text())
        self.assertEqual(result['rank_lower_bound'],1)
        self.assertGreater(state['pool_timeouts'],0)
        self.assertTrue(state['preparations'][0]['pool_vectors'])
        self.assertNotEqual(state['status'],'error')

    def test_exact_odd_division_replaces_the_large_basis_point(self):
        a=list(map(Q,[0,0,0,-25,4]));p=(Q(0),Q(2));triple=multiply(a,p,3);calls=[]
        data={'ainvs':list(map(str,a)),'points':[list(map(str,triple))]}
        with tempfile.TemporaryDirectory(prefix='unified-division-',dir=ROOT/'artifacts/equation-search') as raw:
            root=Path(raw);save(root/'input.json',data)
            with self.fixture([],calls),patch('unified.divide_basis',divide_basis):
                result=unified.search(root/'input.json',root/'search',seconds=1,workers=1,anchors=2,target=2,batch_size=1)
            state=json.loads((root/'search'/'checkpoint.json').read_text())
        self.assertEqual(result['rank_lower_bound'],1)
        self.assertEqual(result['points'],[['0','-2']])
        self.assertGreater(len(state['found']['points']),1)
        self.assertTrue(state['seed_division']['steps'])

    def test_partial_model_preparation_retries_on_resume(self):
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2'],['5','2']]};calls=[]
        def partial(pool,*args,**kwargs):return fixture_models(pool,*args,**kwargs)[:1]
        with tempfile.TemporaryDirectory(prefix='unified-partial-',dir=ROOT/'artifacts/equation-search') as raw:
            root=Path(raw);source=root/'input.json';output=root/'search';save(source,data)
            with self.fixture([],calls):
                # A timed-out model preparation exposes only its first complete
                # model. Interrupt at the next work boundary to inspect resume.
                with patch('unified.prepare_models',partial),patch('seeded.batch_search',side_effect=KeyboardInterrupt):
                    unified.search(source,output,seconds=1,workers=1,anchors=2,target=3,batch_size=1)
                before=json.loads((output/'checkpoint.json').read_text())
                self.assertIsNone(before['basis_prepared'])
                self.assertEqual(len(before['models']),1)
                result=unified.search(source,output,seconds=1,workers=1,anchors=2,target=3,batch_size=1)
                after=json.loads((output/'checkpoint.json').read_text())
        self.assertEqual(result['rank_lower_bound'],2)
        self.assertIsNotNone(after['basis_prepared'])
        self.assertEqual(len(after['models']),2)
        self.assertTrue({m['key'] for m in before['models']}<={m['key'] for m in after['models']})
        self.assertEqual(after['status'],'frontier_exhausted')


if __name__=='__main__':unittest.main()

"""Exact regressions for directions hidden by the local mod-two filter."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from point_arithmetic import add, multiply
from point_models import refine_observed,prepare_covers
from torsion_certificate import select_basis, verify_certificate


class ObservedRefinementTests(unittest.TestCase):
    def test_even_index_observation_exposes_independent_half(self):
        data={'ainvs':['0','0','0','-25','4'],
              'points':[['0','2'],['4/25','-8/125']],'basis_size':1}
        before=select_basis(data,max_prime=2000)
        self.assertEqual(before['rank_lower_bound'],1)
        cache={};extra,stats=refine_observed(data,1,cache)
        self.assertTrue(stats['steps'])
        combined={**data,'points':data['points']+extra['points']}
        after=select_basis(combined,max_prime=2000)
        self.assertEqual(after['rank_lower_bound'],2)
        self.assertTrue(verify_certificate(after,after['certificate'])['all_points_independent_modulo_torsion'])
        for step in stats['steps']:
            a=list(map(Q,data['ainvs']));source=tuple(map(Q,step['positive']))
            for raw in step['negative']:
                from point_arithmetic import negate
                source=add(a,source,negate(a,tuple(map(Q,raw))))
            self.assertEqual(source,tuple(map(Q,step['source'])))
            self.assertEqual(multiply(a,tuple(map(Q,step['point'])),2),source)
        with patch('point_models.gp',side_effect=AssertionError('Completed relation was repeated')):
            repeated,second=refine_observed(data,1,cache)
        self.assertEqual(repeated['points'],[])
        self.assertGreater(second['cache_hits'],0)

    def test_local_dependency_alone_does_not_claim_a_half(self):
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2'],['5','2']],'basis_size':1}
        # Deliberately weak local information: the exact division check must
        # reject the candidate difference between these independent points.
        with patch('certificate.build_certificate',return_value={'independent_rows':[{'bits':'11'}]}):
            extra,stats=refine_observed(data,1,{})
        self.assertEqual(stats['tested_relations'],1)
        self.assertEqual(extra['points'],[])
        self.assertEqual(stats['steps'],[])

    def test_torsion_shift_is_witnessed_exactly(self):
        a=list(map(Q,[0,0,0,-25,0]));p=(Q(-4),Q(6));torsion=(Q(0),Q(0))
        q=add(a,multiply(a,p,2),torsion)
        data={'ainvs':list(map(str,a)),'points':[list(map(str,q))]}
        extra,stats=refine_observed(data,1,{})
        self.assertTrue(extra['points'])
        for step in stats['steps']:
            shift=tuple(map(Q,step['torsion_shift'])) if step['torsion_shift'] else None
            self.assertIsNone(multiply(a,shift,step['torsion_order_multiple']))
            self.assertEqual(multiply(a,tuple(map(Q,step['point'])),2),
                             add(a,tuple(map(Q,step['source'])),shift))

    def test_interrupted_relation_remains_retryable(self):
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2'],['4/25','-8/125']]}
        cache={}
        with patch('point_models.gp',return_value=('',{'timed_out':True,'script_sha256':'interrupted'})):
            extra,stats=refine_observed(data,1,cache)
        self.assertEqual(extra['points'],[])
        self.assertEqual(cache['observed_refinement'],{})
        extra,stats=refine_observed(data,1,cache)
        self.assertTrue(extra['points'])

    def test_partial_cover_cache_extends_without_losing_checked_maps(self):
        data={'ainvs':['0','0','0','-36','0'],'points':[['12','36']]};cache={}
        transform=['1',['1','0','0','1'],['0','0','0']]
        first=['0','1',True,['-576','0','0','0','1'],['0','0','0'],transform]
        second=['24','1',True,['1152','0','72','0','1'],['0','0','0'],transform]
        def output(rows,complete=False):
            return ''.join('COVER '+json.dumps(r)+'\n' for r in rows)+('COVERS_END\n' if complete else '')
        with patch('point_models.gp',side_effect=[
                (output([first]),{'timed_out':True}),('',{'timed_out':True}),
                (output([second],True),{'timed_out':False}),
                (output([first,second],True),{'timed_out':False})]) as native:
            a=prepare_covers(data,.05,cache,1)
            self.assertEqual(len(a),1)
            self.assertFalse(cache['isogeny_covers']['complete'])
            self.assertEqual(prepare_covers(data,.05,cache,1),a)
            self.assertEqual(native.call_count,1)
            self.assertEqual(prepare_covers(data,.1,cache,1),a)
            self.assertEqual(native.call_count,2)
            b=prepare_covers(data,.2,cache,2)
            self.assertEqual({m['key'] for m in b}|{a[0]['key']},{m['key'] for m in b})
            self.assertEqual(len(b),2)
            self.assertTrue(cache['isogeny_covers']['complete'])
            self.assertEqual(len(prepare_covers(data,1,cache,1)),1)
            self.assertEqual(native.call_count,3)
            self.assertEqual(len(prepare_covers(data,.1,cache,3)),2)
            self.assertEqual(native.call_count,4)


if __name__=='__main__':unittest.main()

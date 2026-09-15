"""Exact geometry and provenance for balanced torsion translations."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from point_arithmetic import add,multiply,on_curve
from torsion_certificate import translate_pool


class TorsionBalanceTests(unittest.TestCase):
    def setUp(self):
        self.a=list(map(Q,[0,0,0,-36,0]));self.p=(Q(12),Q(36))
        points=[multiply(self.a,self.p,n) for n in (1,2,3)]
        self.pool={'ainvs':list(map(str,self.a)),
                   'points':[list(map(str,p)) for p in points],
                   'vectors':[[1],[2],[3]],'approximate_heights':[1,4,9]}

    def test_original_series_survives_a_tight_limit(self):
        result=translate_pool(self.pool,3)
        self.assertEqual(result['points'],self.pool['points'])
        self.assertEqual(result['vectors'],self.pool['vectors'])
        self.assertEqual(result['approximate_heights'],self.pool['approximate_heights'])
        self.assertEqual(result['torsion_translations'],[None,None,None])

    def test_each_shift_visits_all_free_anchors_before_next_shift(self):
        result=translate_pool(self.pool,15)
        self.assertEqual(result['points'][:3],self.pool['points'])
        self.assertEqual(result['vectors'][3:6],[[1],[2],[3]])
        self.assertEqual(result['approximate_heights'][3:6],[1,4,9])
        first=result['torsion_translations'][3]
        self.assertIsNotNone(first)
        self.assertEqual(result['torsion_translations'][3:6],[first]*3)
        self.assertEqual(len(result['points']),15)
        self.assertEqual(len({tuple(p) for p in result['points']}),15)
        for raw,vector,shift in zip(result['points'],result['vectors'],result['torsion_translations']):
            p=tuple(map(Q,raw));t=None if shift is None else tuple(map(Q,shift))
            self.assertTrue(on_curve(self.a,p))
            self.assertIsNone(multiply(self.a,t,2))
            self.assertEqual(p,add(self.a,multiply(self.a,self.p,vector[0]),t))


if __name__=='__main__':unittest.main()

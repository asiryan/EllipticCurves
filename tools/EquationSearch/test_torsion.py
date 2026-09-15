"""Adversarial checks for torsion-aware independence, without a rank oracle."""
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'RankHunt'))
from point_arithmetic import add,multiply
import certificate as old
from torsion_certificate import select_basis,verify_certificate,torsion_points,singleton_certificate,translate_pool


class TorsionTests(unittest.TestCase):
    def setUp(self):
        self.a=['0','0','0','-36','0'];self.p=(Q(12),Q(36))
        self.torsion=((2,('0','0')),(2,('6','0')))

    def test_torsion_is_not_counted_as_free_rank(self):
        self.assertIsNone(select_basis({'ainvs':self.a,'points':[['0','0'],['6','0']]},self.torsion))
        shifted=add(list(map(Q,self.a)),self.p,(Q(0),Q(0)))
        data={'ainvs':self.a,'points':[list(map(str,self.p)),list(map(str,shifted))]}
        proof=select_basis(data,self.torsion)
        self.assertEqual(len(proof['points']),1)
        self.assertEqual(verify_certificate(proof,proof['certificate'])['rank_lower_bound'],1)
        # P and P+T really are dependent modulo torsion.
        self.assertEqual(multiply(list(map(Q,self.a)),self.p,2),multiply(list(map(Q,self.a)),shifted,2))

    def test_tampered_torsion_and_rank_are_rejected(self):
        data={'ainvs':self.a,'points':[list(map(str,self.p))]}
        proof=select_basis(data,self.torsion)
        for mutation in ('order','point','row','rank','missing'):
            c=copy.deepcopy(proof['certificate'])
            if mutation=='order':c['torsion_witnesses'][0][0]=3
            elif mutation=='point':c['torsion_witnesses'][0][1]=list(map(str,self.p))
            elif mutation=='row':c['augmented_certificate']['independent_rows'][0]['bits']='000'
            elif mutation=='rank':c['rank_lower_bound']=2
            else:c['torsion_witnesses'].pop()
            with self.assertRaises(ValueError,msg=mutation):verify_certificate(data,c)

    def test_order_four_generator_is_needed_not_only_its_double(self):
        a=['0','1','0','5365964600284239371851791','14610806763317668249280146321957009215']
        data={'ainvs':a,'points':[['-481561549686','3451824108575669073']]}
        ts=torsion_points(tuple(a));self.assertEqual([m for m,_ in ts],[4])
        proof=select_basis(data,ts)
        self.assertEqual(verify_certificate(proof,proof['certificate'])['rank_lower_bound'],1)
        doubled=multiply(list(map(Q,a)),tuple(map(Q,ts[0][1])),2)
        self.assertIsNone(select_basis(data,[(2,list(map(str,doubled)))]))

    def test_old_certificates_still_verify(self):
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2'],['5','2']]}
        c=old.build_certificate(data,2000)
        self.assertEqual(verify_certificate(data,c),old.verify_certificate(data,c))

    def test_even_multiple_remains_a_valid_single_seed(self):
        doubled=multiply(list(map(Q,self.a)),self.p,2)
        data={'ainvs':self.a,'points':[list(map(str,doubled))]}
        self.assertIsNone(select_basis(data,self.torsion))
        self.assertEqual(verify_certificate(data,singleton_certificate(data))['rank_lower_bound'],1)
        with self.assertRaises(ValueError):singleton_certificate({'ainvs':self.a,'points':[['0','0']]})

    def test_torsion_coset_anchors_have_exact_provenance(self):
        pool={'ainvs':self.a,'points':[list(map(str,self.p))],
              'vectors':[[1]],'approximate_heights':[1.0]}
        result=translate_pool(pool,16)
        self.assertEqual(len(result['points']),7)
        for point,vector,shift in zip(result['points'],result['vectors'],result['torsion_translations']):
            expected=add(list(map(Q,self.a)),multiply(list(map(Q,self.a)),self.p,vector[0]),
                         None if shift is None else tuple(map(Q,shift)))
            self.assertEqual(tuple(map(Q,point)),expected)

    def test_reduced_only_model_keeps_the_exact_inverse(self):
        from models import prepare_models
        from seeded import batch_search
        doubled=multiply(list(map(Q,self.a)),self.p,2)
        pool={'ainvs':self.a,'points':[list(map(str,doubled))],'approximate_heights':[1.0]}
        models=prepare_models(pool,1,{},minimal=False)
        self.assertTrue(models);self.assertEqual(models[0]['reduction'],'reduced_only')
        # Polynomial identities are checked in GP; every emitted point is
        # checked again against the original equation by batch_search.
        result=batch_search(pool,models,256,16,1)
        self.assertTrue(result['observations'])

    def test_gp_exponent_normalization_is_limited_to_lattice(self):
        from seeded import normalize_lattice_output
        raw='LATTICE [[[1,0],[0,1]],[[1.0,-2.34 E-95],[-2.34 E-95,4.0]]]\nBASIS ["1", "2"]\nPOOL_END\n'
        result=normalize_lattice_output(raw)
        parsed=json.loads(result.splitlines()[0][8:])
        self.assertEqual(parsed[1][0][1],-2.34e-95)
        self.assertEqual(result.splitlines()[1:],raw.splitlines()[1:])


if __name__=='__main__':unittest.main()

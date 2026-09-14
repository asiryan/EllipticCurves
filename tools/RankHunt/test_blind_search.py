"""Small checks of the seeded search, data boundary and exact independence."""
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from blind_search import ROOT, clean, certify, boxes, order_models
from point_arithmetic import add
from point_search import save
from bounded_anchor_pool import bounded_vectors
from anchor_diversity import diverse_vectors


class BlindSearchTests(unittest.TestCase):
    def setUp(self):
        parent=ROOT/'artifacts/blind-search-tests';parent.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=parent);self.root=Path(self.temp.name)

    def tearDown(self):
        if not self.root.resolve().is_relative_to((ROOT/'artifacts/blind-search-tests').resolve()):raise RuntimeError('Bad test path')
        self.temp.cleanup()

    def test_exact_basis_drops_dependent_multiple(self):
        a=list(map(Q,[0,0,0,-25,4]));p=(Q(0),Q(2));twice=add(a,p,p)
        data={'ainvs':list(map(str,a)),'points':[list(map(str,p)),list(map(str,twice))]}
        source=self.root/'input.json';save(source,data)
        result=certify(source)
        self.assertTrue(result['all_selected_independent'])
        self.assertEqual(result['LowerBound'],1)
        self.assertEqual(len(result['points']),1)

    def test_lattice_enumeration_has_a_real_candidate_limit(self):
        gram=[[4.0 if i==j else .2 for j in range(24)] for i in range(24)]
        rows,stats=bounded_vectors(gram,200,50)
        self.assertEqual(stats['vectors_considered'],200)
        self.assertEqual(len(rows),50)
        for score,v in rows:
            exact=sum(v[i]*gram[i][j]*v[j] for i in range(24) for j in range(24))
            self.assertAlmostEqual(score,exact)
            self.assertLessEqual(max(map(abs,v)),3)
            self.assertLessEqual(sum(map(abs,v)),12)

    def test_two_torsion_not_promoted_to_independent_basis(self):
        source=self.root/'torsion.json'
        save(source,{'ainvs':['0','0','0','-1','0'],'points':[['0','0']]})
        result=certify(source)
        self.assertEqual(result['LowerBound'],0)
        self.assertFalse(result['all_selected_independent'])
        self.assertEqual(result['points'],[])

    def test_diverse_vectors_are_bounded_reproducible_and_include_wider_support(self):
        gram=[[4.0 if i==j else .2 for j in range(24)] for i in range(24)]
        rows,stats=diverse_vectors(gram,1024,64)
        self.assertEqual((rows,stats),diverse_vectors(gram,1024,64))
        self.assertEqual(len(rows),64)
        self.assertEqual(len({v for _,v in rows}),64)
        self.assertLessEqual(stats['vectors_considered'],1024)
        self.assertLessEqual(stats['random_attempts'],stats['random_attempt_limit'])
        self.assertTrue(any(sum(x!=0 for x in v)==8 for _,v in rows))
        for score,v in rows:
            self.assertAlmostEqual(score,sum(v[i]*gram[i][j]*v[j] for i in range(24) for j in range(24)))

    def test_untrusted_metadata_is_not_input(self):
        result=clean({'ainvs':['0','0','0','-25','4'],'points':[['0','2']],
            'reference':'somewhere.json','anchor_source':'other.json','mode':'rank'})
        self.assertEqual(set(result),{'ainvs','points'})
        with self.assertRaises(ValueError):clean({'ainvs':['0','0','0','-25','4'],'points':[['1','1']]})

    def test_bounded_search_finds_new_independent_point_without_reference(self):
        source=self.root/'input.json'
        save(source,{'ainvs':['0','0','0','-25','4'],'points':[['0','2']]})
        out=self.root/'run'
        args=[sys.executable,str(ROOT/'tools/RankHunt/blind_search.py'),'--input',str(source),
            '--output',str(out),'--seconds','10','--anchors','16','--workers','1','--target','2']
        p=subprocess.run(args,capture_output=True,text=True,timeout=20)
        self.assertEqual(p.returncode,0,p.stdout+'\n'+p.stderr)
        summary=json.loads((out/'summary.json').read_text())
        self.assertGreaterEqual(summary['lower_bound'],2)
        self.assertEqual(summary['status'],'target_reached')
        self.assertFalse(summary['reference_points_loaded'])
        config=json.loads((out/'config.json').read_text())
        before=summary['total_jobs']
        p=subprocess.run(args+['--resume'],capture_output=True,text=True,timeout=20)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(json.loads((out/'summary.json').read_text())['total_jobs'],before)
        self.assertEqual(json.loads((out/'config.json').read_text()),config)
        from audit_blind_search import audit
        report=audit(out,self.root/'audited')
        self.assertTrue(report['archive_is_exact_seed_union_search_outputs'])
        self.assertGreaterEqual(report['certified_lower_bound'],2)

    def test_workspace_reference_read_is_blocked(self):
        source=self.root/'input.json';save(source,{})
        secret=self.root/'reference.json';save(secret,{'secret':1})
        out=self.root/'guard';out.mkdir()
        script=('import sys;from pathlib import Path;'
            'sys.path.insert(0,sys.argv[1]);from blind_search import file_guard;'
            'file_guard(Path(sys.argv[2]),Path(sys.argv[3]));Path(sys.argv[4]).read_text()')
        p=subprocess.run([sys.executable,'-c',script,str(ROOT/'tools/RankHunt'),str(source),str(out),str(secret)],
            capture_output=True,text=True,timeout=10)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('Unexpected workspace data access',p.stderr)

    def test_audit_requires_an_exact_unimodular_change(self):
        from audit_blind_search import integer_inverse
        self.assertEqual(integer_inverse([[1,3],[0,1]]),[[1,-3],[0,1]])
        self.assertEqual(integer_inverse([[0,1],[1,0]]),[[0,1],[1,0]])
        for bad in ([[2,0],[0,1]],[[1,1],[1,1]]):
            with self.assertRaises(ValueError):integer_inverse(bad)

    def test_batched_search_matches_individual_models(self):
        from fast_search import batch_search
        from blind_search import search_job
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2']]}
        models=[{'key':str(i),'anchor':p,'pool_index':i+1}
                for i,p in enumerate([['0','2'],['-5','2'],['5','2']])]
        expected=[]
        for i,m in enumerate(models):
            expected.extend(search_job(data,m,128,16,self.root/str(i),10)['points'])
        actual=batch_search(data,models,128,16,10)
        self.assertEqual(len(actual['complete']),len(models))
        self.assertEqual(clean({**data,'points':expected}),
            clean({**data,'points':[o['point'] for o in actual['observations']]}))

    def test_batch_timeout_preserves_points_and_only_completed_models(self):
        from fast_search import batch_search
        from unittest.mock import patch
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2']]}
        models=[{'key':str(i),'anchor':['0','2']} for i in range(3)]
        timeout=subprocess.TimeoutExpired('gp',1,output=b'ANCHOR_BEGIN 1\nPOINT [0, 2]\nANCHOR_DONE 1 1\nANCHOR_BEGIN 2\n')
        with patch('fast_search.subprocess.run',side_effect=timeout):
            result=batch_search(data,models,128,16,1)
        self.assertEqual(result['complete'],['0:128:16'])
        self.assertEqual(result['observations'][0]['point'],['0','-2'])
        self.assertTrue(result['timed_out'])

    def test_fast_search_finds_independent_points_and_resumes(self):
        from fast_search import search,independent_result
        source=self.root/'fast-input.json'
        save(source,{'ainvs':['0','0','0','-25','4'],'points':[['0','2']]})
        out=self.root/'fast'
        found=search(source,out,10,1,16,2)
        self.assertGreaterEqual(found['rank_lower_bound'],2)
        proof=independent_result(found,found['rank_lower_bound'])
        self.assertTrue(proof['verification']['all_points_independent_modulo_torsion'])
        jobs=json.loads((out/'checkpoint.json').read_text())['attempted_models']
        search(source,out,1,1,16,2)
        self.assertEqual(json.loads((out/'checkpoint.json').read_text())['attempted_models'],jobs)
        self.assertLessEqual(len(list(out.rglob('*'))),5)


if __name__=='__main__':unittest.main()

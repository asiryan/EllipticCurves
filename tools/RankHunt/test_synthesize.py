"""Exact identities, birational maps and independent checks of constructed points."""
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from synthesize import ROOT, canonical_roots, construct, evaluate, identity, model, run, save, select_basis, trial, verifier
from point_arithmetic import on_curve


class SynthesisTests(unittest.TestCase):
    def test_identity_and_quartic_inverse(self):
        roots = [0,1,2,4,10,12,16,18,28,32]
        roots, g, f, scale = identity(roots)
        # Equality at > degree many different inputs checks the entire identity.
        for x in range(-4,10):
            expected = scale*scale
            for root in roots: expected *= x-root
            self.assertEqual(evaluate(g,x)**2-evaluate(f,x),expected)
        curve, mapping = model(f,g)
        e,d,c,b,a = f; q = g[0]
        self.assertIsNone(mapping(0,q))
        self.assertTrue(on_curve(curve,mapping(0,-q)))
        for x in roots[1:]:
            for y in (evaluate(g,x),-evaluate(g,x)):
                X,Y = mapping(x,y)
                self.assertTrue(on_curve(curve,(X,Y)))
                if X*X != 4*a*e:
                    back = (2*q*Y+d*X+2*e*b)/(X*X-4*a*e)
                    self.assertEqual(back,x)
                    self.assertEqual((X*back*back-d*back-2*e)/(2*q),y)
        # Independent classical binary-quartic invariants give the same j.
        I = 12*a*e-3*b*d+c*c
        J = 72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
        expected_j = Q(6912*I**3,4*I**3-J*J)
        self.assertEqual(Q(construct(roots,0)['j_invariant']),expected_j)

    def test_affine_root_duplicates_and_degenerate_outputs(self):
        roots = [0,3,9,11,15,17,27,31]
        self.assertEqual(canonical_roots(roots),canonical_roots([13-7*r for r in roots]))
        with self.assertRaises(ValueError): construct(list(range(8)),0)
        with self.assertRaises(ValueError): construct([1]*10,0)

    def test_certified_points_are_observed_and_independently_verified(self):
        dll = ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'
        parent = ROOT/'artifacts/synthesis-tests'; parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as folder:
            for count, bound in [(8,7),(10,9)]:
                result = trial((0,20260917,32,count,128))
                self.assertEqual(result['rank_lower_bound'],bound)
                self.assertEqual(result,trial((0,20260917,32,count,128)))
                construction = result['construction']
                rebuilt = construct(construction['roots'],128)
                self.assertTrue(set(map(tuple,result['points'])) <= set(map(tuple,rebuilt['points'])))
                self.assertTrue(verifier().verify_certificate(result,result['certificate'])['all_points_independent_modulo_torsion'])
                path = Path(folder)/f'curve{count}.json'; path.write_text(json.dumps(result))
                process = subprocess.run(['dotnet',str(dll),'verify','--input',str(path)],
                    capture_output=True,text=True,timeout=20)
                self.assertEqual(process.returncode,0,process.stderr)
                self.assertEqual(json.loads(process.stdout)['LowerBound'],bound)
        self.assertIsNone(select_basis({'ainvs':['0','0','0','-1','0'],'points':[['0','0']]}))

    def test_bounded_run_resumes_without_new_work(self):
        parent = ROOT/'artifacts/synthesis-tests'; parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as folder:
            output = Path(folder)/'run'; result = Path(folder)/'best.json'
            command = [sys.executable,str(ROOT/'tools/RankHunt/synthesize.py'),'--samples','8',
                '--workers','2','--scan','16','--search-count','0','--seconds','0',
                '--output',str(output),'--result',str(result)]
            process = subprocess.run(command,capture_output=True,text=True,timeout=30)
            self.assertEqual(process.returncode,0,process.stderr)
            before = (output/'pool.json').read_bytes(); best = result.read_bytes()
            process = subprocess.run(command,capture_output=True,text=True,timeout=30)
            self.assertEqual(process.returncode,0,process.stderr)
            self.assertEqual((output/'pool.json').read_bytes(),before)
            self.assertEqual(result.read_bytes(),best)
            changed = command.copy(); changed[changed.index('--samples')+1] = '9'
            process = subprocess.run(changed,capture_output=True,text=True,timeout=30)
            self.assertNotEqual(process.returncode,0)
            self.assertIn('Synthesis settings changed',process.stderr)

    def test_mutation_rebuilds_valid_points_on_its_own_curve(self):
        parent = [0,3,4,9,12,18,23,24,27,32]
        result = trial((3,20260917,48,10,64,parent))
        self.assertEqual(result['status'],'certified')
        c = result['construction']
        self.assertNotEqual(list(c['roots']),parent)
        self.assertEqual(c['mutation_parent_roots'],parent)
        rebuilt = construct(c['roots'],64)
        self.assertEqual(result['ainvs'],rebuilt['ainvs'])
        self.assertTrue(set(map(tuple,result['points'])) <= set(map(tuple,rebuilt['points'])))

    def test_parallel_search_receives_prepared_directories_and_exact_seeds(self):
        parent = ROOT/'artifacts/synthesis-tests'; parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as folder:
            args = SimpleNamespace(output=str(Path(folder)/'run'),result=str(Path(folder)/'best.json'),
                samples=8,heights=[24],roots=[10],scan=16,seed=20260917,keep=4,workers=2,
                adaptive=False,search_count=4,seconds=1,point_workers=1,parallel_curves=4,anchors=16)
            seen=[]
            def search(source,destination,seconds,*unused,**kwargs):
                self.assertTrue(destination.is_dir())
                data=json.loads(source.read_text()); seen.append(data['trial'])
                save(destination/'checkpoint.json',{'wall_seconds':seconds})
                return {'ainvs':data['ainvs'],'points':data['points'],'rank_lower_bound':len(data['points'])}
            with patch('fast_search.search',side_effect=search): run(args)
            self.assertEqual(len(seen),4)
            # A stale or modified file must never replace the generated seed curve.
            path=next((Path(args.output)/'inputs').glob('*.json'))
            data=json.loads(path.read_text()); data['construction']['roots'][0] = -999
            save(path,data)
            with self.assertRaisesRegex(ValueError,'Candidate input differs'): run(args)


if __name__ == '__main__': unittest.main()

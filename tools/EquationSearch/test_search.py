"""Regressions for equation-only input, exact maps, isolation and the full pipeline."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools/RankHunt'))
import bootstrap as b


class EquationSearchTests(unittest.TestCase):
    def setUp(self):
        directory=ROOT/'artifacts/equation-search-tests';directory.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=directory);self.folder=Path(self.temp.name)

    def tearDown(self):
        if not self.folder.resolve().is_relative_to((ROOT/'artifacts/equation-search-tests').resolve()):
            raise RuntimeError('Test cleanup escaped its temporary workspace')
        self.temp.cleanup()

    def test_input_rejects_points_metadata_and_singular_equations(self):
        a=['0','0','0','-25','4']
        self.assertEqual(b.equation({'ainvs':a}),{'ainvs':a})
        for data in ({'ainvs':a,'points':[]},{'ainvs':a,'family':'anything'},
                     {'ainvs':['0','0','0','0','0']},{'ainvs':['0','0','0','1/2','1']}):
            with self.assertRaises(ValueError):b.equation(data)
        for call in ('ellrank(E)','ell2cover(E)','read("points.json")','ellsearch(11)'):
            with self.assertRaises(ValueError):b.checked_script(call)

    def test_nontrivial_quartic_map_and_finite_points_at_infinity(self):
        data={'ainvs':['0','0','0','-25','4']}
        prepared,_=b.prepare(data,5)
        out,_=b.gp(b.model_script(data,prepared,3,7),5)
        reduced=json.loads(next(line[6:] for line in out.splitlines() if line.startswith('MODEL ')))
        self.assertNotEqual(reduced[2][1][2],'0')
        model={'id':0,'center':'3','stride':'7','reduced':reduced}
        out,_=b.gp(b.search_script(data,prepared,model,256,16),5)
        points=b.parse_points(out,list(map(b.Q,data['ainvs'])))
        self.assertTrue({('0','-2'),('5','-2'),('-5','-2')}<=set(points))
        self.assertGreaterEqual(b.select_basis({**data,'points':points})['rank_lower_bound'],2)

    def test_first_hit_only_when_it_cannot_be_the_original_infinity(self):
        data={'ainvs':['0','0','0','-25','4']}
        prepared,_=b.prepare(data,5)
        affine={'id':0,'center':'0','stride':'1',
                'reduced':[['4','-25','0','1','0'],['0','0','0'],
                           ['1',['1','0','0','1'],['0','0','0']]]}
        script=b.search_script(data,prepared,affine,256,16)
        self.assertIn('hyperellratpoints(C,[256,16],1)',script)
        out,_=b.gp(script,5)
        self.assertEqual(len(b.parse_points(out,list(map(b.Q,data['ainvs'])))),1)
        out,_=b.gp(b.model_script(data,prepared,3,7),5)
        affine['reduced']=json.loads(next(line[6:] for line in out.splitlines() if line.startswith('MODEL ')))
        self.assertIn('hyperellratpoints(C,[256,16],0)',b.search_script(data,prepared,affine,256,16))

    def test_large_equation_without_seed_reaches_target_and_blocks_other_data(self):
        u=10**12;r=10**27+17
        a=[0,-3*r,0,3*r*r-25*u**4,4*u**6+25*u**4*r-r**3]
        source=self.folder/'curve.json';b.save(source,{'ainvs':list(map(str,a))})
        out=self.folder/'run'
        process=subprocess.run([sys.executable,str(HERE/'run.py'),'--input',str(source),'--output',str(out),
            '--bootstrap-seconds','10','--search-seconds','10','--target','2','--workers','2','--anchors','16','--seed-limit','1'],
            text=True,capture_output=True,timeout=35)
        self.assertEqual(process.returncode,0,process.stdout+'\n'+process.stderr)
        result=json.loads((out/'result.json').read_text())
        self.assertTrue(result['target_reached']);self.assertTrue(result['equation_only'])
        self.assertEqual(result['initial_points'],0)
        self.assertTrue(result['verification']['all_points_independent_modulo_torsion'])
        self.assertEqual(result['bootstrap_errors'],[])
        self.assertEqual(len(json.loads((out/'seed.json').read_text())['points']),1)
        secret=self.folder/'reference.json';b.save(secret,{'points':[['0','2']]})
        script=('import sys;from pathlib import Path;sys.path.insert(0,sys.argv[1]);'
                'from run import data_boundary;data_boundary(Path(sys.argv[2]),Path(sys.argv[3]));'
                'Path(sys.argv[4]).read_text()')
        denied=subprocess.run([sys.executable,'-c',script,str(HERE),str(source),str(out),str(secret)],
                              text=True,capture_output=True,timeout=5)
        self.assertNotEqual(denied.returncode,0);self.assertIn('Equation-only data boundary',denied.stderr)

    def test_local_expansion_engine_finds_and_verifies_an_extra_point(self):
        from seeded import search,independent_result
        source=self.folder/'seed.json'
        b.save(source,{'ainvs':['0','0','0','-25','4'],'points':[['0','2']]})
        result=search(source,self.folder/'expand',10,1,16,2)
        self.assertGreaterEqual(result['rank_lower_bound'],2)
        verified=independent_result(result,result['rank_lower_bound'])
        self.assertTrue(verified['verification']['all_points_independent_modulo_torsion'])

    def test_projection_inverse_and_companion_group_relation(self):
        from geometry import to_quartic,from_quartic
        from point_arithmetic import add,multiply,negate
        a=list(map(b.Q,[0,0,0,-25,4]));p=(b.Q(0),b.Q(2));q=(b.Q(5),b.Q(2))
        checked=0
        for i in range(-3,4):
            for j in range(-3,4):
                candidate=add(a,multiply(a,p,i),multiply(a,q,j))
                if candidate is None or candidate[0]==p[0]:continue
                t,z=to_quartic(a,p,candidate)
                self.assertEqual(from_quartic(a,p,(t,z)),candidate)
                self.assertEqual(from_quartic(a,p,(t,-z)),negate(a,add(a,p,candidate)))
                checked+=1
        self.assertGreater(checked,30)

    def test_denominator_square_search_reaches_nonintegral_points(self):
        from geometry import small_points
        from point_arithmetic import multiply,negate
        a=list(map(b.Q,[0,0,0,-25,4]));twice=multiply(a,(b.Q(0),b.Q(2)),2)
        observed={p for p,_ in small_points(a,700,4)}
        self.assertIn(min(twice,negate(a,twice)),observed)
        self.assertTrue(all(b.on_curve(a,p) for p in observed))

    def test_fixed_single_anchor_finds_new_direction(self):
        from seeded import search,independent_result
        source=self.folder/'single.json';b.save(source,{'ainvs':['0','0','0','-25','4'],'points':[['0','2']]})
        result=search(source,self.folder/'fixed',5,1,16,2,anchor_mode='fixed')
        self.assertGreaterEqual(result['rank_lower_bound'],2)
        self.assertTrue(independent_result(result,result['rank_lower_bound'])['verification']['all_points_independent_modulo_torsion'])


if __name__=='__main__':unittest.main()

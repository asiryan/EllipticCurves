import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from point_search import ROOT,save
from record_hunt import candidate_order,promoted,j_invariant,check_engine


class RecordHuntTests(unittest.TestCase):
    def setUp(self):
        parent=ROOT/'artifacts/record-hunt-tests';parent.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=parent);self.root=Path(self.temp.name)

    def tearDown(self):
        if not self.root.resolve().is_relative_to((ROOT/'artifacts/record-hunt-tests').resolve()):raise RuntimeError('Bad test path')
        self.temp.cleanup()

    def test_j_ignores_scaled_equation_and_rejects_singular_curve(self):
        self.assertEqual(j_invariant({'ainvs':[0,0,0,-1,0]}),1728)
        self.assertEqual(j_invariant({'ainvs':[0,0,0,-16,0]}),1728)
        with self.assertRaises(ValueError):j_invariant({'ainvs':[0,0,0,0,0]})

    def test_allocation_uses_certificates_before_scores(self):
        rows=[{'id':'a','score':100,'tail_score':10,'lower_bound':18,'status':'budget_completed'},
              {'id':'b','score':1,'tail_score':1,'lower_bound':24,'status':'budget_completed'},
              {'id':'c','score':200,'tail_score':100,'lower_bound':25,'status':'error'}]
        self.assertEqual(promoted(rows,1)[0]['id'],'b')
        ordered=candidate_order(rows,3)
        self.assertEqual(len({r['id'] for r in ordered}),3)
        self.assertEqual(ordered,candidate_order(list(reversed(rows)),3))

    def test_preserved_engine_is_unchanged(self):
        self.assertEqual(len(check_engine()),5)

    def test_multiscale_generation_counts_deduplicates_and_resumes(self):
        from record_hunt import prepare_multiscale
        calls=[]
        def sample(command,**kwargs):
            count=int(command[command.index('--samples')+1]);calls.append(count)
            folder=Path(command[command.index('--output')+1]);folder.mkdir()
            save(folder/'candidates.json',[{'file':'curve.json','j_invariant':'1728','u':1,'v':1}])
            save(folder/'complete.json',{'PrimitiveDraws':count,'exported':1,'seconds':.1})
            return SimpleNamespace(returncode=0,stdout='',stderr='')
        out=self.root/'multiscale'
        with patch('record_hunt.subprocess.run',side_effect=sample):
            prepare_multiscale(out,11,[100,200],1,10)
            prepare_multiscale(out,11,[100,200],1,10)
            with self.assertRaises(ValueError):prepare_multiscale(out,12,[100,200],1,10)
        self.assertEqual(sorted(calls),[5,6])
        meta=json.loads((out/'generation.json').read_text())
        self.assertEqual(meta['primitive_draws'],11)
        rows=json.loads((out/'candidates.json').read_text())
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['file'],'h100/curve.json')

    def test_crt_sampler_repeats_and_keeps_valid_specializations(self):
        dll=os.environ.get('RANKHUNT_TEST_DLL',str(ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'))
        outputs=[]
        for name in ['crt-a','crt-b']:
            folder=self.root/name
            command=['dotnet',dll,'sample','--selection','crt','--samples','100','--height','1000',
                '--keep','8','--refine-keep','8','--final-keep','2','--prime-bound','16382',
                '--workers','2','--seed','101','--output',str(folder)]
            process=subprocess.run(command,capture_output=True,text=True,timeout=45)
            self.assertEqual(process.returncode,0,process.stdout+process.stderr)
            meta=json.loads((folder/'complete.json').read_text())
            self.assertTrue(meta['congruence_sampling'])
            self.assertGreater(meta['CongruencePrimitiveDraws'],0)
            rows=json.loads((folder/'candidates.json').read_text());outputs.append(rows)
            from blind_search import clean,certify
            for row in rows:
                path=folder/row['file'];data=clean(json.loads(path.read_text()))
                self.assertEqual(str(j_invariant(data)),row['j_invariant'])
                self.assertEqual(certify(path)['LowerBound'],row['seed_lower_bound'])
        self.assertEqual(outputs[0],outputs[1])

    def test_adaptive_selection_keeps_best_and_reserves_exploration(self):
        from record_hunt import adaptive_promoted
        rows=[{'id':str(i),'lower_bound':25-i,'tail_score':float(i),
               'score':30-i,'status':'budget_completed','quartic_bits_p10':100-i}
              for i in range(8)]
        chosen=adaptive_promoted(rows,3)
        self.assertEqual(chosen[0]['id'],'0')
        self.assertIn('7',{r['id'] for r in chosen})
        rows[-1]['status']='error'
        self.assertNotIn('7',{r['id'] for r in adaptive_promoted(rows,3)})

    def test_later_prime_confirmation_overrides_early_tail_noise(self):
        from record_hunt import adaptive_promoted,arithmetic_score
        rows=[{'id':'best','lower_bound':24,'score':20,'tail_score':2,'confirmation_score':2.8,'status':'budget_completed'},
              {'id':'noise','lower_bound':17,'score':20,'tail_score':4,'confirmation_score':1.4,'status':'budget_completed'},
              {'id':'steady','lower_bound':17,'score':19,'tail_score':3,'confirmation_score':2.7,'status':'budget_completed'}]
        self.assertEqual([r['id'] for r in adaptive_promoted(rows,2)],['best','steady'])
        # Independent direct enumeration for the exact same good-reduction primes.
        import math
        data={'ainvs':['0','0','0','-25','4']};expected=0
        for p in range(5,100):
            if any(p%d==0 for d in range(2,math.isqrt(p)+1)):continue
            if (-16*(4*(-25)**3+27*4**2))%p==0:continue
            count=1+sum((y*y-x*x*x+25*x-4)%p==0 for x in range(p) for y in range(p))
            expected+=math.log(count/p)
        self.assertAlmostEqual(arithmetic_score(data,3,99),expected,places=12)

    def test_compact_campaign_certifies_and_completed_resume_is_idle(self):
        pool=self.root/'compact-pool';pool.mkdir()
        save(pool/'curve.json',{'ainvs':['0','0','0','-25','4'],'points':[['0','2']]})
        save(pool/'candidates.json',[{'id':'toy','u':0,'v':1,'file':'curve.json',
                                     'score':0,'tail_score':0}])
        out=self.root/'compact-run'
        command=[sys.executable,str(ROOT/'tools/RankHunt/record_hunt.py'),'--fast',
            '--candidates',str(pool),'--output',str(out),'--counts','1','--seconds','5',
            '--target','2','--point-workers','1','--parallel-curves','1','--anchors','16']
        p=subprocess.run(command,capture_output=True,text=True,timeout=20)
        self.assertEqual(p.returncode,0,p.stderr)
        result=json.loads((out/'best.json').read_text())
        self.assertGreaterEqual(result['rank_lower_bound'],2)
        before=(out/'campaign.json').read_bytes()
        p=subprocess.run(command,capture_output=True,text=True,timeout=20)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(before,(out/'campaign.json').read_bytes())
        self.assertLessEqual(len([p for p in out.rglob('*') if p.is_file()]),6)

    def test_live_monitor_does_not_open_the_workers_replaceable_checkpoint(self):
        import record_hunt
        out=self.root/'monitor';(out/'processes').mkdir(parents=True)
        folder=out/'searches/c000';folder.mkdir(parents=True)
        summary=folder/'summary.json';active=[False]
        real_read=record_hunt.read
        def guarded_read(path):
            if Path(path)==summary and active[0]:raise PermissionError('simulated Windows sharing violation')
            return real_read(path)
        class Process:
            pid=12345;returncode=0
            def __init__(self,*args,stdout,**kwargs):
                self.polls=0;active[0]=True
                save(summary,{'wall_seconds':2,'lower_bound':2,'status':'budget_completed'})
                stdout.write('{"lower_bound":2,"seconds":1}\n');stdout.flush()
            def poll(self):
                self.polls+=1
                if self.polls<4:return None
                active[0]=False;return 0
        args=SimpleNamespace(anchors=16,point_workers=1,job_seconds=1,target=3)
        row={'id':'c000','lower_bound':1}
        with patch.object(record_hunt.subprocess,'Popen',Process),patch.object(record_hunt,'read',guarded_read),\
             patch.object(record_hunt.time,'sleep'),patch.object(record_hunt,'kill_worker'):
            result=record_hunt.search_one(row,out,5,args,0)
        self.assertEqual(result['lower_bound'],2)
        self.assertEqual(result['status'],'budget_completed')

    def test_campaign_finds_points_audits_them_and_completed_resume_is_idle(self):
        pool=self.root/'pool';pool.mkdir()
        data={'ainvs':['0','0','0','-25','4'],'points':[['0','2']]}
        save(pool/'curve.json',data)
        save(pool/'candidates.json',[{'id':'toy','file':'curve.json','score':0,'tail_score':0,
                                     'j_invariant':str(j_invariant(data))}])
        out=self.root/'run'
        command=[sys.executable,str(ROOT/'tools/RankHunt/record_hunt.py'),'--candidates',str(pool),
            '--output',str(out),'--counts','1','1','--seconds','5','10','--target','2',
            '--point-workers','1','--parallel-curves','1','--anchors','16']
        p=subprocess.run(command,capture_output=True,text=True,timeout=30)
        self.assertEqual(p.returncode,0,p.stdout+'\n'+p.stderr)
        report=json.loads((out/'best-verified.json').read_text())
        self.assertGreaterEqual(report['certified_lower_bound'],2)
        self.assertTrue(report['different_j_from_reference'])
        self.assertEqual(json.loads((out/'summary.json').read_text())['status'],'target_found')
        count=len(list((out/'processes').glob('*.command.json')))
        p=subprocess.run(command+['--resume'],capture_output=True,text=True,timeout=30)
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(count,len(list((out/'processes').glob('*.command.json'))))
        changed=command.copy();changed[changed.index('--target')+1]='3'
        p=subprocess.run(changed+['--resume'],capture_output=True,text=True,timeout=30)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('unchanged inputs, code and policy',p.stderr)

    def test_sampler_is_repeatable_and_exported_seed_points_are_valid(self):
        outputs=[]
        command=['dotnet',str(ROOT/'tools/RankHunt/bin/Release/net8.0/RankHunt.dll'),'sample',
            '--samples','100','--height','1000','--keep','8','--refine-keep','4','--final-keep','2',
            '--prime-bound','16382','--workers','2','--seed','73']
        for name in ['one','two']:
            out=self.root/name
            p=subprocess.run(command+['--output',str(out)],capture_output=True,text=True,timeout=30)
            self.assertEqual(p.returncode,0,p.stdout+'\n'+p.stderr)
            rows=json.loads((out/'candidates.json').read_text());outputs.append(rows)
            self.assertEqual(len({r['j_invariant'] for r in rows}),len(rows))
            from blind_search import clean,certify
            for row in rows:
                data=clean(json.loads((out/row['file']).read_text()))
                self.assertEqual(str(j_invariant(data)),row['j_invariant'])
                self.assertEqual(certify(out/row['file'])['LowerBound'],row['seed_lower_bound'])
        self.assertEqual(outputs[0],outputs[1])


if __name__=='__main__':unittest.main()

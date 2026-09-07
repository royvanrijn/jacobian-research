#!/usr/bin/env python3
"""Sequential focused intake, isolated seed replay, point exposure and proofs."""
import sys,shutil
from pathlib import Path
import curve302_focused_point_exposure_v2 as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
cert=batch.cert;D=batch.BATCH

def main():
    assert not (D/'run-ledger.json').exists()
    ledger=dict(status='RUNNING',stages=[],proofs=[]);checkpoint(D/'run-ledger.json',ledger)
    checkpoint(D/'run-protocol.json',dict(source_sha256=cert.hashed(Path(__file__)),
        maximum_workers=1,rss_bytes=2147483648,intake_seconds=120,seed_replay_seconds=120,
        geometry_seconds_per_curve=180,worker_seconds_per_curve=1200,replay_seconds_per_curve=1200,
        proof_driver_seconds_per_curve=1800,maximum_curves=4,maximum_boxes=196))
    def stage(name,cmd,seconds,cwd=batch.ROOT):
        s=run(cmd,limits=Limits(seconds,2147483648),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=cwd)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(D/'run-ledger.json',ledger)
        print(name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'run-ledger.json',ledger);raise ArithmeticError('preserve failed focused stage')
    stage('prepare',[sys.executable,str(batch.CAS/'curve302_focused_point_exposure_v2.py'),'prepare'],120)
    fresh=D/'seed-replay';fresh.mkdir(exist_ok=False)
    for n in ['intake.json','generic-parent.json']:shutil.copy2(D/n,fresh/n)
    for n in ['verify_curve302_focused_seeds.sage','verify_factor_free_rank.sage']:shutil.copy2(batch.CAS/n,fresh/n)
    for row in cert.read(D/'intake.json')['rows']:
        for n,out in [('seed.json','-seed.json'),('seed-cloud.json','-cloud.json')]:shutil.copy2(D/row['id']/n,fresh/(row['id']+out))
    checkpoint(fresh/'protocol.json',dict(files={p.name:cert.hashed(p) for p in fresh.iterdir()},seconds=120,rss_bytes=2147483648,maximum_workers=1))
    stage('seed-replay',[batch.SAGE,str(fresh/'verify_curve302_focused_seeds.sage'),'--directory',str(fresh)],120,fresh)
    stage('freeze',[sys.executable,str(batch.CAS/'curve302_focused_point_exposure_v2.py'),'freeze'],60)
    batch.launch()
    for row in batch.protocol()['rows']:
        folder=D/row['id'];name='proof-'+row['id']
        stage(name,[sys.executable,str(batch.CAS/'certify_factor_free_exposure_v3.py'),
            '--run',str(folder),'--protocol',str(D/'protocol.json'),
            '--prefix','curve302_focused_'+row['id'].replace('-','_')],1800)
        proof=cert.read(folder/'certification-ledger.json')
        assert proof['status']=='PASS' and all(v==proof['rank_lower_bound'] for v in proof['odd_modulus_ranks'].values())
        ledger['proofs'].append(dict(id=row['id'],rank_lower_bound=proof['rank_lower_bound'],initial_rank=row['initial_rank'],
            completed_boxes=proof['completed_boxes'],certification_sha256=cert.hashed(folder/'certification-ledger.json')))
        checkpoint(D/'run-ledger.json',ledger)
    ledger['status']='PASS';checkpoint(D/'run-ledger.json',ledger);print('PASS196 boxes and four independently certified clouds',flush=True)

if __name__=='__main__':main()

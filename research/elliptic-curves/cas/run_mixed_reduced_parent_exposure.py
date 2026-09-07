#!/usr/bin/env python3
"""Sequential bounded intake, full map freeze, exposure and independent certificates."""
import sys,shutil
from pathlib import Path
import mixed_reduced_parent_exposure as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
cert=batch.cert;D=batch.BATCH
def main():
 assert not (D/'run-ledger.json').exists();ledger=dict(status='RUNNING',stages=[]);checkpoint(D/'run-ledger.json',ledger)
 checkpoint(D/'run-protocol.json',dict(driver_sha256=cert.hashed(Path(__file__)),maximum_workers=1,rss_bytes=2147483648,intake_seconds=180,seed_replay_seconds=120,certificate_curve_seconds=1800))
 def stage(name,cmd,seconds,cwd=batch.ROOT):
  s=run(cmd,limits=Limits(seconds,2147483648),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=cwd);ok=s['outcome']=='completed' and s['returncode']==0;ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(D/'run-ledger.json',ledger);print(name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'run-ledger.json',ledger);raise ArithmeticError('preserve failed or censored stage')
 stage('prepare',[sys.executable,str(batch.CAS/'mixed_reduced_parent_exposure.py'),'prepare'],60)
 stage('intake',[batch.SAGE,str(batch.CAS/'prepare_mixed_reduced_parent_seeds.sage')],180)
 fresh=D/'seed-replay';fresh.mkdir(exist_ok=False)
 for name in ['intake-input.json','intake-result.json']:shutil.copy2(D/name,fresh/name)
 shutil.copy2(batch.CAS/'verify_mixed_reduced_parent_seeds.sage',fresh/'verify.sage')
 for row in cert.read(D/'intake-result.json')['rows']:shutil.copy2(D/row['id']/'seed.json',fresh/(row['id']+'-seed.json'))
 checkpoint(fresh/'protocol.json',dict(files={p.name:cert.hashed(p) for p in fresh.iterdir()},seconds=120,rss_bytes=2147483648,maximum_workers=1))
 stage('seed-replay',[batch.SAGE,str(fresh/'verify.sage'),'--directory',str(fresh)],120,fresh)
 stage('freeze',[sys.executable,str(batch.CAS/'mixed_reduced_parent_exposure.py'),'freeze'],60)
 batch.launch();post=dict(status='RUNNING',rows=[]);checkpoint(D/'post-ledger.json',post)
 for row in batch.protocol()['rows']:
  folder=D/row['id'];cmd=[sys.executable,str(batch.CAS/'certify_factor_free_exposure_v3.py'),'--run',str(folder),'--protocol',str(D/'protocol.json'),'--prefix','mixed_reduced_parent_'+row['id'].replace('-','_')]
  s=run(cmd,limits=Limits(1800,2147483648),log_path=folder/'certification-driver.log',checkpoint_path=folder/'certification-driver.supervisor.json',cwd=batch.ROOT);ok=s['outcome']=='completed' and s['returncode']==0;post['rows'].append(dict(id=row['id'],status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(D/'post-ledger.json',post);print('PROOFS',row['id'],s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:post['status']=ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'post-ledger.json',post);checkpoint(D/'run-ledger.json',ledger);raise ArithmeticError('preserve failed or censored certificate')
 post['status']=ledger['status']='PASS';checkpoint(D/'post-ledger.json',post);checkpoint(D/'run-ledger.json',ledger);print('PASS all441 boxes and54 independent certificate stages',flush=True)
if __name__=='__main__':main()

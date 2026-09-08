#!/usr/bin/env python3
"""Run bounded post-search full-cloud, exact geometry and portable rank checks."""
import json,sys,shutil
from pathlib import Path
import retained26_gap_trial as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

def main():
 p=batch.protocol();d=batch.BATCH;out=d/'post-verification-ledger.json';assert cert.read(d/'ledger.json')['status']=='PASS' and not out.exists();ledger={'status':'RUNNING','stages':[]};checkpoint(out,ledger)
 jobs=[('cloud-audits',[sys.executable,str(batch.CAS/'audit_retained26_gap_trial.py')],600),('geometry',[batch.SAGE,str(batch.CAS/'replay_retained26_gap_geometry.sage')],120)]
 for name,cmd,seconds in jobs:
  s=run(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=d/(name+'.log'),checkpoint_path=d/(name+'.supervisor.json'),cwd=batch.ROOT);ok=s['outcome']=='completed' and s['returncode']==0;ledger['stages'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger);print(name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('post verification failed')
 fresh=d/'standalone';fresh.mkdir(exist_ok=False);files=[batch.CAS/'verify_retained26_gap_rank.sage']+[batch.ART/('retained26_gap_'+r['id'].replace('-','_')+'_mod2_v1.json') for r in p['rows']]
 for path in files:shutil.copy2(path,fresh/path.name)
 checkpoint(fresh/'protocol.json',{'files':{path.name:cert.hashed(path) for path in files},'seconds':120,'rss_bytes':p['rss_bytes'],'workers':1})
 cmd=[batch.SAGE,str(fresh/files[0].name),'--input',*[str(fresh/path.name) for path in files[1:]]]
 s=run(cmd,limits=Limits(120,p['rss_bytes']),log_path=fresh/'replay.log',checkpoint_path=fresh/'supervisor.json',cwd=fresh);ok=s['outcome']=='completed' and s['returncode']==0;ledger['stages'].append({'name':'independent-rank','status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});ledger['status']='PASS' if ok else 'FAILED_OR_CENSORED';checkpoint(out,ledger);print('independent-rank',s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
 if not ok:raise ArithmeticError('standalone rank proof failed')
 print('PASS all retained26 post checks',flush=True)
if __name__=='__main__':main()

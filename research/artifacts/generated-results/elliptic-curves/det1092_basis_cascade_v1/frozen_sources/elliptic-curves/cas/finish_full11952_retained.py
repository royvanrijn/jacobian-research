#!/usr/bin/env python3
"""One retained-bank, cached-score and64-finalist validation pipeline after the full scan."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_full11952_h131072 as parent
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=parent.ROOT;CAS=parent.CAS;D=ROOT/'artifacts/local/elliptic-curves/full11952-retained-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve one full11952 score controller')
    p=parent.protocol();paths=[Path(__file__).resolve(),CAS/'score_full11952_retained.py',parent.D/'protocol.json']
    if p['retained_rows']!=1048576 or len(p['rows'])!=2048:raise ArithmeticError('fixed complete scan required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-retained-controller.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_deadline_seconds':7800,'rss_bytes':4294967296,'jobs':[{'name':n,'seconds':s} for n,s in [('prepare',120),('bank',180),('bank-check',180),('short',300),('short-check',60),('extended',300),('extended-check',60),('select',180),('selection-check',180),('validate-selected',300),('validation-check',180)]],'scope':'Wait for the already frozen2048-slice full11952 short scan. Require every slice to complete within its original bounds. Verify the million-row retained bank and every short score, apply the unchanged cached extended score, freeze exactly64 distinct equations after excluding all789 previous equations, then require fresh scalar selection-score agreement and disjoint validation on all64. Stop on any failed/censored stage without retries, score changes, roster refills or point searches.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen retained controller changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve retained pipeline')
    ledger={'status':'WAITING_FOR_FULL_SCAN','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_deadline_seconds']
    while True:
        path=parent.D/'scan.supervisor.json'
        if path.exists():
            s=cert.read(path)
            if s['outcome']!='running':
                if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('full scan failed/censored; no score pipeline')
                break
        if time.monotonic()>deadline:raise TimeoutError('fixed full-scan wait deadline')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
    for job in p['jobs']:
        name=job['name'];s=run([sys.executable,str(CAS/'score_full11952_retained.py'),name],limits=Limits(job['seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('FULL11952 RETAINED',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('fixed retained stage failed/censored; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','launch']);v=a.parse_args();globals()[v.stage]()

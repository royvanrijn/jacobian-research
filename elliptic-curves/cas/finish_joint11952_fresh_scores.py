#!/usr/bin/env python3
"""Wait for exact full-score slices and terminal prior equations, then select64."""
import argparse,sys,time
from pathlib import Path
import score_joint11952_fresh as scores
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=scores.ROOT;CAS=scores.CAS;D=scores.D/'controller'
UPSTREAM=[scores.scan.D/'controller',ROOT/'artifacts/local/elliptic-curves/fresh60-mw16-pari-v1/post-batch',ROOT/'artifacts/local/elliptic-curves/annulus64-r17-pari-v2/post-batch']
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve joint score controller')
    paths=[Path(__file__).resolve(),Path(scores.__file__),Path(scores.engine.__file__),Path(scores.models.__file__),Path(scores.scoring.__file__),*(p/'protocol.json' for p in UPSTREAM),CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_seconds':7200,'rss_bytes':2147483648,'stages':[('prepare',120),('freeze',120),('run',1800),('check',180),('select',60),('selection-check',120),('validate',300),('validation-check',180)],'scope':'Wait for all eight full-score slice proofs and both terminal earlier curve cohorts before freezing the1101-equation exclusion snapshot. No prior point outcome controls the new slice definition or score weights. Bind32768 models, freeze4096 new scalar candidates, pass a first8 runtime gate, replay scores and exclusions, freeze64 and validate on disjoint primes. No point search or automatic retry.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('joint controller sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve joint score launch')
    ledger={'status':'WAITING_FOR_UPSTREAM_PROOFS','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            values=[cert.read(q/'ledger.json')['status'] for q in UPSTREAM]
            if all(v=='PASS' for v in values):break
            if any(v not in ('PASS','RUNNING','WAITING_FOR_POINT_PROOFS') for v in values) or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait elapsed')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,seconds in p['stages']:
            s=run([sys.executable,str(Path(scores.__file__)),name],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('JOINT11952 SCORE CONTROLLER',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('joint score stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

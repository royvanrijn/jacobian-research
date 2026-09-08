#!/usr/bin/env python3
"""Population-gated MW16 outer-band scoring and disjoint finalist validation."""
import argparse,sys,time
from pathlib import Path
import score_mw16_fresh_retention as scores
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=scores.ROOT;D=scores.D/'controller';CAS=scores.CAS

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve outer score controller')
    paths=[Path(__file__).resolve(),Path(scores.__file__),Path(scores.engine.__file__),Path(scores.models.__file__),Path(scores.scoring.__file__),scores.scan.D/'controller/protocol.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_seconds':1200,'stages':[('prepare',120),('run',3600),('check',600),('select',60),('selection_check',60),('validate',300),('validation-check',180)],'rss_bytes':2147483648,'scope':'Wait for all ten fresh disjoint outer slices and exact40960 short scores. Freeze40960 scalar traces, pass first10 runtime gate, use four workers and80-row checkpoints, replay all scores, select12 distinct curves per family and validate the frozen60 with disjoint primes. No point searches or automatic retries.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('score controller sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve launch')
    ledger={'status':'WAITING_FOR_OUTER_SCAN','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            upstream=cert.read(scores.scan.D/'controller/ledger.json')
            if upstream['status']=='PASS':break
            if upstream['status']!='RUNNING' or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait exhausted')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,seconds in p['stages']:
            s=run([sys.executable,str(Path(scores.__file__)),name],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('MW16 FRESH SCORE CONTROLLER',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('fixed score stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

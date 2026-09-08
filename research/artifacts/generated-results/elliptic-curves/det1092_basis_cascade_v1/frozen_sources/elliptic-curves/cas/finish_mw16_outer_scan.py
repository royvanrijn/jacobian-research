#!/usr/bin/env python3
"""Bounded controller for the fixed thirty-slice MW16 outer population."""
import argparse,sys
from pathlib import Path
import scan_mw16_outer_bands as scan
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
D=scan.D/'controller';ROOT=scan.ROOT;CAS=scan.CAS

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve outer scan controller')
    paths=[Path(__file__).resolve(),Path(scan.__file__),scan.CPP,CAS/'research_runtime/supervisor.py',CAS/'research_runtime/store.py']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'stages':[('prepare',120),('benchmark',180),('run',900),('replay',180)],'rss_bytes':2147483648,'scope':'One immutable scanner build and protocol, forty small signed/frame regression cases, at most thirty broad annular slices after three fixed runtime gates, then retained-score and population replay. No retries or point searches.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('controller sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve launch')
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    for name,seconds in p['stages']:
        supervision=run([sys.executable,str(Path(scan.__file__)),name],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=supervision['outcome']=='completed' and supervision['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':supervision})
        if not ok:ledger['status']='FAILED_OR_CENSORED'
        checkpoint(D/'ledger.json',ledger);print('OUTER MW16 CONTROLLER',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:raise ArithmeticError('preserved failed/censored stage; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

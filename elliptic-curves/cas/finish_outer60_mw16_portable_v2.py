#!/usr/bin/env python3
"""Wait for the fixed point proofs, then package and independently replay once."""
import argparse,sys,time
from pathlib import Path
import outer60_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=ROOT/'artifacts/local/elliptic-curves/outer60-mw16-point-portable-controller-v2'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve portable controller')
    paths=[Path(__file__).resolve(),CAS/'package_outer60_mw16_points_v2.py',CAS/'verify_outer60_mw16_points_portable_v2.py',batch.D/'post-batch/protocol.json']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.outer60-mw16-portable-controller.v2','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'wait_seconds':35000,'jobs':[{'name':'package','script':'package_outer60_mw16_points_v2.py','seconds':600},{'name':'isolated-replay','script':'verify_outer60_mw16_points_portable_v2.py','seconds':20000}],'rss_bytes':4294967296,'scope':'Wait for all fixed60 point searches, admission-history proofs and final exact geometry/mod2/mod3/mod5/equation comparisons to finish successfully. Require no unresolved odd-prime improvement. Then create one standalone point-only bundle and perform182 isolated replay stages. No scanner, score calculation, point search, adaptive follow-up or automatic retry.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('portable controller inputs changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve portable controller attempt')
    ledger={'status':'WAITING_FOR_POINT_FINALIZATION','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    while True:
        q=cert.read(batch.D/'post-batch/ledger.json')
        if q['status']=='PASS':break
        if q['status'] not in ('WAITING_FOR_POINT_PROOFS','RUNNING'):raise ArithmeticError('upstream point proof failed/censored or needs review')
        if time.monotonic()>deadline:raise TimeoutError('declared finalization wait elapsed')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
    for job in p['jobs']:
        n=job['name'];s=run([sys.executable,str(CAS/job['script'])],limits=Limits(job['seconds'],p['rss_bytes']),log_path=D/(n+'.log'),checkpoint_path=D/(n+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':n,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('OUTER60 PORTABLE',n,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('portable stage failed/censored; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','launch']);v=a.parse_args();globals()[v.stage]()

#!/usr/bin/env python3
"""Finish live finite MW16 cache proofs, encoding and saved-score regressions."""
import argparse,sys,time
from pathlib import Path
import build_mw16_extended_projective_caches as parent
import encode_mw16_joint_caches as encoder
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=parent.ROOT;CAS=parent.CAS
D=ROOT/'artifacts/local/elliptic-curves/mw16-joint-cache-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve cache proof controller')
    paths=[Path(__file__).resolve(),Path(parent.__file__),Path(encoder.__file__),
           parent.D/'protocol.json',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':8000,'rss_bytes':2147483648,
        'stages':[('table-replay','build_mw16_extended_projective_caches.py','check',3600),
                  ('prepare','encode_mw16_joint_caches.py','prepare',120),
                  ('encode','encode_mw16_joint_caches.py','encode',2400),
                  ('encoding-check','encode_mw16_joint_caches.py','encoding-check',2400),
                  ('score','encode_mw16_joint_caches.py','score',900),
                  ('score-check','encode_mw16_joint_caches.py','score-check',300)],
        'scope':'Wait for the already launched fixed14740-table five-family build. Replay all tables, encode both finite prime bands for every family, replay every byte, then compare all40960 saved target-free short and extended scores and independently replay those comparisons. Preserve failure or censoring without retry. No new parameter selection or point search.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen cache proof sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve one cache proof attempt')
    ledger={'status':'WAITING_FOR_CACHE_BUILD','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            status=cert.read(parent.D/'ledger.json')['status']
            if status=='PASS':break
            if status!='RUNNING' or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait expired')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,script,stage,seconds in p['stages']:
            s=run([sys.executable,str(CAS/script),stage],limits=Limits(seconds,p['rss_bytes']),
                  log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            print('FIVE MW16 CACHE PROOFS',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('cache proof stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','launch']);a=parser.parse_args();globals()[a.stage]()

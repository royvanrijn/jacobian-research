#!/usr/bin/env python3
"""Frozen new-annulus11952 point experiment; all maps precede points."""
import argparse,sys,time
from pathlib import Path
import joint64_r17_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=ROOT/'artifacts/local/elliptic-curves/joint64-controller-v1'
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve joint64 controller')
    names=['finish_joint64_points.py','joint64_r17_pari_batch.py','prepare_joint64_r17_pari_batch.sage','verify_joint64_r17_pari_batch.py','replay_joint64_geometry.py','audit_joint64_clouds_modl.py','score_joint11952_fresh.py','audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py','research_runtime/supervisor.py','research_runtime/store.py']
    paths=[*(CAS/n for n in names),batch.extension.D/'controller/protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_seconds':10800,'rss_bytes':4294967296,'stages':[('freeze','joint64_r17_pari_batch.py',['freeze'],120),('maps','joint64_r17_pari_batch.py',['maps'],1800),('points','joint64_r17_pari_batch.py',['batch'],20000),('verify','verify_joint64_r17_pari_batch.py',[],3600),('odd','audit_joint64_clouds_modl.py',[],1800),('geometry','replay_joint64_geometry.py',[],1200)],'scope':'After full3510-prime slice scoring,32768-survivor replay,4096-scalar-score and64-selection and disjoint validation proofs pass, freeze exactly64 new-annulus11952 point attempts,49 generic17 parity maps each, height125000,10 seconds per chart. All64 map files pass before any point attempt. Two point workers,600 seconds per curve, provisional28 stop; no adaptive wave or retries. Replay all terminal admission histories and retained point clouds modulo2,3,5 and rational geometry. Any catalogue/previous-equation comparison, inventory promotion or portable proof packaging follows separate terminal proofs.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('fixed joint64 sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve point launch')
    ledger={'status':'WAITING_FOR_FRESH_SCORES','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            d=cert.read(batch.extension.D/'controller/ledger.json')
            if d['status']=='PASS':break
            if d['status'] not in ('RUNNING','WAITING_FOR_UPSTREAM_PROOFS') or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait exhausted')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,script,args,seconds in p['stages']:
            s=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('JOINT64 POINT CONTROLLER',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('fixed point/proof stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

#!/usr/bin/env python3
"""Frozen new-annulus11952 point experiment; all maps precede points."""
import argparse,sys,time
from pathlib import Path
import annulus64_r17_pari_batch_v2 as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=ROOT/'artifacts/local/elliptic-curves/annulus64-controller-v2'
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve annulus64 controller')
    names=['finish_annulus64_points_v2.py','annulus64_r17_pari_batch_v2.py','prepare_annulus64_r17_pari_batch_v2.sage','verify_annulus64_r17_pari_batch_v2.py','replay_annulus64_geometry_v2.py','audit_annulus64_clouds_modl_v2.py','score_11952_new_annulus_v2.py','audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py','research_runtime/supervisor.py','research_runtime/store.py']
    paths=[*(CAS/n for n in names),batch.extension.D/'controller/protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_seconds':4800,'rss_bytes':4294967296,'stages':[('freeze','annulus64_r17_pari_batch_v2.py',['freeze'],120),('maps','annulus64_r17_pari_batch_v2.py',['maps'],1800),('points','annulus64_r17_pari_batch_v2.py',['batch'],20000),('verify','verify_annulus64_r17_pari_batch_v2.py',[],3600),('odd','audit_annulus64_clouds_modl_v2.py',[],1800),('geometry','replay_annulus64_geometry_v2.py',[],1200)],'scope':'After32768-cache,4096-score,64-selection and disjoint validation proofs pass, freeze exactly64 new-annulus11952 point attempts,49 generic17 parity maps each, height125000,10 seconds per chart. All64 map files pass before any point attempt. Two point workers,600 seconds per curve, provisional28 stop; no adaptive wave or retries. Replay all terminal admission histories and retained point clouds modulo2,3,5 and rational geometry. Any catalogue/previous-equation comparison, inventory promotion or portable proof packaging follows separate terminal proofs.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('fixed annulus64 sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve point launch')
    ledger={'status':'WAITING_FOR_FRESH_SCORES','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            d=cert.read(batch.extension.D/'controller/ledger.json')
            if d['status']=='PASS':break
            if d['status'] not in ('RUNNING','WAITING_FOR_OUTER_SCAN') or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait exhausted')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,script,args,seconds in p['stages']:
            s=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('ANNULUS64 POINT CONTROLLER',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('fixed point/proof stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

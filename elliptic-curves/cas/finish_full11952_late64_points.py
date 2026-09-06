#!/usr/bin/env python3
"""Frozen selection-gated64 point trial and complete local exact proofs."""
import argparse,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import full11952_late64_r17_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=ROOT/'artifacts/local/elliptic-curves/full11952-late64-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve one late64 controller')
    names=['finish_full11952_late64_points.py','full11952_late_band_selection.py','verify_full11952_late_band_selection.py',
        'full11952_late64_r17_pari_batch.py','prepare_full11952_late64_r17_pari_batch.sage','verify_full11952_late64_r17_pari_batch.py',
        'replay_full11952_late64_geometry.py','audit_full11952_late64_clouds_modl.py','stream_full11952_late64_verification.py']
    paths=[*(CAS/n for n in names),batch.extension.D/'protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':3600,'map_seconds':1800,'point_seconds':20000,'stream_seconds':22000,
        'odd_seconds':1800,'geometry_seconds':1200,'rss_bytes':4294967296,
        'scope':'After the live frozen4096 score process and all fresh validation stages pass, independently replay the fresh equation exclusions and final64 ordering. Freeze64 generic17-only point attempts at125000 and10 seconds per chart, at most3136 boxes. All64 maps must precede any points. Use two point workers and two independent proof workers, stop each curve at provisional28 pending proof. Then check every retained cloud modulo3,5 and all rational geometry. No refill, follow-up wave or automatic retry; any failure/censoring is preserved. Catalogue and model/inventory promotion remain separate exact stages.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen controller sources differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve one late64 execution')
    ledger={'status':'WAITING_FOR_FROZEN_SCORE_TRIAL','rows':[]};checkpoint(out,ledger)
    deadline=time.monotonic()+p['wait_seconds']
    while True:
        upstream=cert.read(batch.extension.D/'controller/ledger.json')
        if upstream['status']=='PASS':break
        if upstream['status']!='RUNNING' or time.monotonic()>deadline:raise ArithmeticError('upstream failed/censored or fixed wait deadline reached')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(out,ledger)
    def stage(name,script,args,seconds,record=True,folder=D):
        path=folder/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve stage '+name)
        result=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=folder/(name+'.log'),checkpoint_path=path,cwd=ROOT)
        ok=result['outcome']=='completed' and result['returncode']==0
        row={'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':result}
        if record:ledger['rows'].append(row);checkpoint(out,ledger)
        print('LATE64 POINT CONTROLLER',name,row['status'],flush=True)
        if not ok:raise ArithmeticError('fixed stage failed/censored; no retry')
        return row
    stage('selection-replay','verify_full11952_late_band_selection.py',[],180,folder=batch.extension.D)
    stage('freeze','full11952_late64_r17_pari_batch.py',['freeze'],120)
    stage('maps','full11952_late64_r17_pari_batch.py',['maps'],p['map_seconds'])
    stage('stream-freeze','stream_full11952_late64_verification.py',['prepare'],120)
    with ThreadPoolExecutor(max_workers=1) as pool:
        point=pool.submit(stage,'points','full11952_late64_r17_pari_batch.py',['batch'],p['point_seconds'],False)
        point_deadline=time.monotonic()+30
        while not (batch.D/'ledger.json').exists():
            if point.done():point.result();raise ArithmeticError('batch exited without ledger')
            if time.monotonic()>point_deadline:raise TimeoutError('point ledger startup gate')
            time.sleep(0.1)
        stage('stream','stream_full11952_late64_verification.py',['launch'],p['stream_seconds'])
        ledger['rows'].append(point.result());checkpoint(out,ledger)
    stage('odd','audit_full11952_late64_clouds_modl.py',[],p['odd_seconds'])
    stage('geometry','replay_full11952_late64_geometry.py',[],p['geometry_seconds'])
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args()
    try:globals()[a.stage]()
    except Exception as exc:
        if a.stage=='launch' and (D/'ledger.json').exists():
            d=cert.read(D/'ledger.json');d.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',d)
        raise

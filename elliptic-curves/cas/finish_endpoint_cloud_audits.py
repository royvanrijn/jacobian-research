#!/usr/bin/env python3
"""Replay completed fixed endpoint trials as they become available, without searches."""
import sys,time,argparse
from pathlib import Path
import endpoint_specialized_parity_trial as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.BATCH/'cloud-verification-v2'

def prepare():
    p=batch.protocol();out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve21-cloud verification protocol')
    files=[Path(__file__).resolve(),CAS/'audit_endpoint_specialized_trial_v2.py',CAS/'audit_recorded_point_mod2_rank_v3.py',CAS/'audit_retained_cloud_modl.py',batch.BATCH/'protocol.json',batch.BATCH/'completion-v1/protocol.json',D/'00.json',D/'index00.supervisor.json']
    checkpoint(out,{'schema':'elliptic-curves.endpoint-cloud-verification.v2','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in files},'rows':p['rows'],'seconds_per_curve':1200,'overall_wall_seconds':19000,'rss_bytes':2147483648,'scope':'Audit each of the same21 endpoints only after its original fixed worker and selected exact history replay pass. Preserve the existing first-curve proof. Each remaining successful curve gets one full retained-cloud build and check modulo2,3,5. No point search is run here. Failed/censored workers or proof stages are recorded and not retried.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen endpoint cloud verification sources differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve cloud driver ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger);deadline=time.monotonic()+p['overall_wall_seconds']
    for index,row in enumerate(p['rows']):
        done=None
        while time.monotonic()<deadline:
            c=cert.read(batch.BATCH/'completion-v1/ledger.json');done=next((r for r in c['rows'] if r['id']==row['id']),None)
            if done:break
            if c['status']!='RUNNING':raise ArithmeticError('terminal completion ledger omits frozen curve')
            time.sleep(5)
        if done is None:raise TimeoutError('declared overall cloud verification deadline')
        if done['status']!='PASS':item={'id':row['id'],'status':'UPSTREAM_FAILED_OR_CENSORED'}
        else:
            path=D/f'index{index:02}.supervisor.json'
            if index==0:s=cert.read(path)
            else:
                if path.exists():raise FileExistsError('preserve single cloud audit attempt')
                s=run([sys.executable,str(CAS/'audit_endpoint_specialized_trial_v2.py'),'--index',str(index)],limits=Limits(p['seconds_per_curve'],p['rss_bytes']),log_path=D/f'index{index:02}.log',checkpoint_path=path,cwd=ROOT)
            item={'id':row['id'],'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':str(path.relative_to(ROOT)),'supervision_sha256':cert.hashed(path)}
        ledger['rows'].append(item);checkpoint(out,ledger);print('ENDPOINT CLOUD',row['id'],item['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES_OR_CENSORING';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

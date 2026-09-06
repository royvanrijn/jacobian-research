#!/usr/bin/env python3
"""Exact history and point proofs as the fixed192 cohort finishes."""
import argparse,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import compact192_r17_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.D/'stream-verification-v1'

def prepare():
    p=batch.protocol();out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve fixed192 streaming verification')
    paths=[Path(__file__).resolve(),CAS/'verify_compact192_r17_pari_batch.py',CAS/'replay_compact192_geometry.py',CAS/'audit_recorded_point_mod2_rank_v3.py',CAS/'audit_retained_cloud_modl.py',batch.D/'protocol.json',batch.D/'maps-ledger.json']
    checkpoint(out,{'schema':'elliptic-curves.compact192-stream-verification.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':p['rows'],'maximum_workers':2,'seconds_per_curve':1500,'overall_wall_seconds':60000,'rss_bytes_per_worker':2147483648,'gate':'All192 exact map files precede the fixed point batch. Verify each completed worker independently as its immutable result becomes available. This changes verification timing only, not selection, points, charts, stopping conditions or search budgets.','early_high_rank_check':'For any returned lower bound at least28, also replay rational map/point provenance and build/check a separate complete-cloud modulo3/5 certificate, with300 seconds per stage. Catalogue novelty and model normalization remain separate before inventory promotion.','scope':'At most192 exact history/mod2 cloud proofs and one early odd-prime check per provisional>=28 curve. Failed or censored source workers and proof failures remain explicit; no point search or automatic retry occurs here.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen streaming verification sources changed')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve streaming proof ledger')
    data={'status':'RUNNING','rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(out,data);deadline=time.monotonic()+p['overall_wall_seconds']
    def one(index):
        row=p['rows'][index];terminal=None
        while time.monotonic()<deadline:
            result=cert.read(batch.D/'ledger.json');candidate=result['rows'][index]
            if candidate['status']!='PENDING':terminal=candidate;break
            if result['status']!='RUNNING':raise ArithmeticError('terminal batch has a pending row')
            time.sleep(5)
        if terminal is None:return {'id':row['id'],'status':'WAIT_DEADLINE_CENSORED'}
        if terminal['status'] not in ('COMPLETE_DECLARED_POINT_ATTEMPT','TARGET_REACHED_PENDING_REPLAY'):return {'id':row['id'],'status':'UPSTREAM_FAILED_OR_CENSORED','upstream_status':terminal['status']}
        folder=D/row['id'];v=batch.D/row['id']/'verification.json'
        s=run([sys.executable,str(CAS/'verify_compact192_r17_pari_batch.py'),'--index',str(index)],limits=Limits(p['seconds_per_curve'],p['rss_bytes_per_worker']),log_path=folder/'verification.log',checkpoint_path=folder/'verification.supervisor.json',cwd=ROOT)
        if s['outcome']!='completed' or s['returncode']!=0:return {'id':row['id'],'status':'FAILED_OR_CENSORED','supervision':s}
        result=cert.read(v);cloud=cert.read(ROOT/result['cloud_certificate'])
        if result['status']!='PASS' or cert.hashed(batch.D/row['id']/'result.json')!=result['input_sha256'] or cert.hashed(ROOT/result['cloud_certificate'])!=result['cloud_sha256']:raise ArithmeticError('streamed immutable point proof differs')
        if cloud['rank_lower_bound']>=28:
            import replay_compact192_geometry as geometry
            raw=cert.read(batch.D/row['id']/'result.json');maps=cert.read(batch.D/row['id']/'maps.json');initial=geometry.tuples(raw['initial_state']['state']['reductions']['points'])
            geometry.geometry(raw,maps,initial,{**batch.protocol(),'maps_path':batch.D/row['id']/'maps.json'});geometry.cloud_check(raw,raw['charts'],ROOT/result['cloud_certificate'],batch.D/row['id']/'result.json')
            odd=batch.ART/('compact192_r17_'+row['id'].replace('-','_')+'_early_modl_v1.json')
            for label,args in [('odd-build',['--input',str(ROOT/result['cloud_certificate']),'--output',str(odd)]),('odd-check',['--check',str(odd)])]:
                q=run([sys.executable,str(CAS/'audit_retained_cloud_modl.py'),*args],limits=Limits(300,p['rss_bytes_per_worker']),log_path=folder/(label+'.log'),checkpoint_path=folder/(label+'.supervisor.json'),cwd=ROOT)
                if q['outcome']!='completed' or q['returncode']!=0:return {'id':row['id'],'status':'EARLY_ODD_CHECK_FAILED_OR_CENSORED','mod2_proof':result,'supervision':q}
            result={**result,'early_geometry':'PASS','early_odd_certificate':str(odd.relative_to(ROOT)),'early_odd_sha256':cert.hashed(odd)}
        return {'id':row['id'],**result,'rank_lower_bound':cloud['rank_lower_bound']}
    with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f];data['rows'][i]=f.result();checkpoint(out,data);print('COMPACT192 EXACT PROOF',p['rows'][i]['id'],data['rows'][i]['status'],data['rows'][i].get('rank_lower_bound'),flush=True)
    data['status']='PASS' if all(r['status']=='PASS' for r in data['rows']) else 'COMPLETE_WITH_FAILURES_OR_CENSORING';checkpoint(out,data)
    if data['status']=='PASS':
        final=batch.D/'verification-ledger.json'
        if final.exists():raise FileExistsError('preserve final192 verification ledger')
        checkpoint(final,data)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

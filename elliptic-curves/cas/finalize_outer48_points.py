#!/usr/bin/env python3
"""Close the fixed outer48 point cohort after its streaming exact proofs finish."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BASE=ROOT/'artifacts/local/elliptic-curves/outer48-r17-pari-v1';D=BASE/'post-batch'

def prepare():
    out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve outer48 finalization protocol')
    names=['finalize_outer48_points.py','certify_outer48_r17_results.py','replay_outer48_geometry.py','audit_outer48_clouds_modl.py','report_outer48_experiment.py']
    paths=[*(CAS/n for n in names),BASE/'protocol.json',BASE/'maps-ledger.json',BASE/'stream-verification-v1/protocol.json']
    checkpoint(out,{'schema':'elliptic-curves.outer48-finalization.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_deadline_seconds':17000,'point_proof_seconds':1200,'geometry_seconds':1200,'odd_cloud_seconds':15000,'aggregate_seconds':600,'rss_bytes':4294967296,'scope':'Wait for all48 point attempts and streaming exact history/mod2 cloud proofs. Only after they pass perform the pinned catalogue comparison and cohort point-proof replay, complete rational geometry/provenance checks, and fixed modulo3/5 audits of all48 retained clouds. Report any odd-prime bound stronger than the mod2 bound explicitly before promotion. No point search, new parameter, adaptive wave, automatic retry or inventory promotion occurs here.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen outer48 finalization inputs differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve outer48 finalization ledger')
    ledger={'status':'WAITING_FOR_FIXED_JOBS','rows':[]};checkpoint(out,ledger);deadline=time.monotonic()+p['wait_deadline_seconds']
    while True:
        a=cert.read(BASE/'ledger.json');b=cert.read(BASE/'stream-verification-v1/ledger.json')
        if a['status'] not in ('RUNNING','COMPLETE_FIXED_BATCH_ATTEMPTS') or b['status'] not in ('RUNNING','PASS'):raise ArithmeticError('upstream non-success retained; finalization stopped')
        if a['status']=='COMPLETE_FIXED_BATCH_ATTEMPTS' and b['status']=='PASS':
            if len(a['rows'])!=48 or len(b['rows'])!=48 or any(r['status']!='PASS' for r in b['rows']):raise ArithmeticError('all48 exact point proofs required')
            break
        if time.monotonic()>deadline:raise TimeoutError('declared upstream wait deadline')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(out,ledger)
    jobs=[('certify','certify_outer48_r17_results.py',[],p['point_proof_seconds']),('proof-replay','certify_outer48_r17_results.py',['--check'],p['point_proof_seconds']),('geometry','replay_outer48_geometry.py',[],p['geometry_seconds']),('odd-clouds','audit_outer48_clouds_modl.py',[],p['odd_cloud_seconds']),('aggregate-build','report_outer48_experiment.py',[],p['aggregate_seconds']),('aggregate-check','report_outer48_experiment.py',['--check'],p['aggregate_seconds'])]
    for name,script,args,seconds in jobs:
        path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve outer48 finalization stage')
        s=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger);print('OUTER48 FINALIZATION',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:raise ArithmeticError('finalization stage failed/censored; no retry')
        if name=='odd-clouds':
            odd=cert.read(BASE/'odd-cloud-audit/ledger.json')
            if odd['status']!='PASS' or len(odd['rows'])!=48:raise ArithmeticError('all48 odd-prime cloud audits required')
            ledger['stronger_odd_prime_bounds']=[{'id':r['id'],'mod2_lower_bound':r['mod2_lower_bound'],'odd_lower_bounds':r['odd_lower_bounds']} for r in odd['rows'] if max(r['odd_lower_bounds'].values())>r['mod2_lower_bound']];checkpoint(out,ledger)
    ledger['status']='PASS_STRONGER_ODD_BOUNDS_REQUIRE_REVIEW' if ledger['stronger_odd_prime_bounds'] else 'PASS';checkpoint(out,ledger)
    print('OUTER48 TERMINAL',ledger['status'],'INVENTORY/NOVELTY PROMOTION STILL SEPARATE',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

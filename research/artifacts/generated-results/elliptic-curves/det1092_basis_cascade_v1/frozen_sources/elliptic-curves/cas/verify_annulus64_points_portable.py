#!/usr/bin/env python3
"""194 isolated exact replay stages for the completed full11952 point cohort."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/annulus64-point-portable-v1'
OUT=ART/'annulus64_point_portable_replay_v1.json'

def main():
    path=ART/'annulus64_point_evidence_v1.json';manifest=cert.read(path);workspace=D/'workspace';folder=D/'verification'
    if workspace.exists() or OUT.exists():raise FileExistsError('preserve isolated full64 point replay')
    archive=ROOT/manifest['archive']
    if manifest['required_base_archives'] or cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('standalone point archive differs')
    workspace.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('extracted member differs')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves'
    p=cert.read(workspace/'artifacts/local/elliptic-curves/annulus64-r17-pari-v2/protocol.json');jobs=[]
    if len(p['rows'])!=64:raise ArithmeticError('fixed64 roster required')
    for i,row in enumerate(p['rows']):
        jobs.append((row['id']+'-history',[str(cas/'annulus64_r17_pari_batch_v2.py'),'replay','--index',str(i)],600))
        for suffix,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
            proof=art/('annulus64_r17_'+row['id'].replace('-','_')+'_'+suffix+'_v1.json')
            jobs.append((row['id']+'-'+suffix,[str(cas/script),'--check',str(proof)],300))
    jobs.extend([('geometry',[str(cas/'replay_annulus64_geometry_v2.py')],1200),('points-and-equations',[str(cas/'certify_annulus64_r17_results.py'),'--check'],1200)])
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'seconds':s} for n,a,s in jobs],'rss_bytes':2147483648,'maximum_workers':1,'scope':manifest['claim_boundary']})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        s=run([sys.executable,*args],limits=Limits(seconds,2147483648),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'ledger.json',ledger);print('ISOLATED ANNULUS64',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(folder/'ledger.json',ledger);raise ArithmeticError('isolated proof failure; no retry')
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('immutable member changed during replay')
    ledger['status']='PASS';checkpoint(folder/'ledger.json',ledger)
    paths=[Path(__file__).resolve(),path,folder/'protocol.json',folder/'ledger.json']
    checkpoint(OUT,{'schema':'elliptic-curves.annulus64-point-portable-replay.v1','status':'PASS','logical_stages':len(jobs),'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'archive_sha256':manifest['archive_sha256'],'ledger':ledger,'claim_boundary':manifest['claim_boundary']})
    print('ISOLATED ANNULUS64',len(jobs),'STAGES PASS',flush=True)

if __name__=='__main__':main()

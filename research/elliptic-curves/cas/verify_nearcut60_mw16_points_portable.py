#!/usr/bin/env python3
"""Isolated exact proofs for certified rows, preserving all unresolved allocations."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/nearcut60-mw16-point-portable-v1'
OUT=ART/'nearcut60_mw16_point_portable_replay_v1.json'


def main():
    path=ART/'nearcut60_mw16_point_evidence_v1.json';manifest=cert.read(path);workspace=D/'workspace';folder=D/'verification'
    if workspace.exists() or OUT.exists():raise FileExistsError('preserve isolated comparison replay')
    archive=ROOT/manifest['archive']
    if manifest['required_base_archives'] or cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('standalone archive differs')
    workspace.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        names=z.namelist()
        if len(names)!=len(set(names)) or any(Path(n).is_absolute() or '..' in Path(n).parts for n in names):raise ArithmeticError('unsafe or duplicate member')
        if set(names)!={r['path'] for r in manifest['files']}:raise ArithmeticError('complete member roster differs')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('extracted member differs')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves'
    report=cert.read(art/'nearcut60_mw16_experiment_v1.json')
    p=cert.read(workspace/'artifacts/local/elliptic-curves/nearcut60-mw16-pari-v1/protocol.json')
    if len(p['rows'])!=60 or len(manifest['cases'])!=60 or [r['id'] for r in manifest['cases']]!=[r['id'] for r in report['rows']] or [r['id'] for r in p['rows']]!=[r['id'] for r in report['rows']]:raise ArithmeticError('all60 allocations required')
    jobs=[];ranks=[];unresolved=[]
    for i,(case,row) in enumerate(zip(manifest['cases'],report['rows'])):
        if case['index']!=i or any(case[k]!=row[k] for k in ('id','arm','family','parameter','rank_lower_bound','certified_gain','worker_status','verification_status','completed_boxes','attempted_boxes')):raise ArithmeticError('manifest outcome differs')
        if case['replay_point_proof']:
            if row['verification_status']!='PASS' or row['rank_lower_bound'] is None:raise ArithmeticError('unverified row scheduled as point proof')
            ident=case['id'];cloud=workspace/case['cloud_path'];odd=workspace/case['odd_path']
            verified=cert.read(workspace/'artifacts/local/elliptic-curves/nearcut60-mw16-pari-v1'/ident/'verification.json')
            if verified['status']!='PASS' or verified['cloud_path']!=case['cloud_path'] or verified['odd_path']!=case['odd_path'] or cert.hashed(cloud)!=verified['cloud_sha256'] or cert.hashed(odd)!=verified['odd_sha256']:raise ArithmeticError('bound standalone certificates differ')
            lower=max(cert.read(cloud)['rank_lower_bound'],*(a['finite_column_rank'] for a in cert.read(odd)['audits']))
            if row['rank_lower_bound']!=lower or row['certified_gain']!=lower-16:raise ArithmeticError('standalone rank binding differs')
            jobs.extend([
                (ident+'-history',[str(cas/'nearcut60_mw16_pari_batch.py'),'replay','--index',str(i)],600),
                (ident+'-geometry',[str(cas/'verify_nearcut60_mw16_points.py'),'geometry-check','--index',str(i)],300),
                (ident+'-mod2',[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(cloud)],300),
                (ident+'-modl',[str(cas/'audit_retained_cloud_modl.py'),'--check',str(odd)],300)])
            ranks.append({'id':ident,'rank_lower_bound':lower,'certified_gain':lower-16})
        else:
            if row['verification_status']=='PASS' or row['rank_lower_bound'] is not None or row['certified_gain'] is not None:raise ArithmeticError('unresolved allocation has certified rank')
            unresolved.append(case['id'])
    if len(ranks)!=manifest['certified_rows'] or len(unresolved)!=manifest['unresolved_rows']:raise ArithmeticError('certified/unresolved coverage differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),
        'jobs':[{'name':name,'args':args,'seconds':seconds} for name,args,seconds in jobs],
        'rss_bytes':2147483648,'maximum_workers':1,'allocated_rows':60,'certified_rows':len(ranks),
        'unresolved_ids':unresolved,'scope':manifest['claim_boundary']})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        s=run([sys.executable,*args],limits=Limits(seconds,2147483648),cwd=workspace,log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'))
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'ledger.json',ledger)
        print('ISOLATED NEARCUT60',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(folder/'ledger.json',ledger);raise ArithmeticError('isolated point replay failed; no retry')
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('immutable member changed during replay')
    ledger['status']='PASS';checkpoint(folder/'ledger.json',ledger)
    sources=[Path(__file__).resolve(),path,folder/'protocol.json',folder/'ledger.json']
    checkpoint(OUT,{'schema':'elliptic-curves.nearcut60-point-portable-replay.v1','status':'PASS',
        'logical_stages':len(jobs),'allocated_curves':60,'certified_rows':ranks,'unresolved_ids':unresolved,
        'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in sources},
        'archive_sha256':manifest['archive_sha256'],'ledger':ledger,'claim_boundary':manifest['claim_boundary']})
    print('ISOLATED NEARCUT60',len(jobs),'PROOF STAGES PASS;',len(unresolved),'UNRESOLVED ALLOCATIONS',flush=True)


if __name__=='__main__':main()

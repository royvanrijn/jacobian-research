#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/endpoint-point-trial-portable-v1'

def main():
    manifest=cert.read(ART/'endpoint_point_trial_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
    if workspace.exists():raise FileExistsError('preserve isolated paired replay')
    workspace.mkdir(parents=True)
    for m in [*manifest['required_base_archives'],manifest]:
        archive=ROOT/m['archive']
        if cert.hashed(archive)!=m['archive_sha256']:raise ArithmeticError('archive changed')
        with zipfile.ZipFile(archive) as z:
            if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
            z.extractall(workspace)
    for r in [*manifest['files'],*manifest['inherited_exact_members']]:
        if cert.hashed(workspace/r['path'])!=r['sha256']:raise ArithmeticError('isolated member hash differs')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';local=workspace/'artifacts/local/elliptic-curves'
    jobs=[('endpoint-geometry',[str(cas/'replay_endpoint_specialized_geometry.sage')],300),('endpoint-summary',[str(cas/'report_endpoint_specialized_trial.py'),'--check'],180),('six-small-prime-saturation',[str(cas/'certify_new27_small_prime_saturation.py'),'--check'],180)]
    for i,row in enumerate(cert.read(local/'endpoint-specialized-parity-v1/protocol.json')['rows']):
        jobs.append((row['id']+'-history',[str(cas/'endpoint_specialized_parity_trial.py'),'replay','--index',str(i)],600))
        for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
            jobs.append((row['id']+'-'+label,[str(cas/script),'--check',str(art/('endpoint_specialized_'+row['id'].replace('-','_')+'_'+label+'_v1.json'))],300))
    if len(jobs)!=66:raise ArithmeticError('fixed66-stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'endpoint_point_trial_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated replay of21 fixed endpoint admission histories, complete retained point clouds modulo2,3,5, exact sampled parity/norm and rational map/provenance identities, endpoint aggregate, and six own27-point subgroup2,3,5-saturation proofs. No point searches or whole-curve rank upper bounds.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        executable='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if args[0].endswith('.sage') else sys.executable
        r=run([executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE ENDPOINT POINT TRIAL',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

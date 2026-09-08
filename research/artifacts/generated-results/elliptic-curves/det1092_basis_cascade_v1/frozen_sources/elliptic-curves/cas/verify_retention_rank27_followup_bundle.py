#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/retention-rank27-followup-portable-v1'

def main():
    manifest=cert.read(ART/'retention_rank27_followup_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';p=cert.read(workspace/'artifacts/local/elliptic-curves/retention24-r17-pari-v1/protocol.json')
    jobs=[('geometry602',[str(cas/'replay_retention_rank27_followup_geometry.py')],180)]
    for identifier in ('103b2_733','11952_113'):
        for modulus,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
            jobs.append((identifier+'-'+modulus,[str(cas/script),'--check',str(art/('retention_rank27_'+identifier+'_all_retained_'+modulus+'_v1.json'))],180))
    if len(jobs)!=5:raise ArithmeticError('fixed5 followup portable roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'retention_rank27_followup_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'All602 adaptive chart centres, parity rosters, exact maps, raw points and full-cloud provenance, plus both1426/1928-point clouds modulo2,3,5. No point search. All admission histories passed separately and are retained.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE RETENTION RANK27',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

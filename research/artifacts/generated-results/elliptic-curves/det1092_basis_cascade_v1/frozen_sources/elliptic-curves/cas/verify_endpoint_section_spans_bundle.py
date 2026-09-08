#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/endpoint-section-spans-portable-v1'

def main():
    manifest=cert.read(ART/'endpoint_section_spans_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('exact-section-relations',[str(cas/'certify_endpoint_section_relations.py'),'check'],120),('known-subgroup-indices',[str(cas/'verify_endpoint_section_lattice_indices.py')],60),('span-summary',[str(cas/'report_endpoint_section_spans.py'),'--check'],120)]
    if len(jobs)!=3:raise ArithmeticError('fixed three-stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'endpoint_section_spans_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact section relations and finite independent subsets on21 endpoints, independent integer matrix checks of all21 known-section indices, and the aggregate. Numerical proposals are not trusted as proofs. No point search or whole-curve upper bound.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE ENDPOINT SECTION SPANS',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/skew-endpoint-portable-v1'

def main():
    manifest=cert.read(ART/'skew_endpoint_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('coefficient-bounds',[str(cas/'audit_r17_parameter_box_skew.py'),'--check'],60),('skew-population',[str(cas/'scan_skew_r17_boxes.py'),'replay'],180),('extended2048',[str(cas/'extend_skew_r17_scores.py'),'replay'],180),('trace-benchmark',[str(cas/'extend_skew_r17_scores.py'),'benchmark-check'],60),('point-proofs',[str(cas/'certify_skew8_r17_results.py'),'--check'],180),('point-geometry',[str(cas/'replay_skew8_geometry.py')],180),('skew-aggregate',[str(cas/'report_skew_r17_experiment.py'),'--check'],180),('endpoint-audit',[str(cas/'audit_compact_atlas_endpoints_v2.py'),'check'],180),('endpoint-summary',[str(cas/'report_compact_endpoint_audit.py'),'--check'],60)]
    for i,row in enumerate(cert.read(local/'skew8-r17-pari-v1/protocol.json')['rows']):
        jobs.append((row['id']+'-history',[str(cas/'skew8_r17_pari_batch.py'),'replay','--index',str(i)],300))
        for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
            jobs.append((row['id']+'-'+label,[str(cas/script),'--check',str(art/('skew8_r17_'+row['id'].replace('-','_')+'_'+label+'_v1.json'))],180))
    for row in cert.read(local/'compact-endpoint-odd-primes-v1/ledger.json')['rows']:
        jobs.append((row['family']+'-'+row['endpoint']+'-modl',[str(cas/'audit_endpoint_section_cloud_modl.py'),'--check',str(workspace/row['output'])],180))
    if len(jobs)!=54:raise ArithmeticError('fixed54 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'skew_endpoint_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact coefficient bounds, four skew scans and2048 extended scores, eight admission histories and full clouds modulo2,3,5,368 rational maps and point provenance,22 zero/infinity evaluations and21 full generic-section odd-modulus checks. The failed endpoint-format assumption remains archived. No new point search, optimality theorem or rank upper bound.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE SKEW ENDPOINTS',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

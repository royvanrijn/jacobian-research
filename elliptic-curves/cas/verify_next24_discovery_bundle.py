#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/next24-portable-v1'

def main():
    manifest=cert.read(ART/'next24_discovery_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';p=cert.read(workspace/'artifacts/local/elliptic-curves/next24-r17-extended-pari-v1/protocol.json')
    jobs=[('geometry-and-provenance',[str(cas/'replay_next24_point_geometry.py')],180),('minimal-models',[str(cas/'certify_next24_high_rank_minimal_v4.py'),'--check'],120),('batch-certificates',[str(cas/'certify_next24_r17_results.py'),'--check'],120),('inventory62',[str(cas/'replay_inventory_v5_memory.py'),'--output',str(folder/'inventory62-result.json')],120),('published-visibility',[str(cas/'audit_native11952_published_visibility_v2.py'),'--check'],120)]
    for row in p['rows']:
        name='next24_r17_extended_'+row['id'].replace('-','_')+'_mod2_v1.json';jobs.append(('cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/name)],120))
    jobs.extend([('adaptive-mod2',[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/'next24_rank27_all_retained_mod2_v1.json')],120),('adaptive-mod3-mod5',[str(cas/'audit_retained_cloud_modl.py'),'--check',str(art/'next24_rank27_all_retained_modl_v1.json')],180)])
    jobs.append(('local-minimality-regression',['-m','unittest','discover','-s',str(workspace/'elliptic-curves/tests'),'-p','test_next24_local_minimality.py'],60))
    if len(jobs)!=32:raise ArithmeticError('fixed portable roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'next24_discovery_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact centres, maps, raw points and complete-cloud provenance for1080 initial and301 adaptive charts; all24 standalone cloud certificates, three minimal model proofs,62-curve inventory andCSV, adaptive modulo2/3/5 proofs and separate retrospective visibility. All admission/archive histories passed local checks and are retained; this isolated run does not repeat them or run a point search.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE NEXT24',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

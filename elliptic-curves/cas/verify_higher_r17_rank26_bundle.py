#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/higher-r17-rank26-portable-v1'

def main():
    manifest=cert.read(ART/'higher_r17_rank26_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('scanner12',[str(cas/'replay_higher_r17_stratified_portable.py')],120),('trace-benchmark',[str(cas/'extend_higher_r17_stratified.py'),'benchmark-check'],120),('traces6144',[str(cas/'extend_higher_r17_stratified.py'),'replay'],180),('geometry1080',[str(cas/'replay_higher24_geometry.py')],120),('cohort24',[str(cas/'certify_higher24_r17_results.py'),'--check'],120),('initial-minimal25',[str(cas/'certify_higher24_rank25_minimal.py'),'--check'],120),('geometry301',[str(cas/'replay_higher_rank25_followup_geometry.py')],120),('adaptive26',[str(cas/'certify_higher25_followup_results.py'),'--check'],120),('minimal26',[str(cas/'certify_higher26_minimal.py'),'--check'],120),('inventory100',[str(cas/'replay_inventory_v11_memory.py'),'--output',str(folder/'inventory-replay.json')],180),('sage-export',[str(cas/'export_higher26_sage.py'),'--check'],30),('incidence24',[str(cas/'replay_higher24_incidence.py'),'--input',str(art/'higher24_cross_family_j_incidence_v1.json'),'--output',str(folder/'incidence-replay.json')],120),('incidence1200',[str(cas/'certify_inventory100_incidence_v2.py'),'--check',str(art/'inventory100_cross_family_incidence_v2.json')],120),('visibility-cost',[str(cas/'audit_higher24_visibility_cost.py'),'--check'],120)]
    p=cert.read(local/'higher24-r17-pari-v1/protocol.json')
    for row in p['rows']:jobs.append(('cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/('higher24_r17_'+row['id'].replace('-','_')+'_mod2_v1.json'))],120))
    for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:jobs.append(('adaptive-cloud-'+label,[str(cas/script),'--check',str(art/('higher_rank25_11952_069_all_retained_'+label+'_v1.json'))],180))
    if len(jobs)!=40:raise ArithmeticError('fixed40 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'higher_r17_rank26_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated saved-score, exact geometry and point-cloud provenance checks for1080 initial boxes and301 adaptive boxes;24 initial point proofs, one improved26-point minimal equation,100 IDs/CSV and1200 incidence binding, and descriptive coefficient/visibility comparisons. No new scanner, point search or prime trace. Completed admission histories remain available.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE HIGHER R17 RANK26',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

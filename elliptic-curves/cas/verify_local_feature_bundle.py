#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/local-feature-portable-v1'

def main():
    manifest=cert.read(ART/'local_feature_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('scaled-selector',[str(cas/'select_higher13_scaled_candidates.py'),'replay'],120),('split-node-selector',[str(cas/'score_higher_split_nodes.py'),'replay'],180),('scaled-proof',[str(cas/'certify_scaled13_24_r17_results.py'),'--check'],120),('scaled-geometry',[str(cas/'replay_scaled13_24_geometry.py')],120),('split-node-proof',[str(cas/'certify_splitnode24_r17_results.py'),'--check'],120),('split-node-geometry',[str(cas/'replay_splitnode24_geometry.py')],120),('followup-proof',[str(cas/'certify_scaled13_25_followup_results.py'),'--check'],120),('followup-geometry',[str(cas/'replay_scaled13_rank25_followup_geometry.py')],120),('minimal25',[str(cas/'certify_scaled13_24_rank25_minimal.py'),'--check'],120),('sage-export',[str(cas/'export_scaled13_rank25_sage.py'),'--check'],60),('sage-load',[str(art/'new_scaled13_rank25_curve_11952.sage')],60),('inventory101',[str(cas/'replay_inventory_v12_memory.py'),'--output',str(folder/'inventory101.json')],300),('incidence12',[str(cas/'replay_inventory101_added_incidence.py'),'--input',str(art/'inventory101_added_cross_family_j_incidence_v1.json'),'--output',str(folder/'incidence12.json')],120),('incidence1212',[str(cas/'certify_inventory101_incidence.py'),'--check',str(art/'inventory101_incidence_v1.json')],60),('aggregate',[str(cas/'report_local_feature_experiments.py'),'--check'],120)]
    for cohort in cert.read(art/'local_feature_experiments_v1.json')['cohorts']:
        for row in cohort['rows']:
            name=row['id'].replace('-','_')
            for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
                jobs.append((cohort['cohort']+'-'+row['id']+'-'+label,[str(cas/script),'--check',str(art/(cohort['cohort']+'_'+name+'_'+label+'_v1.json'))],180))
    for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
        jobs.append(('followup-'+label,[str(cas/script),'--check',str(art/('scaled13_rank25_11952_300_all_retained_'+label+'_v1.json'))],300))
    if len(jobs)!=113:raise ArithmeticError('fixed113 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'local_feature_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated two-selector and48-curve replay:2160 initial chart maps and raw-point provenance, full clouds modulo2,3,5,301 adaptive maps and3260-point union, new25 minimality and Sage export,101-curve inventory,12 new incidence pairs and1212-pair aggregate binding. Original admission histories remain retained. No new point search or rank upper bound.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='sage-load' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE LOCAL FEATURES',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

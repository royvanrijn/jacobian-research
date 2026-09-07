#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/product22-comparison-portable-v1'

def main():
    manifest=cert.read(ART/'product22_comparison_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('saved6144-scores',[str(cas/'compare_higher_r17_product_score.py'),'replay'],180),('curve-proofs',[str(cas/'certify_product22_r17_results.py'),'--check'],120),('geometry988',[str(cas/'replay_product22_geometry.py')],120)]
    for row in cert.read(art/'product22_comparison_v1.json')['rows']:
        name=row['id'].replace('-','_')
        for label,script,prefix in [('mod2','audit_recorded_point_mod2_rank_v3.py','product24'),('modl','audit_retained_cloud_modl.py','product22')]:
            jobs.append((row['id']+'-'+label,[str(cas/script),'--check',str(art/(prefix+'_r17_'+name+'_'+label+'_v1.json'))],180))
    jobs.extend([('visibility48',[str(cas/'audit_generic_point_box_visibility.py'),'replay'],180),('visibility-known28',[str(cas/'audit_native28_generic_visibility.py'),'--check'],120),('aggregate',[str(cas/'report_product22_comparison.py'),'--check'],60)])
    if len(jobs)!=50:raise ArithmeticError('fixed50 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'product22_comparison_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated saved6144 score replay,22 exact curve proofs,988 map and raw-point geometry checks,all22 complete clouds modulo2,3,5,and both exact generic-representative visibility audits. The completed admission histories remain retained. No new point search or rank upper bound.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE PAIRED22',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/paired-rank27-portable-v1'

def main():
    manifest=cert.read(ART/'paired_rank27_discovery_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';p=cert.read(workspace/'artifacts/local/elliptic-curves/fresh-r17-paired-pari-v1/protocol.json')
    jobs=[('trace-extension',[str(cas/'extend_retained_r17_prime_scores.py'),'replay'],180),('geometry-and-provenance',[str(cas/'replay_paired_point_geometry.py')],300),('minimal-models',[str(cas/'certify_paired_high_rank_minimal_v2.py'),'--check'],180),('first-rank27',[str(cas/'certify_paired_rank27.py'),'--check'],120),('batch-certificates',[str(cas/'certify_fresh_r17_paired_results.py'),'--check'],180),('inventory47',[str(cas/'export_new_high_rank_curve_index_v4.py'),'--check',str(art/'new_high_rank_curve_index_v4.json')],180)]
    for row in p['rows']:
        name='fresh_r17_paired_'+row['id'].replace('-','_')+'_mod2_v1.json';jobs.append(('cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v2.py'),'--check',str(art/name)],120))
    for prefix in ('paired_rank27','paired_second27'):
        jobs.extend([(prefix+'-mod2',[str(cas/'audit_recorded_point_mod2_rank_v2.py'),'--check',str(art/(prefix+'_all_retained_mod2_v1.json'))],180),(prefix+'-mod3-mod5',[str(cas/'audit_retained_cloud_modl.py'),'--check',str(art/(prefix+'_all_retained_modl_v1.json'))],180)])
    for name in ('memory_rank_certificate','r17_prime_extension'):jobs.append(('test-'+name,['-m','unittest','discover','-s',str(workspace/'elliptic-curves/tests'),'-p','test_'+name+'.py'],60))
    if len(jobs)!=35:raise ArithmeticError('fixed portable roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'paired_rank27_discovery_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated pure-Python proofs:768 frozen trace rosters/scores,1639 exact rational centres/geometries/raw point lists with complete-cloud provenance, all standalone rank certificates and47-curve inventory, both adaptive clouds modulo2/3/5, targeted tests. Initial1037 and adaptive602 admission/archive histories passed separate retained local replays; this isolated run does not repeat those histories or invoke either point-search engine.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE PAIRED RANK27',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='COMPLETE_DECLARED_REPLAY_ATTEMPTS';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

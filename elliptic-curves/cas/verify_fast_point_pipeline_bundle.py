#!/usr/bin/env python3
"""Isolated witness, continuation, rank and unit checks for the fast pipeline."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/fast-point-pipeline-portable-v1'

def main():
    manifest=cert.read(ART/'fast_point_pipeline_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
    if workspace.exists():raise FileExistsError('preserve isolated replay')
    workspace.mkdir(parents=True)
    for m in [*manifest['required_base_archives'],manifest]:
        archive=ROOT/m['archive']
        if cert.hashed(archive)!=m['archive_sha256']:raise ArithmeticError('archive changed')
        with zipfile.ZipFile(archive) as z:
            if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
            z.extractall(workspace)
    for r in [*manifest['files'],*manifest['inherited_exact_members']]:
        if cert.hashed(workspace/r['path'])!=r['sha256']:raise ArithmeticError('isolated member hash differs')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves'
    jobs=[('all602-geometries',[str(cas/'replay_all_fast_backend_geometry.py')],300),('strict-calibrations',[str(cas/'replay_backend_calibrations_strict.py'),'--check'],180)]
    for case in ('rank26','small-conductor'):jobs.append(('continuation-'+case,[str(cas/'continue_fixed_pari_search.py'),'replay','--case',case],900))
    for prefix in ('rank26','small_conductor'):
        jobs.append((prefix+'-mod2',[str(cas/'audit_recorded_point_mod2_rank_v2.py'),'--check',str(art/(prefix+'_all_current_retained_mod2_v1.json'))],180))
        jobs.append((prefix+'-mod3-mod5',[str(cas/'audit_retained_cloud_modl.py'),'--check',str(art/(prefix+'_all_retained_modl_v1.json'))],180))
    for name in ('pari_pointed_backend','preloaded_prime_state','rotated_observation_state','retained_cloud_modl','strict_calibration_roster'):
        jobs.append(('test-'+name,['-m','unittest','discover','-s',str(workspace/'elliptic-curves/tests'),'-p','test_'+name+'.py'],60))
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'fast_point_pipeline_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated checks without point search or external CAS: all602 new adaptive geometries/raw point sets, all278 continuation admission/archive histories, complete point-cloud rank proofs, strict calibration rosters and targeted tests. The285+39 original slow admission histories were independently replayed locally and are retained; they are not redundantly re-executed here.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE FAST PIPELINE',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='COMPLETE_DECLARED_REPLAY_ATTEMPTS';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

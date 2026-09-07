#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/height-and-discarded-portable-v1'

def main():
    manifest=cert.read(ART/'height_and_discarded_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';p=cert.read(workspace/'artifacts/local/elliptic-curves/discarded12-r17-pari-v1/protocol.json')
    jobs=[('geometry828',[str(cas/'replay_height_and_discarded_geometry.py')],180),('oracle-translations196',[str(cas/'replay_native11952_translated_visibility_v2.py'),'--output',str(folder/'oracle-replay.json')],120),('height-control',[str(cas/'certify_native29_height_control_v2.py'),'--check'],120),('saved-trace-rosters768',[str(cas/'extend_r17_discarded_shard_scores.py'),'replay'],120),('batch-certificates12',[str(cas/'certify_discarded12_r17_results.py'),'--check'],120),('new-minimal-rank26',[str(cas/'certify_discarded_rank26_minimal.py'),'--check'],120),('inventory70',[str(cas/'replay_inventory_v6_memory.py'),'--output',str(folder/'inventory70-result.json')],120)]
    for row in p['rows']:
        jobs.append(('cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/('discarded12_r17_'+row['id'].replace('-','_')+'_mod2_v1.json'))],120))
    for height in (100000,125000,1000000):jobs.append(('control-cloud-'+str(height),[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/f'native11952_height{height}_mod2_v1.json')],120))
    for identifier in ('paired_074d9','paired_11952','next24_11952'):
        for modulus,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:jobs.append(('followup-'+identifier+'-'+modulus,[str(cas/script),'--check',str(art/f'new27_height125_{identifier}_{modulus}_v1.json')],180))
    if len(jobs)!=28:raise ArithmeticError('fixed28 portable roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'height_and_discarded_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact geometry and raw-point/cloud provenance for828 charts,196 oracle group words, three known-control point proofs,768 saved trace rosters and fixed12 selection, all twelve new-curve proofs, one minimal rank26 model,70-curve inventory andCSV, and three prospective follow-up clouds modulo2,3,5. No point search. All admission/archive histories passed local replay and are retained but not repeated here.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE HEIGHT/DISCARDED',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

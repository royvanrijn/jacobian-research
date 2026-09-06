#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/retention24-discovery-portable-v1'

def main():
    manifest=cert.read(ART/'retention24_discovery_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('geometry1080',[str(cas/'replay_retention24_geometry.py')],180),('retained-short-scores6144',[str(cas/'retain512_compact_r17.py'),'replay'],180),('extended-trace-rosters4608',[str(cas/'extend_retention512_r17_scores.py'),'replay'],180),('oracle-translations98',[str(cas/'replay_native11952_rank28_coset_visibility.py'),'--output',str(folder/'oracle-replay.json')],120),('batch-certificates24',[str(cas/'certify_retention24_r17_results.py'),'--check'],120),('high-rank-minimal',[str(cas/'certify_retention_high_rank_minimal.py'),'--check'],120),('inventory-v7',[str(cas/'replay_inventory_v7_memory.py'),'--output',str(folder/'inventory-v7-result.json')],180)]
    for row in p['rows']:
        jobs.append(('cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/('retention24_r17_'+row['id'].replace('-','_')+'_mod2_v1.json'))],120))
    for name in ('compact','latest7','latest8','latest23','retention'):
        jobs.append(('incidence-'+name,[str(cas/('replay_'+name+'_cross_family_incidence.py')),'--input',str(art/(name+'_cross_family_j_incidence_v1.json')),'--output',str(folder/(name+'-incidence-replay.json'))],120))
    jobs += [('incidence-aggregate89',[str(cas/'certify_inventory89_incidence.py'),'--check',str(art/'inventory89_cross_family_incidence_v1.json')],120),('catalogue593',[str(cas/'refresh_retention_catalogue.py'),'--check'],120)]
    # The generic transport is a separate Sage subprocess; the outer command below selects its interpreter.
    jobs.append(('generic-transport',[str(cas/'audit_compact_published_r17_transport_v3.sage'),'--check',str(art/'compact_published_r17_generic_transport_v1.json')],120))
    if len(jobs)!=39:raise ArithmeticError('fixed39 portable roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'retention24_discovery_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact geometry and raw-point/cloud provenance for1080 charts,6144 short scores and all old128 prefixes,4608 extended trace rosters and fixed24 selection,98 retrospective oracle translations,24 exact point proofs, all new>=27 minimal models, the89-curve V7 inventory/CSV, all1068 incidence checks and generic transport, and the593-equation refreshed comparison. No point search. All admission/archive histories passed local replay and are retained but not repeated here.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE RETENTION24',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/million-height-mw16-portable-v1'

def main():
    manifest=cert.read(ART/'million_height_mw16_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('diagnostic-'+stage,[str(cas/'replay_million_and_mw16_diagnostics.py'),stage],180) for stage in ('native','sieve','traces','neighbours')]
    jobs += [('scores1280',[str(cas/'extend_mw16_retained_prime_scores.py'),'replay'],180),('geometry860',[str(cas/'replay_extended20_mw16_geometry.py')],120),('cohort20',[str(cas/'certify_extended20_mw16_results.py'),'--check'],120),('unions20',[str(cas/'combine_extended20_mw16_retained_points.py'),'--check'],180),('minimal26s',[str(cas/'certify_extended20_mw16_minimal.py'),'--check'],120),('inventory98',[str(cas/'replay_inventory_v8_memory.py'),'--output',str(folder/'inventory-replay.json')],180),('sage-exports',[str(cas/'export_extended20_mw16_sage.py'),'--check'],30),('incidence108',[str(cas/'replay_extended20_mw16_incidence_v2.py'),'--input',str(art/'extended20_mw16_cross_family_j_incidence_v1.json'),'--output',str(folder/'incidence-replay.json')],120),('incidence1176',[str(cas/'certify_inventory98_incidence.py'),'--check',str(art/'inventory98_cross_family_incidence_v1.json')],120)]
    for f in range(1,6):
        family=f'a1-fibration-{f:02}';jobs.append(('parity-'+family,[str(cas/'replay_mw16_exact_maximum_parities.py'),'--input',str(local/'mw16-exact-maximum-parities-v2'/family/'result.json'),'--output',str(folder/(family+'-parity.json'))],180))
    for kind,index,stem in [('mw16',0,'mw16_top25_mw16_04'),('mw16',1,'mw16_top25_mw16_05'),('million',0,'new27_million_retention_11952')]:
        jobs.append((stem+'-geometry',[str(cas/'audit_mw16_and_million_followups.py'),'--kind',kind,'--index',str(index),'--check'],120))
        for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:jobs.append((stem+'-'+label,[str(cas/script),'--check',str(art/(stem+'_all_retained_'+label+'_v1.json'))],180))
    p=cert.read(local/'extended20-mw16-pari-v1/protocol.json')
    for row in p['rows']:jobs.append(('union-cloud-'+row['id'],[str(cas/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(art/('extended20_mw16_'+row['id'].replace('-','_')+'_union_mod2_v1.json'))],120))
    if len(jobs)!=47:raise ArithmeticError('fixed47 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'million_height_mw16_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact diagnostic, score, geometry, point proof/union, minimality,98-curve inventory/CSV, all327680 parity upper witnesses and42 exact maximum minima,108 new incidence pairs and aggregate1176, plus old25 and million-height point clouds modulo2,3,5. No new point search or prime trace call. Admission histories passed separately and remain available.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE MILLION MW16',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

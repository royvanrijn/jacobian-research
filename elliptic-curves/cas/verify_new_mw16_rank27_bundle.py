#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/new-mw16-rank27-portable-v1'

def main():
    manifest=cert.read(ART/'new_mw16_rank27_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[('geometry903',[str(cas/'replay_new_mw16_rank26_followup_geometry.py')],120),('geometry301',[str(cas/'replay_new_mw16_rank27_direction_geometry.py')],120),('followup-three',[str(cas/'certify_new_mw16_followup_results.py'),'--check'],120),('minimal27',[str(cas/'certify_new_mw16_rank27_minimal.py'),'--check'],120),('inventory98',[str(cas/'replay_inventory_v9_memory.py'),'--output',str(folder/'inventory-replay.json')],180),('sage-export',[str(cas/'export_new_mw16_rank27_sage.py'),'--check'],30),('incidence1176',[str(cas/'certify_inventory98_incidence_v2.py'),'--check',str(art/'inventory98_cross_family_incidence_v2.json')],120),('public-height',[str(cas/'audit_public_compact_parameter_heights.py'),'--check'],60),('periodic-benchmark',[str(cas/'benchmark_periodic_nagao_scanner.py'),'replay'],120),('periodic-strict',[str(cas/'replay_periodic_nagao_scanner_portable.py')],120)]
    for stem in ['mw16_new26_a1_fibration_01_015','mw16_new26_a1_fibration_01_052','mw16_new26_a1_fibration_02_014','mw16_new27_direction_a1_fibration_01_052']:
        for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:jobs.append((stem+'-'+label,[str(cas/script),'--check',str(art/(stem+'_all_retained_'+label+'_v1.json'))],180))
    if len(jobs)!=18:raise ArithmeticError('fixed18 portable stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'new_mw16_rank27_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated exact geometry and point-cloud provenance for1204 boxes; three follow-up rank proofs and one minimal27 equation;98 IDs and CSV,1176 incidence binding,67 exact public coordinate comparisons with two explicit unknowns; unchanged Nagao outputs and strict full small-population scores. No new point search, scanner invocation or prime trace. Completed admission histories remain available in the archive.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if name=='generic-transport' else sys.executable),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE NEW MW16 RANK27',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

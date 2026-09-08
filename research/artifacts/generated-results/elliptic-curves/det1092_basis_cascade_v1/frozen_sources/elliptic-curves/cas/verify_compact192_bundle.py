#!/usr/bin/env python3
"""Isolated exact proof and point-provenance replay, without a new point search."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/compact192-portable-v1'

def main():
    manifest=cert.read(ART/'compact192_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
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
    jobs=[
        ('selection',[str(cas/'select_compact192_unsearched.py'),'replay'],300),
        ('point-proofs',[str(cas/'certify_compact192_r17_results.py'),'--check'],1200),
        ('geometry',[str(cas/'replay_compact192_geometry.py')],1200),
        ('experiment',[str(cas/'report_compact192_experiment.py'),'--check'],600),
        ('inventory',[str(cas/'replay_inventory_v13_memory.py'),'--output',str(art/'new_high_rank_curve_index_v13_isolated_replay_v1.json')],1200),
        ('universal13',[str(cas/'report_r17_13_scaling_geometry_v2.py'),'--check'],120),
        ('symbolic13',[str(cas/'verify_r17_integral13_charts.sage')],120),
        ('first26-model',[str(cas/'export_compact192_first26_candidate.py'),'--check'],180),
        ('first26-export',[str(art/'compact192_first26_candidate.sage')],180),
    ]
    for i,row in enumerate(cert.read(local/'compact192-r17-pari-v1/protocol.json')['rows']):
        jobs.append((row['id']+'-history',[str(cas/'compact192_r17_pari_batch.py'),'replay','--index',str(i)],600))
        for label,script in [('mod2','audit_recorded_point_mod2_rank_v3.py'),('modl','audit_retained_cloud_modl.py')]:
            jobs.append((row['id']+'-'+label,[str(cas/script),'--check',str(art/('compact192_r17_'+row['id'].replace('-','_')+'_'+label+'_v1.json'))],300))
    if len(jobs)!=585:raise ArithmeticError('fixed585-stage roster differs')
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'compact192_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'rss_bytes':1610612736,'scope':'Isolated replay of192 fixed compact point admission histories and complete retained clouds modulo2,3,5, exact pre-search maps and rational geometry, selection and pinned catalogue comparisons, V13 inventory and CSV, first26 minimal model, and universal13 scaling/bad-reduction proofs. At most585 stages, one verifier at a time, no new point search or automatic retry. No whole-curve rank upper bound or universal novelty.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        executable='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python' if args[0].endswith('.sage') else sys.executable
        r=run([executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace);ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE COMPACT192',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

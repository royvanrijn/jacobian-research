#!/usr/bin/env python3
"""Extract pinned bases plus supplement and replay the new mathematics in isolation."""
from hashlib import sha256
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.supervisor import run,Limits
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/exact-parity-coordinate-portable-v1'

def main():
    manifest=cert.read(ART/'exact_parity_coordinate_evidence_v1.json');workspace=D/'workspace';folder=D/'verification'
    if workspace.exists():raise FileExistsError('preserve isolated verification')
    workspace.mkdir(parents=True);archives=[*manifest['required_base_archives'],manifest]
    for row in archives:
        archive=ROOT/row['archive']
        if cert.hashed(archive)!=row['archive_sha256']:raise ArithmeticError('archive digest changed')
        with zipfile.ZipFile(archive) as z:
            if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
            z.extractall(workspace)
    for row in [*manifest['files'],*manifest['inherited_exact_members']]:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated member bytes differ')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves';local=workspace/'artifacts/local/elliptic-curves'
    jobs=[('geometry',[str(cas/'export_r17_exact_maximum_parity_classes.py'),'--check',str(art/'r17_exact_maximum_parity_classes_v1.json')],180),
          ('rank26-tails',[str(cas/'replay_new_rank26_tails_v2.py'),'--input',str(local/'new-rank26-fixed-tails-v1/candidate-00/result.json'),'--output',str(workspace/'replayed-tail-coverage.json')],300),
          ('native-metric',[str(cas/'replay_native11952_metric49.py'),str(local/'native11952-metric49-control-v1/candidate-00/result.json')],300),
          ('native-tails',[str(cas/'replay_native11952_tails.py'),'--input',str(local/'native11952-metric49-tails-v1/candidate-00/result.json'),'--output',str(workspace/'replayed-native-coverage.json')],300),
          ('native-pari',[str(cas/'replay_native11952_pari49_v2.py'),str(local/'native11952-pari49-control-v1/candidate-00/result.json')],300),
          ('new26-pari',[str(cas/'replay_new_rank26_pari43.py'),str(local/'new-rank26-pari43-v1/candidate-00/result.json')],300),
          ('visibility',[str(cas/'audit_native11952_visibility.py'),'--check'],180),('engine-hit-sets',[str(cas/'compare_native11952_search_engines.py'),'--check'],180)]
    for family in ('11952','08f72'):
        for i in range(4):jobs.append((f'omitted-{family}-{i}',[str(cas/'replay_r17_omitted_classes.py'),str(local/f'r17-omitted-generic-classes-v1/{family}/candidate-{i:02}/result.json')],180))
    names=['new_rank26_all_retained_mod2_v1.json','new_rank26_pari43_recorded_mod2_v1.json','native11952_metric49_recorded_mod2_v1.json','native11952_all_retained_mod2_v1.json','native11952_pari49_recorded_mod2_v1.json','native11952_old43_control_mod2_v1.json','native11952_union_control_mod2_v1.json']
    names += [p.name for p in sorted(art.glob('r17_omitted_classes_recorded_mod2_*_v1.json'))]
    for name in names:jobs.append((name.removesuffix('.json'),[str(cas/'audit_recorded_point_mod2_rank_v2.py'),'--check',str(art/name)],180))
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(ART/'exact_parity_coordinate_evidence_v1.json'),'source_sha256':cert.hashed(Path(__file__).resolve()),'rss_bytes':1610612736,'jobs':[{'name':n,'args':a,'wall_seconds':s} for n,a,s in jobs],'scope':'Exact replays only in an isolated extraction; no point search, CVP search, parameter sweep or external catalogue call.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,args,seconds in jobs:
        r=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace)
        ok=r['outcome']=='completed' and r['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':r});checkpoint(folder/'ledger.json',ledger);print('PORTABLE PARITY/COORDINATE',name,ledger['rows'][-1]['status'],flush=True)
    ledger['status']='COMPLETE_DECLARED_REPLAY_ATTEMPTS';checkpoint(folder/'ledger.json',ledger)
if __name__=='__main__':main()

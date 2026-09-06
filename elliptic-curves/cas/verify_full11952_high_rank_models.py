#!/usr/bin/env python3
"""Two isolated arithmetic and Sage-export checks of the standalone new27 proof."""
from pathlib import Path
import sys
import zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/full11952-high-rank-models-portable-v1'
OUT=ART/'full11952_high_rank_models_portable_replay_v1.json'


def main():
    path=ART/'full11952_high_rank_models_evidence_v1.json';manifest=cert.read(path);workspace=D/'workspace';folder=D/'verification'
    if workspace.exists() or OUT.exists():raise FileExistsError('preserve standalone isolated new27 replay')
    if manifest['required_base_archives']:raise ArithmeticError('standalone evidence must need no base archive')
    archive=ROOT/manifest['archive']
    if cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('standalone archive hash differs')
    workspace.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        if any(Path(name).is_absolute() or '..' in Path(name).parts for name in z.namelist()):raise ArithmeticError('unsafe archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated member differs')
    jobs=[('models',sys.executable,[str(workspace/'elliptic-curves/cas/export_full11952_high_rank_models.py'),'--check']),
          ('sage-export','/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',[str(workspace/'artifacts/generated-results/elliptic-curves/new_full11952_high_rank_curves.sage')])]
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'jobs':[{'name':n,'executable':e,'args':a} for n,e,a in jobs],'wall_seconds_per_stage':120,'rss_bytes':2147483648,'scope':manifest['claim_boundary']})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,executable,args in jobs:
        s=run([executable,*args],limits=Limits(120,2147483648),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'ledger.json',ledger)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(folder/'ledger.json',ledger);raise ArithmeticError('standalone proof replay failed; no retry')
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated proof member changed during replay')
    ledger['status']='PASS';checkpoint(folder/'ledger.json',ledger)
    paths=[Path(__file__).resolve(),path,folder/'protocol.json',folder/'ledger.json']
    checkpoint(OUT,{'schema':'elliptic-curves.full11952-high-rank-models-portable-replay.v1','status':'PASS','logical_stages':2,'curves_checked':1,'independent_point_witnesses':27,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'archive_sha256':manifest['archive_sha256'],'ledger':ledger,'claim_boundary':manifest['claim_boundary']})
    print('ISOLATED NEW27 MODELS AND SAGE EXPORT: TWO STAGES PASS',flush=True)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Three isolated checks of all185 incidence and exact generic transport."""
from pathlib import Path
import sys
import zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory185-incidence-portable-v1'
OUT=ART/'inventory185_incidence_portable_replay_v1.json'


def main():
    path=ART/'inventory185_incidence_evidence_v1.json';manifest=cert.read(path);workspace=D/'workspace';folder=D/'verification'
    if workspace.exists() or OUT.exists():raise FileExistsError('preserve standalone isolated five26 replay')
    if manifest['required_base_archives']:raise ArithmeticError('standalone evidence must need no base archive')
    archive=ROOT/manifest['archive']
    if cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('standalone archive hash differs')
    workspace.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        if any(Path(name).is_absolute() or '..' in Path(name).parts for name in z.namelist()):raise ArithmeticError('unsafe archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated member differs')
    cas=workspace/'elliptic-curves/cas';art=workspace/'artifacts/generated-results/elliptic-curves'
    jobs=[('all2220-incidences',sys.executable,[str(cas/'replay_inventory185_all_incidence.py'),'--output',str(workspace/'artifacts/local/elliptic-curves/all185_incidence_replay_v1.json')]),
          ('aggregate',sys.executable,[str(cas/'certify_inventory185_incidence.py'),'--check',str(art/'inventory185_incidence_v1.json')]),
          ('generic-transport','/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',[str(cas/'audit_compact_published_r17_transport_v3.sage'),'--check',str(art/'compact_published_r17_generic_transport_v1.json')])]
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
    checkpoint(OUT,{'schema':'elliptic-curves.inventory185-incidence-portable-replay.v1','status':'PASS','logical_stages':3,'curves_checked':185,'incidence_pairs_checked':2220,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'archive_sha256':manifest['archive_sha256'],'ledger':ledger,'claim_boundary':manifest['claim_boundary']})
    print('ALL185 INCIDENCE AND GENERIC TRANSPORT: THREE ISOLATED STAGES PASS',flush=True)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""One isolated exact incidence replay from the small self-contained archive."""
from pathlib import Path
import zipfile,sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/latest8-incidence-portable-v1'

def main():
    manifest=cert.read(ART/'latest8_cross_family_evidence_v1.json');workspace=D/'workspace'
    if workspace.exists():raise FileExistsError('preserve isolated incidence replay')
    archive=ROOT/manifest['archive']
    if cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('incidence archive changed')
    with zipfile.ZipFile(archive) as z:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe incidence archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated incidence member changed')
    checkpoint(D/'protocol.json',{'manifest_sha256':cert.hashed(ART/'latest8_cross_family_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'wall_seconds':180,'rss_bytes':1073741824,'scope':'One exact Sage-free96-pair j-incidence replay from isolated extracted sources and inputs.'})
    result=run([sys.executable,str(workspace/'elliptic-curves/cas/replay_latest8_cross_family_incidence.py'),'--input',str(workspace/'artifacts/generated-results/elliptic-curves/latest8_cross_family_j_incidence_v1.json'),'--output',str(workspace/'latest8-replayed.json')],limits=Limits(180,1073741824),log_path=D/'replay.log',checkpoint_path=D/'replay.supervisor.json',cwd=workspace)
    if result['outcome']!='completed' or result['returncode']!=0:raise ArithmeticError('isolated incidence replay failed/censored')
    checkpoint(ART/'latest8_incidence_portable_replay_v1.json',{'schema':'elliptic-curves.latest8-incidence-portable-replay.v1','status':'PASS','manifest_sha256':cert.hashed(ART/'latest8_cross_family_evidence_v1.json'),'protocol':cert.read(D/'protocol.json'),'supervision':result,'log':(D/'replay.log').read_text(),'exact_replay':cert.read(workspace/'latest8-replayed.json')});print('ISOLATED LATEST8 INCIDENCE PASS',flush=True)
if __name__=='__main__':main()

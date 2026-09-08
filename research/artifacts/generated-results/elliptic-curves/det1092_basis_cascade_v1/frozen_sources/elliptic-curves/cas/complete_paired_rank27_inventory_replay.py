#!/usr/bin/env python3
"""Complete the one censored isolated topic, preserving all earlier outcomes."""
from pathlib import Path
import zipfile,sys
from hashlib import sha256
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/paired-rank27-portable-completion-v2'

def main():
    m=cert.read(ART/'paired_rank27_inventory_replay_evidence_v2.json');workspace=D/'workspace'
    if workspace.exists():raise FileExistsError('preserve isolated completion attempt')
    workspace.mkdir(parents=True)
    for b in [*m['required_base_archives'],m]:
        archive=ROOT/b['archive']
        if cert.hashed(archive)!=b['archive_sha256']:raise ArithmeticError('pinned archive changed')
        with zipfile.ZipFile(archive) as z:
            if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive entry')
            z.extractall(workspace)
    for r in m['files']:
        if cert.hashed(workspace/r['path'])!=r['sha256']:raise ArithmeticError('isolated delta member differs')
    prior_path=workspace/'artifacts/generated-results/elliptic-curves/paired_rank27_portable_replay_v1.json';prior=cert.read(prior_path);rows=prior['ledger']['rows'];names=[r['name'] for r in rows]
    if len(rows)!=35 or len(set(names))!=35 or [r['name'] for r in rows if r['status']!='PASS']!=['inventory47'] or names!=[r['name'] for r in prior['protocol']['jobs']]:raise ArithmeticError('preserved original outcomes differ')
    for r in rows:
        if sha256(prior['logs'][r['name']].encode()).hexdigest()!=r['supervision']['log_sha256']:raise ArithmeticError('original log binding differs')
    if next(r for r in rows if r['name']=='inventory47')['supervision']['outcome']!='strict_wall_timeout':raise ArithmeticError('original censor classification differs')
    checkpoint(D/'protocol.json',{'delta_manifest_sha256':cert.hashed(ART/'paired_rank27_inventory_replay_evidence_v2.json'),'prior_replay_sha256':cert.hashed(prior_path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'wall_seconds':180,'rss_bytes':1073741824,'scope':'One new isolated47-curve inventory/CSV replay using explicit MemoryFactStore. Carry the prior34 passed isolated topics by their immutable transcripts; do not claim those topics were repeated.'})
    out=workspace/'inventory47-memory-replayed.json';result=run([sys.executable,str(workspace/'elliptic-curves/cas/replay_inventory_v4_memory.py'),'--output',str(out)],limits=Limits(180,1073741824),log_path=D/'replay.log',checkpoint_path=D/'replay.supervisor.json',cwd=workspace)
    if result['outcome']!='completed' or result['returncode']!=0:raise ArithmeticError('isolated replacement failed/censored')
    replay=cert.read(out)
    if replay['status']!='PASS' or replay['curves_checked']!=47 or not replay['csv_checked']:raise ArithmeticError('replacement inventory scope differs')
    checkpoint(ART/'paired_rank27_portable_completion_v2.json',{'schema':'elliptic-curves.paired-rank27-portable-completion.v2','status':'PASS','original_passed_topics':34,'replacement_inventory_topic':'PASS','total_original_topics_covered':35,'additional_csv_check':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ART/'paired_rank27_inventory_replay_evidence_v2.json',ART/'paired_rank27_portable_replay_v1.json')},'protocol':cert.read(D/'protocol.json'),'supervision':result,'log':(D/'replay.log').read_text(),'replacement_replay':replay,'claim_boundary':'The original35-stage run had34 passes and an inventory timeout. A separately versioned isolated replay now verifies all47 inventory proofs and CSV with the tested explicit memory finite cache. All original outcomes remain retained; no failed attempt is relabelled and no new point search is claimed.'});print('ISOLATED PAIRED INVENTORY COMPLETION PASS:34 PRIOR TOPICS + INVENTORY47 AND CSV',flush=True)
if __name__=='__main__':main()

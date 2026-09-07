#!/usr/bin/env python3
"""Bind the forty completed isolated checks without repeating arithmetic."""
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/higher-r17-rank26-portable-v1';OUT=ART/'higher_r17_rank26_portable_replay_v1.json'
def main():
    if OUT.exists():raise FileExistsError('preserve portable replay report')
    m=cert.read(ART/'higher_r17_rank26_evidence_v1.json');p=cert.read(D/'verification/protocol.json');d=cert.read(D/'verification/ledger.json')
    if d['status']!='PASS' or len(d['rows'])!=40 or [r['name'] for r in d['rows']]!=[r['name'] for r in p['jobs']] or p['manifest_sha256']!=cert.hashed(ART/'higher_r17_rank26_evidence_v1.json'):raise ArithmeticError('complete bound40-stage replay required')
    if p['verifier_sha256']!=cert.hashed(ROOT/'elliptic-curves/cas/verify_higher_r17_rank26_bundle.py'):raise ArithmeticError('verifier changed')
    for r in d['rows']:
        s=r['supervision']
        if r['status']!='PASS' or s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(Path(s['log']))!=s['log_sha256']:raise ArithmeticError('successful isolated transcript differs')
    for r in [*m['files'],*m['inherited_exact_members']]:
        if cert.hashed(D/'workspace'/r['path'])!=r['sha256']:raise ArithmeticError('isolated evidence member differs')
    paths=[Path(__file__).resolve(),ART/'higher_r17_rank26_evidence_v1.json',D/'verification/protocol.json',D/'verification/ledger.json']
    checkpoint(OUT,{'schema':'elliptic-curves.higher-r17-rank26-portable-replay.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'logical_stages':40,'ledger':d,'archive':{'path':m['archive'],'sha256':m['archive_sha256']},'base_archives':m['required_base_archives'],'claim_boundary':p['scope']})
    print('ISOLATED HIGHER R17 RANK26: ALL40 STAGES PASS',flush=True)
if __name__=='__main__':main()

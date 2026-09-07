#!/usr/bin/env python3
"""Freeze/check the bounded det1092 marking and CM classification evidence."""
import argparse
import json
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
STAGES=[('det1092-marking-v1','audit_det1092_marking.sage','det1092_marking_v1.json'),
        ('det1092-marking-independent-v1','verify_det1092_marking.sage','det1092_marking_independent_v1.json'),
        ('det1092-cm-locus-v1','certify_det1092_cm_locus.sage','det1092_cm_locus_v1.json')]

def digest(p):return sha256(p.read_bytes()).hexdigest()
def build():
    records=[]
    for directory,script,artifact in STAGES:
        d=LOCAL/directory;s=ROOT/'elliptic-curves/cas'/script
        protocol=json.loads((d/'protocol.json').read_text())
        ledger=json.loads((d/'supervisor.json').read_text())
        assert protocol['source_sha256']==digest(s)
        assert ledger['outcome']=='completed' and ledger['returncode']==0
        assert ledger['log_sha256']==digest(d/'worker.log')
        assert (ART/artifact).read_bytes()==(d/'result.json').read_bytes()
        records.append(dict(directory=str(d.relative_to(ROOT)),source=str(s.relative_to(ROOT)),
                            source_sha256=digest(s),artifact=artifact,artifact_sha256=digest(ART/artifact),
                            protocol_sha256=digest(d/'protocol.json'),ledger_sha256=digest(d/'supervisor.json'),
                            completed_wall_seconds=ledger['wall_seconds']))
    marking=json.loads((ART/STAGES[0][2]).read_text())
    independent=json.loads((ART/STAGES[1][2]).read_text())
    cm=json.loads((ART/STAGES[2][2]).read_text())
    assert independent['input_sha256']==digest(ART/STAGES[0][2])
    assert independent['checker_sha256']==records[1]['source_sha256']
    assert cm['checker_sha256']==records[2]['source_sha256']
    assert marking['projectively_stable_labels']==independent['stable_labels']==[1,546]
    assert cm['surviving_discriminants']==[-67,-163] and cm['total_rational_CM_points']==16
    return dict(schema='det1092.marking-followup.v1',status='PASS',stages=records,
                total_supervised_seconds=sum(r['completed_wall_seconds'] for r in records),
                full_projective_stable_curve='X(546)/<w546>',
                marking_arithmetic_independently_replayed=True,
                rational_CM_discriminants=[-67,-163],rational_CM_point_count=16,
                CM_completeness_uses_named_external_theorems=True,
                CM_coordinates='UNKNOWN',individual_non_CM_point='UNKNOWN',
                explicit_K3_equation='UNKNOWN',new_searched_fibres=0,
                boundary='No new parent equation, section basis, record recovery, or near-record curve.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    if not a.check:
        for directory,_,artifact in STAGES:
            with (ART/artifact).open('xb') as f:f.write((LOCAL/directory/'result.json').read_bytes())
    result=build();out=ART/'det1092_marking_followup_v1.json'
    if a.check:assert result==json.loads(out.read_text())
    else:
        with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'total_supervised_seconds':result['total_supervised_seconds']}))

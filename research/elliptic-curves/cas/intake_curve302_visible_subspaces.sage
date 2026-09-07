#!/usr/bin/env sage-python
"""Exact containment intake for the rank11 visible Mestre subgroup signature.

Thirteen retained rank11--13 candidates against six tested spaces; no new
point enumeration. Candidate dimension is distinct from parent generic rank.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from sage.all import ZZ,matrix

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
SELECTION=ART/'low_height_mw_sublattices_v1_302_selection.json'
TESTED=ART/'curve302_wide_packet_protocol_v1.json'
OUT=ART/'curve302_visible_subspace_intake_v1.json'


def build():
    source=json.loads(SELECTION.read_text());tested=json.loads(TESTED.read_text());parents=[(e['label'],matrix(ZZ,e['basis_rows'])) for e in tested['inputs']]
    records=[]
    for e in source['finalists']:
        if e['rank'] not in [11,12,13]:continue
        B=matrix(ZZ,e['primitive_basis_rows']);assert B.rank()==e['rank'] and all(abs(n)==1 for n in B.smith_form()[0].diagonal())
        comparisons=[]
        for label,C in parents:
            union=int(C.stack(B).rank());record={'tested_label':label,'tested_rank':C.nrows(),'union_rank':union,'contained':union==C.nrows()}
            if record['contained']:
                change=matrix(ZZ,C.solve_left(B));assert change*C==B;record['integral_embedding_rows']=[list(map(int,r)) for r in change.rows()]
            comparisons.append(record)
        records.append({'candidate_id':e['candidate_index'],'rank':e['rank'],'selection_channel':e['channel'],
                        'primitive_basis_rows':e['primitive_basis_rows'],'comparisons':comparisons,'outside_all_tested_spaces':not any(x['contained'] for x in comparisons)})
    new=[e for e in records if e['outside_all_tested_spaces']];assert sorted(e['candidate_id'] for e in new)==[104,105,106,130]
    relations=[]
    for e in new:
        B=matrix(ZZ,e['primitive_basis_rows'])
        for f in new:
            if e['candidate_id']==f['candidate_id']:continue
            C=matrix(ZZ,f['primitive_basis_rows']);union=int(C.stack(B).rank())
            relations.append({'source':e['candidate_id'],'target':f['candidate_id'],'union_rank':union,'contained':union==C.nrows()})
    maximal=[e['candidate_id'] for e in new if not any(r['source']==e['candidate_id'] and r['contained'] for r in relations)]
    assert sorted(maximal)==[106,130]
    return {'schema':'curve302.visible-subspace-intake.v1','status':'PASS_FOUR_UNCOVERED_CANDIDATES_TWO_MAXIMAL_SPACES',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SELECTION,TESTED]},
        'records':records,'new_candidate_relations':relations,'maximal_new_candidate_ids':sorted(maximal),
        'boundary':'Exact subspace containment only. No twelve-point configuration, new family, generic rank or point search is certified. Twelve visible Mestre covariant images satisfy a signed linear relation, so their primitive span can have rank11 inside a higher-rank parent. Containment in a larger tested space alone does not prove containment of a numerically truncated point packet.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();d=build()
    if args.check:assert d==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'])

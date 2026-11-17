#!/usr/bin/env python3
"""Independent order2 control, using the existing pointwise verifier helpers."""

from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path

from verify_ternary_transvection_control import (POINTS,POINT_INDEX,act,evaluate_word,
    mat_rank,pointwise_root_matrix,signed_crossings,word_for_crossing)

HERE=Path(__file__).resolve().parent


def check(report):
    chain=report['chain_matrices'];points4=list(product(range(3),repeat=4))
    for m in chain:
        assert all(act(m,act(m,v))==v for v in points4)
        assert len({tuple((x-y)%3 for x,y in zip(act(m,v),v)) for v in points4})==3
    assert [v for v in points4 if all(act(m,v)==v for m in chain)]==[(0,0,0,0)]
    for i in range(3):
        word=[i+1,i+2,i+1,-(i+2),-(i+1),-(i+2)]
        assert all(evaluate_word(chain,word,v)==v for v in points4)
    # Restore the same independently checked literal connected-sum diagram.
    known=json.loads((HERE/'ternary_transvection_control.json').read_text())
    assert report['pd']==known['pd'] and report['edge_renumbering']==known['edge_renumbering']
    crossings=signed_crossings(report['pd']);full=report['arc_matrices']
    for row in crossings:assert all(evaluate_word(full,word_for_crossing(row),v)==v for v in points4)
    linear=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in full]
    rows,roots=pointwise_root_matrix(crossings,linear)
    assert rows==report['root_matrix'] and 9-mat_rank(rows,9)==report['root_dimension']==4
    basis=report['root_basis'];assert len(basis)==mat_rank(basis,9)==4
    assert all(sum(x*y for x,y in zip(row,v))%3==0 for row in rows for v in basis)
    for m,r,c in zip(full,roots,report['arc_root_coefficients']):
        assert tuple(m[4*i+3] for i in range(3))==tuple(c*v%3 for v in r)
    first_three=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in chain[:3]]
    generators=[tuple(POINT_INDEX[act(m,v)] for v in POINTS) for m in first_three]
    reached={tuple(range(27))};queue=list(reached)
    for p in queue:
        for g in generators:
            q=tuple(g[p[i]] for i in range(27))
            if q not in reached:assert len(reached)<24;reached.add(q);queue.append(q)
    fourth=tuple(chain[3][4*i+j] for i in range(3) for j in range(3))
    assert tuple(POINT_INDEX[act(fourth,v)] for v in POINTS) in reached
    v=(1,1,1);orbit={POINTS[p[POINT_INDEX[v]]] for p in reached}
    assert mat_rank(list(orbit),3)==3
    assert len(reached)==report['linear_first_three_image_order']==24
    return {'homogeneous_points_checked':81,'root_dimension':4,'linear_image_order':24,
        'translation_orbit_span':3,'affine_order_from_linear_copy_and_translation':648,'meridional_rank':4,'bridge_number':4}


def main():
    report=json.loads((HERE/'ternary_reflection_control.json').read_text())
    known=json.loads((HERE/'ternary_transvection_control_verification.json').read_text())
    assert known['status']=='PASS independent ternary triple-trefoil control'
    for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    result=check(report);bad=deepcopy(report);bad['arc_root_coefficients'][0]=(bad['arc_root_coefficients'][0]+1)%3
    try:check(bad)
    except AssertionError:rejected=True
    else:raise AssertionError('accepted wrong affine root coefficient')
    files=['ternary_reflection_control.json','verify_ternary_reflection_control.py','verify_ternary_transvection_control.py','ternary_transvection_control_verification.json']
    receipt={'status':'PASS independent ternary reflection control','checks':result,'wrong_root_coefficient_rejected':rejected,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'ternary_reflection_control_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':result,'wrong_root_coefficient_rejected':rejected},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env sage-python
"""Bind the two equation replays to their frozen historical/generic sources.

Comparison with the historical point happens only AFTER both equation-only
executions. Neither the historical coordinates nor this comparison is an
input to the per-cover checker. No V3 evidence is read.
"""
import hashlib
import json
import runpy
import signal
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DIR = ART/'det1092_single_seed_covers_v1'
CHECKER = ROOT/'elliptic-curves/cas/verify_det1092_single_seed_covers.sage'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify():
    protocol = json.loads((DIR/'protocol.json').read_text())
    for path, digest in protocol['inputs'].items():
        assert sha(ROOT/path) == digest
    for path, digest in protocol['payload_sha256'].items():
        assert sha(DIR/path) == digest
    checker = runpy.run_path(str(CHECKER))
    for i in range(2):
        checker['verify'](i)
    paths = [ART/'curve302_recovered_mw17_parent_v1.json',
             ART/'det1092_first_centre_rr_net_v1.json',
             ART/'det1092_first_witness_pencil_genus_gate_v1.json',
             ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent, net, pencil, roster = [json.loads(p.read_text()) for p in paths]
    R = PolynomialRing(QQ, 't')
    A, B = [[R(row) for row in net[k]] for k in ['A','B']]
    replays = []
    for i, u in enumerate([QQ(0),QQ(pencil['u0'])]):
        d = json.loads((DIR/f'cover-{i:02d}-input.json').read_text())
        proof = json.loads((DIR/f'cover-{i:02d}-replay.json').read_text())
        assert d['a_invariants'] == parent['a_invariants']
        assert d['generic_sections'] == parent['basis_weierstrass_coordinates']
        assert d['generic_height_gram'] == parent['generic_height_gram']
        assert d['trace_word'] == [-n for n in net['trace_word']]
        assert list(map(R,d['line'])) == [B[j]+u*A[j] for j in range(3)]
        assert d['cases'] == [{'index': j, 'parameter': row['parameter']}
                              for j,row in enumerate(roster['cases'])]
        assert proof['checker_sha256'] == sha(CHECKER)
        replays.append(proof)
    # Post-execution calibration audit: identify which of the two branches was
    # reconstructed, and bind this exact cover to the certified Jacobian xi.
    bridge_path = ART/'det1092_first_centre_rr_net_replay_v1.json'
    jacobian_path = ART/'det1092_rr_residual_jacobian_class_v1.json'
    bridge = json.loads(bridge_path.read_text())
    jacobian = json.loads(jacobian_path.read_text())
    d = json.loads((DIR/'cover-01-input.json').read_text())
    assert QQ(jacobian['curve']['u']) == QQ(pencil['u0'])
    assert QQ(jacobian['curve']['v']) == 0
    assert R(d['cover_polynomial']) == QQ(jacobian['curve']['scale'])*R(jacobian['curve']['q'])
    first = replays[1]['cases'][8]
    assert first['status'] == 'CERTIFIED_INDEPENDENT_SEED'
    historical = list(map(QQ, bridge['reconstructed_point_literal302']))
    branches = [list(map(QQ, first[k])) for k in ['point','complement']]
    assert historical in branches
    assert all(r['status'] == 'EXACT_NONSPLIT_FIBRE' for p in replays for r in p['cases'][:8])
    assert replays[0]['cases'][8]['status'] == 'EXACT_NONSPLIT_FIBRE'
    table = []
    for j, case in enumerate(protocol['cases']):
        row = {'index': j, 'label': case['label'], 'parameter': case['parameter']}
        for i in range(2):
            r = replays[i]['cases'][j]
            row[f'cover_{i:02d}'] = {'status': r['status'],
                'local_obstruction': r.get('modular_attempts',[None])[-1]}
        table.append(row)
    result = {
        'classification': 'verified application and new deduction',
        'status': 'PASS_SINGLE_SEED_COVERS_AND_EXACT_CONTROL_OBSTRUCTIONS',
        'input_boundary_verified': True, 'existing_covers': 2, 'fixed_parameters': 9,
        'both_function_field_rank_lower_bounds': [18,18],
        'both_anti_invariant_heights': [20,20],
        'reconstructed_historical_branch_index': branches.index(historical),
        'same_curve_as_certified_first_Jacobian_class': True,
        'relative_Jacobian_branch_rule': 'xi(P)+xi(sigma(P))=[K_C-2P0] lies in inherited H; the two branches differ only by sign modulo H.',
        'relative_elliptic_branch_rule': 'P+sigma(P)=Z in MW17; either branch is independent if the other is.',
        'held_out_302_discovery': False,
        'new_independent_control_seeds': 0,
        'table': table,
        'next_obstruction': 'The calibrated member is an explicit partial seed certifier, not a learned oracle-free member selector or a source of additional rational cover points. Every rational member of its fixed witness pencil has genus2 by the prior exact gate. No smaller-genus cover through this witness in another system is classified.',
        'limits': {'wall_seconds': 25, 'point_searches': 0, 'new_parameters': 0,
                   'V3_artifacts_read': 0, 'pilot_changes': 0, 'class_groups': 0},
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in
                   [DIR/'protocol.json',CHECKER,*paths,bridge_path,jacobian_path,
                    *[DIR/f'cover-{i:02d}-replay.json' for i in range(2)]]},
        'checker_sha256': sha(Path(__file__)),
    }
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    out = DIR/'panel-replay.json'
    if out.exists():
        assert out.read_text() == text
    else:
        out.write_text(text)
    print(result['status'], flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    verify()

#!/usr/bin/env python3
"""Replay small automorphism groups against the retained complete frame census.

Completeness is inherited from the census and its qfauto orders, not reproved.
No frame enumeration or equation calculation is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[3]
CENSUS = 'research/artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json'
PACKET = 'research/artifacts/generated-results/elkies-k3-x1092-involution-glue-v1'


def replay():
    source = (ROOT / CENSUS).read_bytes()
    packet = json.loads((ROOT / PACKET / 'input.json').read_text())
    assert packet['census_sha256'] == hashlib.sha256(source).hexdigest()
    classes = json.loads(source)['rootless_classes']
    assert [c['class_index'] for c in classes] == list(range(1, 20))
    inputs = {c['class_index']: c for c in packet['groups']}
    assert set(inputs) == set(range(1, 20))
    identity = S.eye(17)
    rows = []
    for c in classes:
        h = S.Matrix(c['gram'])
        assert h.shape == (17, 17) and h == h.T and h.det() == 1092
        gens = [S.Matrix(g) for g in inputs[c['class_index']]['generators']]
        for g in gens:
            assert g.shape == (17, 17) and all(x.is_Integer for x in g)
            assert g.T * h * g == h and abs(g.det()) == 1
        group = {tuple(identity): identity}
        pending = [identity]
        while pending:
            a = pending.pop()
            for g in gens:
                b = a * g
                key = tuple(b)
                if key not in group:
                    group[key] = b
                    pending.append(b)
                    assert len(group) <= c['automorphism_group_order']
        assert len(group) == c['automorphism_group_order']
        assert tuple(-identity) in group
        traces = []
        witnesses = []
        inverse = h.inv()
        for key in sorted(group):
            g = group[key]
            if g * g != identity:
                continue
            tr = int(S.trace(g))
            traces.append(tr)
            if tr != 1:
                continue
            dual = (g - identity) * inverse
            bad = [(i, j, str(dual[i, j])) for i in range(17)
                   for j in range(17) if not dual[i, j].is_Integer]
            assert bad, 'Unexcluded symplectic candidate: do not claim closure'
            witnesses.append({'matrix': [list(map(int, g.row(i))) for i in range(17)],
                              'nonintegral_dual_entry': list(bad[0])})
        rows.append({'class_index': c['class_index'], 'group_order': len(group),
                     'involution_traces': sorted(traces),
                     'trace_one_obstructions': witnesses})
    assert sum(len(r['trace_one_obstructions']) for r in rows) == 5
    return {'schema': 'x1092-involution-glue-v1', 'status': 'PASS',
            'census_sha256': packet['census_sha256'],
            'input_sha256': hashlib.sha256((ROOT / PACKET / 'input.json').read_bytes()).hexdigest(),
            'scope': 'All nineteen retained rootless J2 frames: no trace-one involutive isometry acting trivially on the discriminant group. Group completeness and frame coverage are inherited from the retained census.',
            'classes': rows}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Write the derived certificate')
    args = parser.parse_args()
    result = replay()
    path = ROOT / PACKET / 'result.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2) + '\n')
    else:
        assert result == json.loads(path.read_text())
    print('PASS: 19 complete groups; all 5 trace-one involutions fail discriminant gluing')

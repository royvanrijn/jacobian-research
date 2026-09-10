#!/usr/bin/env sage-python
"""Rank new J2 types by exact shared-core geometry, without fibre outcomes.

The ranking is a prospective equation-cost heuristic, not a proof of nefness,
rational marked-U realization, or a prediction of specialization rank.
"""
import json
from pathlib import Path
import runpy
from sage.all import ZZ, matrix

HERE = Path(__file__).parent
m = runpy.run_path(str(HERE/'classify_det1092_rootless_j2_pruned.sage'))
o, ROOT = m['old'], m['ROOT']
CENSUS = ROOT/'artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json'
PACKET = m['DEFAULT_WORK']/'anchor-16.json'
OUTPUT = CENSUS.with_name('det1092_frame_realization_priority_v1.json')


def build():
    census, packet = json.loads(CENSUS.read_text()), json.loads(PACKET.read_text())
    assert census['status'] == 'PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION'
    assert o['digest'](PACKET) == census['anchor_packet_hashes'][PACKET.name]
    known_class, = census['accounting']['known_frame_matching_class_indices']
    groups = {}
    for x in census['embedding_class_assignments']:
        assert x['anchor_number'] == 16
        groups.setdefault(x['class_index'], []).append(x['embedding_index'])
    cat = json.loads(o['CATALOG'].read_text())['rooted_niemeier_lattices']
    G = matrix(ZZ, next(x for x in cat if x['label'] == packet['anchor']['niemeier'])['gram'])
    out = []
    for klass, targets in sorted(groups.items()):
        if klass == known_class:
            continue
        choices = []
        for i in groups[known_class]:
            left = packet['embeddings'][i]
            B = matrix(ZZ, left['auxiliary_basis_in_ambient'])
            for j in targets:
                right = packet['embeddings'][j]
                if left['sixth_index'] != right['sixth_index']:
                    continue
                D = matrix(ZZ, right['auxiliary_basis_in_ambient'])
                assert B[:6, :] == D[:6, :]
                union = B.stack(D)
                K = (union * G).right_kernel_matrix()
                H = K * G * K.transpose()
                assert K.nrows() == 16 and o['primitive'](K)
                assert K * G * B.transpose() == 0 and K * G * D.transpose() == 0
                # This is a primitive common codimension-one sublattice of
                # both frames, before any identification of the marked NS.
                size = max(abs(int(v)).bit_length() for v in union.list())
                score = (int(H.det()), size, i, j)
                choices.append((score, {'known_embedding_index': i, 'new_embedding_index': j,
                    'shared_sixth_index': left['sixth_index'], 'common_core_rank': 16,
                    'common_core_determinant': int(H.det()), 'common_core_gram': o['rows'](H),
                    'common_core_basis_in_niemeier': o['rows'](K), 'maximum_auxiliary_coordinate_bits': size}))
        if not choices:
            out.append({'class_index': klass, 'status': 'UNKNOWN_SHARED_CORE', 'score': None})
            continue
        score, witness = min(choices, key=lambda item: item[0])
        out.append({'class_index': klass, 'status': 'EXACT_SHARED_RANK16_CORE',
            'score': list(score[:2]), 'eligible_shared_core_pairs': len(choices), 'witness': witness,
            'marked_nef_U': 'UNKNOWN', 'rational_equation': 'UNKNOWN', 'generic_MW17_equation_certificate': 'UNKNOWN'})
    out.sort(key=lambda x: (x['score'] is None, x['score'] or [], x['class_index']))
    return {'schema': 'det1092.frame-realization-priority.v1', 'status': 'PASS_PROSPECTIVE_LATTICE_PRIORITY',
        'bindings': {str(p.relative_to(ROOT)): o['digest'](p) for p in (CENSUS, Path(__file__))},
        'ordering': 'Prefer a common rank16 core; minimize its determinant, then auxiliary-coordinate bits, then class index. This is an equation-tractability heuristic only.',
        'selected_class_index': out[0]['class_index'], 'ranked_new_types': out,
        'next_required_step': 'Transport the selected shared core to the known rational NS marking and construct a primitive nef U with an effective rational zero. Then compile the pencil and certify its generic rational MW17 sections.',
        'specialization_search': 'BLOCKED_UNTIL_EXACT_RATIONAL_EQUATION_AND_GENERIC_MW17_CERTIFICATE',
        'boundary': 'No rank labels, exceptional points, parameters, or specialization outcomes entered this ranking. Shared Niemeier cores do not themselves give marked NS isometries, nef U, equations, or arithmetic fibrations.'}


if __name__ == '__main__':
    result = build()
    raw = o['canonical'](result)
    if OUTPUT.exists():
        assert OUTPUT.read_bytes() == raw
    else:
        with OUTPUT.open('xb') as f:
            f.write(raw)
    print(result['status'], 'selected class', result['selected_class_index'])

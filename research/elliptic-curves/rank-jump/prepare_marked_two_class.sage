#!/usr/bin/env sage-python
"""Freeze the equation-only marking for the bounded MW16-05 carrier bank.

Uses the existing integral U-complement construction, including its A1 root.
No exceptional points, point clouds, or oracle square roots are inputs.
"""
import hashlib
import json
from pathlib import Path
import resource
import time

from sage.all import ZZ, QQ, matrix, vector, block_diagonal_matrix, pari
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'artifacts/local/elliptic-curves/marked-two-class-v1'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, data):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def rows(m):
    return [list(map(int, row)) for row in m.rows()]


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (120, 125))
    start = time.process_time()
    if (OUT / 'protocol.json').exists():
        raise FileExistsError('preserve the frozen run')
    native_path = ROOT / 'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
    template_path = ROOT / 'elliptic-curves/data/a1_mw16_family_template_v1.json'
    atlas_path = ROOT / 'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json'
    native = json.loads(native_path.read_text())
    template = json.loads(template_path.read_text())
    presentation = next(p for p in template['presentations'] if p['presentation_id'] == 'a1-presentation-07')
    family = next(p for p in json.loads(atlas_path.read_text())['families'] if p['fibration_id'] == 'a1-fibration-05')
    assert presentation['fibration_id'] == family['fibration_id']
    assert family['presentation_id'] == presentation['presentation_id']
    marking = presentation['source_marking']
    M = matrix(ZZ, native['sections']['height_gram'])
    C = matrix(ZZ, native['sections']['coordinate_matrix_in_compiled_frame'])
    frame = matrix(ZZ, native['frame_certificate']['frame_gram'])
    assert abs(C.det()) == 1 and M == C * frame * C.transpose()
    assert M.det() == 948 and M.is_positive_definite()
    Q = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -M)
    oldF = vector(ZZ, [1, 0] + [0]*17)
    oldO = vector(ZZ, [-1, 1] + [0]*17)
    def section(w):
        w = vector(ZZ, w)
        return vector(ZZ, [(w*M*w-2)//2, 1] + list(w))
    w = vector(ZZ, marking['trace_section_basis_w'])
    assert w*M*w == 8
    F = vector(ZZ, [2, 2] + list(w))
    O = section(marking['new_zero_source_section_basis_coordinates'])
    assert F*Q*F == 0 and O*Q*O == -2 and F*Q*O == 1
    # At the old lambda=infinity the pencil is oldO + oldP_w.
    components = [oldO, section(w)]
    assert sum(components) == F
    assert sorted(int(c*Q*O) for c in components) == [0, 1]
    root = next(c for c in components if c*Q*O == 0)
    complement = matrix(ZZ, [list(F*Q), list((O+F)*Q)]).right_kernel_matrix()
    G = -(complement*Q*complement.transpose())
    change = G.LLL_gram().transpose()
    basis = change*complement
    G = -(basis*Q*basis.transpose())
    transport = matrix(ZZ, [list(F), list(O+F)] + rows(basis))
    assert abs(transport.det()) == 1
    assert transport*Q*transport.transpose() == block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -G)
    assert G.det() == 948 and G.nrows() == 17 and G.is_positive_definite()
    # Enumerating roots only is a small marking verification, not either shell.
    roots = pari(G).qfminim(2).sage()
    assert int(roots[0]) == 2
    r = vector(ZZ, list(root*transport.inverse())[2:])
    assert r*G*r == 2 and r*basis == root
    generic = [section(v) for v in marking['generic_source_section_basis_coordinates']]
    assert all(S*Q*F == 1 and S*Q*S == -2 for S in generic)
    shioda = [S-O-(2+S*Q*O)*F+QQ(S*Q*root)/2*root for S in generic]
    height = -matrix(QQ, shioda)*Q*matrix(QQ, shioda).transpose()
    assert height == matrix(QQ, family['generic_height_gram']) and height.det() == 474
    # Both signs of the generic sections are effective, exactly by group law.
    negatives = [2*O+2*(2+S*Q*O)*F-(S*Q*root)*root-S for S in generic]
    assert all(S*Q*S == -2 and S*Q*F == 1 for S in negatives)
    geometry = {'schema': 'marked-two-class.geometry.v1', 'sage_version': SAGE_VERSION,
        'family': family, 'presentation': presentation, 'native_generic_model': native,
        'source_ns_gram': rows(Q), 'native_mw_gram': rows(M),
        'fibre': list(map(int, F)), 'zero': list(map(int, O)),
        'nonidentity_I2_component': list(map(int, root)),
        'I2_components': [list(map(int, c)) for c in components],
        'frame_basis_in_source_ns': rows(basis), 'frame_gram': rows(G),
        'unimodular_transport': rows(transport), 'inverse_transport': rows(transport.inverse()),
        'generic_divisors': [list(map(int, S)) for S in generic],
        'effective_wall_divisors': [list(map(int, S)) for S in [O, *components, *generic, *negatives]],
        'marking_checks': {'frame_rank': 17, 'frame_determinant': 948, 'root_count': 2,
            'height_rank': 16, 'height_determinant': 474, 'integral_transport_determinant': int(transport.det())},
        'boundary': 'Full integral marking inherited from the certified generic native model and sanitized MW16 presentation; root, glue, section order and integral transports rechecked. No new full model-reconstruction replay.'}
    write('geometry.json', geometry)
    paths = [native_path, template_path, atlas_path, Path(__file__),
        ROOT/'elkies-k3/scripts/certify_icarm_norm8_low_root_strata.sage',
        ROOT/'elkies-k3/scripts/construct_elkies_2026_bisections.sage',
        ROOT/'elliptic-curves/cas/export_compact_mw16_atlas.sage',
        ROOT/'elliptic-curves/cas/visibility_lattice_fast.py',
        ROOT/'elliptic-curves/cas/visibility_lattice_v2.py',
        ROOT/'elliptic-curves/rank-jump/enumerate_marked_two_class.py']
    sources = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    for p in paths:
        destination = OUT/'source-snapshots'/p.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open('xb') as stream:
            stream.write(p.read_bytes())
    write('protocol.json', {'schema': 'marked-two-class.protocol.v1', 'sources': sources,
        'geometry_sha256': digest(OUT/'geometry.json'), 'control_parameter': '3/17',
        'scope': 'Finite class-labelled smooth rational bisection bank on MW16-05, not a rank correlation study.',
        'separation': 'Equation-only geometry input; frozen columns and generic Kummer data admitted in the later label stage. Instruction and audited input boundaries, no OS isolation claim.',
        'aggregate_cpu_seconds': 14400, 'maximum_workers': 1, 'rss_bytes': 8589934592,
        'enumeration': {'total_nodes': 2000000, 'cpu_seconds': 3600,
            'schedule': 'One exact integer LDL depth-first closed ellipsoid of norm <=14; retain both signed norm10 and norm14 vectors. Increasing integer coordinates, fixed frozen reduced basis. A node stop leaves BOTH shells incomplete.',
            'shells': [10, 14], 'checkpoint_nodes': 100000,
            'candidate_order': 'After known effective-wall rejections, increasing native fibre degree, native O intersection, sum of absolute native NS coordinates, then lexicographic native NS coordinates. This is a fixed generic RR-cost proxy; no splitting or labels enter ordering.',
            'compilation_prefix_per_shell': 128, 'interrupted_order_boundary': 'The ordered first128 of the retained prefix, not the first128 of a completed shell.'},
        'compilation_cpu_seconds_per_candidate': 30,
        'labels_per_nonzero_class': 8, 'maximum_pairs': 192,
        'genus_one_base_probes': 8, 'base_probe_cpu_seconds': 120,
        'fresh_fibres': 16, 'automatic_enlargement': False,
        'preparation_cpu_seconds': time.process_time()-start,
        'status': 'FROZEN_BEFORE_ENUMERATION'})
    print(json.dumps(geometry['marking_checks']), flush=True)


if __name__ == '__main__':
    main()

"""Test containment of the remaining saved order in the joined RR order.

Two fixed cases, no number-field initialization or factorization.
"""
import hashlib
import json
import runpy
import signal
from pathlib import Path
from sage.all import QQ, ZZ, matrix

BASE = Path(__file__).resolve().parent
env = runpy.run_path(str(BASE / 'verify_det1092_rr_global_pair.sage'))
ROOT, ART, R, coords = [env[k] for k in ['ROOT', 'ART', 'R', 'coords']]
x = R.gen()
signal.alarm(60)
inputs = {}


def read(p):
    inputs[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())


results = []
for i in [8, 9]:
    previous = ART / 'det1092_rr_global_pair_v2'
    old = read(previous / ('case-%02d-order.json' % i))
    model = read(previous / ('case-%02d-model.json' % i))
    reduced = read(ART / 'det1092_rr_global_pair_v3' / ('case-%02d-reduced-generator.json' % i))
    joined = read(ART / 'det1092_rr_order_join_v1' / ('case-%02d.json' % i))
    f = R(old['defining_polynomial'])
    theta = R(old['old_root_in_new_field'])
    primitive = R(model['primitive_sextic'])
    env['check_isomorphism'](primitive, f, theta)
    powers = matrix(QQ, [coords(theta**j % f) for j in range(6)]).transpose()
    inverse_root = R(list(powers.inverse()*coords(x)))
    assert inverse_root(theta) % f == x
    target = R(reduced['defining_polynomial'])
    target_theta = R(reduced['monic_root_in_reduced_field']) / primitive[6]
    root = inverse_root(target_theta) % target
    env['check_isomorphism'](f, target, root)
    assert target == R(joined['defining_polynomial'])
    basis = [R(b)(root) % target for b in old['order_basis']]
    env['check_order'](target, basis, old['order_discriminant'])
    new_basis = [R(b) for b in joined['basis']]
    env['check_order'](target, new_basis, joined['order_discriminant'])
    B = matrix(QQ, [coords(b) for b in new_basis]).transpose()
    inclusion = B.inverse()*matrix(QQ, [coords(b) for b in basis]).transpose()
    contained = all(v in ZZ for v in inclusion.list())
    record = {'case_index': i, 'third_order_contained': contained,
              'old_field_root_in_joined_field': list(map(str, root.list())),
              'inclusion_matrix': [[str(v) for v in row] for row in inclusion.rows()]}
    if contained:
        index = abs(ZZ(inclusion.det()))
        assert ZZ(old['order_discriminant']) == ZZ(joined['order_discriminant'])*index**2
        record['index'] = str(index)
    results.append(record)
print(json.dumps({'status': 'PASS_EXACT_CONTAINMENT_AUDIT', 'cases': results,
    'inputs': inputs, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'dependency_sha256': hashlib.sha256((BASE/'verify_det1092_rr_global_pair.sage').read_bytes()).hexdigest(),
    'wall_limit_seconds': 60,
    'boundary': 'No maximality or global squareclass certification. If contained, the third saved order adds no integrality.'}, indent=2))

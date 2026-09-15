"""Verify joined integral orders with Newton traces, without lattice construction."""
import hashlib
import json
import runpy
import signal
from pathlib import Path
from sage.all import QQ, ZZ, matrix

BASE = Path(__file__).resolve().parent
prior = runpy.run_path(str(BASE / 'verify_det1092_rr_global_pair.sage'))
R, coords, check_order = [prior[k] for k in ['R', 'coords', 'check_order']]
ROOT, ART = prior['ROOT'], prior['ART']
OUT = ART / 'det1092_rr_order_join_v1'
signal.alarm(60)
rows = []
for i in [8, 9]:
    path = OUT / ('case-%02d.json' % i)
    d = json.loads(path.read_text())
    assert hashlib.sha256((BASE/'combine_rr_pair_orders.sage').read_bytes()).hexdigest() == d['source_sha256']
    for p, digest in d['inputs'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == digest
    old = ART / 'det1092_rr_global_pair_v3'
    read = lambda suffix: json.loads((old/('case-%02d-%s.json' % (i,suffix))).read_text())
    binary, reduced, partial = [read(s) for s in ['binary-order','reduced-generator','order']]
    f, alpha = R(reduced['defining_polynomial']), R(reduced['monic_root_in_reduced_field'])
    prior['check_isomorphism'](R(binary['monic_polynomial']), f, alpha)
    assert f == R(d['defining_polynomial']) == R(partial['defining_polynomial'])
    basis = [R(b) for b in d['basis']]
    bits = check_order(f, basis, d['order_discriminant'])
    inverse = matrix(QQ, [coords(b) for b in basis]).transpose().inverse()
    for old_basis, disc, index in [
        ([R(b)(alpha)%f for b in binary['basis']], binary['order_discriminant'], d['input_order_indices'][0]),
        ([R(b) for b in partial['order_basis']], partial['order_discriminant'], d['input_order_indices'][1])]:
        check_order(f, old_basis, disc)
        containment = inverse*matrix(QQ, [coords(b) for b in old_basis]).transpose()
        assert all(v in ZZ for v in containment.list())
        assert abs(containment.det()) == ZZ(index)
        assert ZZ(disc) == ZZ(d['order_discriminant'])*ZZ(index)**2
    assert bits == d['order_discriminant_bits']
    rows.append({'case_index': i, 'order_discriminant_bits': bits,
                 'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
print(json.dumps({'status': 'PASS', 'cases': rows,
    'boundary': 'Integral ring closure, Newton trace discriminants and containment only; no maximality or Selmer completeness.',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'dependency_sha256': hashlib.sha256((BASE/'verify_det1092_rr_global_pair.sage').read_bytes()).hexdigest()}, indent=2))

#!/usr/bin/env sage-python
"""Independent source/incidence/frame/cycle replay for the corrected four cases.

No constructor, admission producer, campaign, or V3 imports. --write creates
immutable replay receipts; the default verifies exact reproducibility.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from sage.all import QQ, PolynomialRing
from sage.env import SAGE_VERSION
from replay_split_admission import ReplayFrame, add, multiply, word_point, require

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2'
ROSTER = [(4, '2', 82931), (17, '-4/3', 65035), (43, '-5/7', 82931), (46, '2/7', 30223)]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(index):
    protocol = read(OUT/'protocol.json')
    require(protocol['roster'] == [dict(index=i, parameter=t, mask=m) for i, t, m in ROSTER], 'exact four-case roster')
    require(protocol['max_successful_halves'] == 8, 'declared halving bound')
    for name, info in protocol['inputs'].items():
        require(sha(OUT/name) == info['sha256'], 'immutable input seal')
    _, parameter, mask = next(row for row in ROSTER if row[0] == index)
    inputs = OUT/'inputs'; old = inputs/str(index); folder = OUT/'cases'/str(index)
    parent = read(inputs/'parent.json')
    plan, manifest = read(inputs/'plan.json'), read(inputs/'manifest.json')
    require(sha(inputs/'plan.json') == manifest['plan_sha256'], 'original plan seal')
    require(sha(inputs/'parent.json') == manifest['files']['seed-inputs/parent.json'], 'original parent seal')
    require(dict(index=index, parameter=parameter, control=index % plan['control_every'] == 0)
            in plan['parameters'], 'original parameter selection')
    seed, request, admission = map(read, [old/'seed.json', old/'request.json', old/'admission.json'])
    require(seed['parameter'] == request['parameter'] == parameter and len(seed['points']) == 17, 'M17 parameter/prefix')
    require(seed['generic_rank'] == seed['rank_lower_bound'] == seed['proof']['rank_lower_bound'] == 17, 'original M17 rank')
    require(request['parent_sha256'] == sha(inputs/'parent.json')
            and request['incidence_sha256'] == sha(old/'incidence.json'), 'original request seals')
    require(admission['request_sha256'] == sha(old/'request.json')
            and admission['seed_sha256'] == sha(old/'seed.json'), 'original admission seals')
    prior_supervisor = read(old/'supervision-receipt.json')
    require(prior_supervisor['outcome'] == 'completed' and prior_supervisor['returncode'] == 0, 'original execution')
    for name, digest in prior_supervisor['output_sha256'].items():
        require(sha(old/name) == digest, 'original completed output seal')
    require(request['candidates'] == [row for row in read(old/'incidence.json')['cells']
                                    if row['status'] == 'SPLIT_REQUIRES_ADMISSION'], 'sealed split roster')
    require(len(request['candidates']) == len(admission['candidates']) == 1, 'single split')
    require(admission['candidates'][0]['status'] == 'UNKNOWN_FINITE_COLUMN_IN_SPAN', 'original UNKNOWN status')

    R = PolynomialRing(QQ, 't'); K = R.fraction_field(); t = QQ(parameter); d = t.denominator()
    A, B = [R(parent[k]) for k in ('A_coefficients_low_to_high', 'B_coefficients_low_to_high')]

    def rational(row):
        return K(R(row['numerator_coefficients_low_to_high']))/R(row['denominator_coefficients_low_to_high'])

    sections = [(rational(row['X']), rational(row['Y'])) for row in parent['sections']]
    require(len(sections) == parent['generic_rank_lower_bound'] == 17, 'generic-only parent')
    require([row['basis_index'] for row in parent['sections']] == list(range(17)), 'generic ordering')
    require(all(y*y == x**3+A*x+B for x, y in sections), 'generic section equations')
    curve = ['0', '0', '0', str(A(t)*d**8), str(B(t)*d**12)]
    basis = [[str(x(t)*d**4), str(y(t)*d**6)] for x, y in sections]
    require(curve == seed['curve'] and basis == seed['points'], 'independent generic specialization')

    trace_folder = inputs/'traces'/str(mask)
    trace, trace_request = map(read, [trace_folder/'result.json', trace_folder/'request.json'])
    trace_supervisor = read(trace_folder/'supervision-receipt.json')
    require(trace_supervisor['outcome'] == 'completed' and trace_supervisor['returncode'] == 0, 'cold trace completion')
    require(trace_supervisor['output_sha256']['result.json'] == sha(trace_folder/'result.json'), 'cold trace seal')
    require(trace['parent_sha256'] == trace_request['parent_sha256'] == sha(inputs/'parent.json')
            and trace['request_sha256'] == sha(trace_folder/'request.json'), 'cold generic input binding')
    require(trace['mask'] == trace_request['mask'] == mask
            and trace['word'] == trace_request['word'], 'selected trace word')
    require(dict(mask=mask, word=trace['word']) in plan['orbits'], 'original orbit selection')
    conic = trace['conic']
    h, nx, ny, m, g, b, k, q = [R(conic[key]) for key in ('h', 'nx', 'ny', 'm', 'g', 'b', 'k', 'q')]
    require(h.degree() == 3 and h.is_monic() and nx.degree() <= 10 and ny.degree() <= 15, 'norm-ten pole chart')
    require(nx.gcd(h*h).degree() == 0 and m.degree() < 6, 'unique Euclidean residue')
    generic_trace = word_point(K(A), sections, [-n for n in trace['word']])
    require(generic_trace == (K(nx)/(h*h), K(ny)/(h**3)), 'independent generic trace group law')
    require(h*h*g == m*nx+ny and h*h*b == m*m-nx and h*h*k == m*b-2*g
            and h*h*q == 4*m*k-3*b*b-4*A, 'Euclidean identities')
    require(q.degree() in (1, 2) and (q.degree() == 1 or q[1]**2 != 4*q[0]*q[2]), 'smooth conic')
    maps = list(map(R, conic['maps']))
    require(maps == [b/2, h/2, -h*k/2, -m/2], 'polynomial lift')
    x0, x1, y0, y1 = maps
    require(y0*y0+y1*y1*q == x0**3+3*x0*x1*x1*q+A*x0+B
            and 2*y0*y1 == 3*x0*x0*x1+x1**3*q+A*x1, 'complete lift identity')
    cell = request['candidates'][0]; square = cell['square_certificate']
    w = QQ(square['root']); value = q(t)
    require(cell['mask'] == mask and cell['parameter'] == parameter
            and square['status'] == 'SQUARE' and w > 0 and w*w == value, 'nonzero exact square split')
    require(str(value) == square['value'] and str(value.numerator()) == square['numerator']
            and str(value.denominator()) == square['denominator']
            and QQ(square['numerator_floor_sqrt'])**2 == value.numerator()
            and QQ(square['denominator_floor_sqrt'])**2 == value.denominator(), 'square receipt')
    split = [x0(t)+x1(t)*w, y0(t)+y1(t)*w]
    require(list(map(str, split)) == cell['point'], 'original conic branch')
    candidate = list(map(str, [split[0]*d**4, split[1]*d**6]))
    require(candidate == admission['candidates'][0]['point'], 'archived short-model candidate')

    frame, receipt = read(folder/'frame.json'), read(folder/'receipt.json')
    require(frame['curve'] == curve and frame['basis'] == basis, 'generic frame equation/prefix')
    signatures = seed['proof']['signatures']
    require(frame['finite_rows'] == [row for sig in signatures for row in sig['rows']]
            and frame['no_two_torsion_prime'] == seed['proof']['no_rational_2_torsion_prime'], 'original finite footprint')
    require([(r['prime'], r['order'], r['dimension']) for r in frame['records']] ==
            [(s['prime'], s['group_order'], s['quotient_dimension']) for s in signatures], 'same sealed places and order')
    require(all(s['doubled_subgroup_order']*2**s['quotient_dimension'] == s['group_order'] for s in signatures),
            'saved doubled subgroup orders')
    require(receipt['original'] == candidate and receipt['index'] == index
            and receipt['parameter'] == parameter and receipt['mask'] == mask, 'case identity')
    require(receipt['protocol_sha256'] == sha(OUT/'protocol.json')
            and receipt['frame_sha256'] == sha(folder/'frame.json')
            and receipt['seed_sha256'] == sha(old/'seed.json'), 'new receipt bindings')
    require(receipt['max_successful_halves'] == protocol['max_successful_halves'], 'case proof bound')
    for name, digest in receipt['source_hashes'].items():
        require(sha(ROOT/name) == digest, 'producer source seal')
    replay = ReplayFrame(frame)
    report = replay.replay(receipt)
    # The other split branch also lies in the same rational span when this
    # one does: check the trace identity directly on this exact fibre.
    other = tuple(F(str(v)) for v in [(x0(t)-x1(t)*w)*d**4, (y0(t)-y1(t)*w)*d**6])
    require(add(replay.A, replay.checked_point(candidate), other) == replay.word(trace['word']), 'split trace identity')
    report.update(schema='independent-euclidean-admission-replay.v2', index=index, parameter=parameter,
        mask=mask, receipt_sha256=sha(folder/'receipt.json'), frame_sha256=sha(folder/'frame.json'),
        protocol_sha256=sha(OUT/'protocol.json'), sage_version=SAGE_VERSION,
        checker_sha256=sha(Path(__file__)), arithmetic_checker_sha256=sha(Path(__file__).with_name('replay_split_admission.py')),
        checks=['source seals', 'generic equations and specialization', 'generic trace word',
                'Euclidean identities and full lift', 'exact square and archived branch',
                'complete sealed finite groups', 'unique parity and subtraction at every layer',
                'manual rational doubling', 'terminal classification', 'companion trace identity'],
        new_parameters=0, new_orbits=0, new_finite_places=0, later_V3_point_inputs=0)
    if receipt['status'] == 'INHERITED_RATIONAL_SPAN':
        report.update(relation_multiplier=receipt['relation_multiplier'], relation_word=receipt['relation_word'])
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', required=True, type=int, choices=[r[0] for r in ROSTER])
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = verify(args.index)
    output = OUT/'cases'/str(args.index)/'independent-replay.json'
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.write and not output.exists():
        with output.open('x') as stream:
            stream.write(text)
    else:
        require(output.read_text() == text, 'immutable independent replay differs')
    print('PASS_INDEPENDENT_EUCLIDEAN_ADMISSION', args.index, result['status'], flush=True)

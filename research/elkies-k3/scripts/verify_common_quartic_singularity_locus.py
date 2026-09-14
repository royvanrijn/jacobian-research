#!/usr/bin/env python3
"""Portable replay: integer resultants, Fraction identities, no Sage/producer.

The general delta, height, and specialization arguments are written proofs.
This checks their exact arithmetic hypotheses, not a full quartic-chart exclusion.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1'
SOURCE = 'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
NAMES = ['published-r17', 'alternate-q80', 'curve302-parent', 'x1092-class1']
SPEC = importlib.util.spec_from_file_location(
    'common_quartic_arithmetic', ROOT / 'elkies-k3/scripts/verify_r17_mestre_shared_twist.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
V, require = M.V, M.require


def fscale(a, c, p):
    return M.trim([x*c % p for x in a])


def fadd(a, b, p):
    return M.ff_sub(a, fscale(b, -1, p), p)


def fpow(a, n, p):
    out = [1]
    for _ in range(n):
        out = M.ff_mul(out, a, p)
    return out


def fderivative(a, p):
    return M.trim([i*c % p for i, c in enumerate(a)][1:])


def resultant(a, b, p):
    """Euclidean resultant, retaining leading coefficients and swap signs."""
    a, b, answer = M.trim(a), M.trim(b), 1
    if not a or not b:
        return 0
    while len(b) > 1:
        m, n = len(a)-1, len(b)-1
        r = M.ff_rem(a, b, p)
        if not r:
            return 0
        answer = answer * (-1 if m*n % 2 else 1) * pow(b[-1], m-len(r)+1, p) % p
        a, b = b, r
    return answer * pow(b[0], len(a)-1, p) % p


def interpolate_consecutive(values, p):
    """Newton interpolation at 0,...,n-1, with n<p."""
    require(len(values) < p, 'interpolation collisions')
    ds = list(values)
    for step in range(1, len(ds)):
        inv = pow(step, -1, p)
        for i in range(len(ds)-1, step-1, -1):
            ds[i] = (ds[i]-ds[i-1])*inv % p
    answer, basis = [], [1]
    for i, c in enumerate(ds):
        answer = fadd(answer, fscale(basis, c, p), p)
        basis = M.ff_mul(basis, [(-i) % p, 1], p)
    return answer


def verify_parent(parent, row):
    require(row['name'] == parent['name'], 'parent attachment')
    p = row['prime']
    require(p == 1009, 'fixed witness prime')
    require([len(parent[k]) for k in ('A', 'B')] == [9, 13], 'parent degrees')
    A, B = [M.reduction(V.poly(parent[k]), p) for k in ('A', 'B')]
    Ap, Bp = fderivative(A, p), fderivative(B, p)
    H = M.ff_sub(M.ff_sub(M.ff_mul(B, fpow(Ap, 3, p), p),
                         M.ff_mul(M.ff_mul(A, Bp, p), fpow(Ap, 2, p), p), p),
                 fpow(Bp, 3, p), p)
    require(len(H) == 34, 'critical polynomial degree drop')
    require(H == row['critical_t_polynomial_mod_p'], 'critical polynomial identity')
    require(len(M.ff_gcd(H, Ap, p)) == 1, 'critical derivative denominator')
    require(len(M.ff_gcd(H, fderivative(H, p), p)) == 1, 'multiple critical t')
    # Res_t(H, c*A'+B') has degree <=33 in c over Q before reduction.
    # Thus 34 exact residues determine its entire reduction, not a parameter sample.
    values = [resultant(H, fadd(fscale(Ap, c, p), Bp, p), p) for c in range(34)]
    N = interpolate_consecutive(values, p)
    require(len(N) == 34, 'critical value norm degree drop')
    N = fscale(N, pow(N[-1], -1, p), p)
    require(N == row['critical_value_norm_monic_mod_p'], 'critical value norm identity')
    require(len(M.ff_gcd(N, fderivative(N, p), p)) == 1, 'critical values collide')
    require(row['constant_x_genus_lower_bound'] == 4 and
            row['critical_values_over_Qbar'] == 33, 'genus-gate conclusion')
    return {'name': parent['name'], 'prime': p, 'critical_values': 33,
            'all_algebraic_constant_x_genus_at_least': 4}


def genus_from_multiplicities(finite, degree=12):
    """A geometric multiplicity roster, with infinity restored as degree-deg(F)."""
    require(all(isinstance(m, int) and m > 0 for m in finite), 'invalid multiplicities')
    infinity = degree-sum(finite)
    require(infinity >= 0, 'homogeneous degree overflow')
    roster = list(finite) + ([infinity] if infinity else [])
    branch = sum(m % 2 for m in roster)
    require(branch > 0 and branch % 2 == 0, 'split or constant quadratic field')
    return {'infinity_multiplicity': infinity,
            'delta': sum(m//2 for m in roster), 'branch_degree': branch,
            'genus': (branch-2)//2}


def smooth(f, p, message):
    a = M.reduction(f, p)
    require(len(M.ff_gcd(a, fderivative(a, p), p)) == 1, message)


def coprime(a, b, p, message):
    require(len(M.ff_gcd(M.reduction(a, p), M.reduction(b, p), p)) == 1, message)


def elliptic_add(P, R, A):
    if P is None:
        return R
    if R is None:
        return P
    x, y = P; u, v = R
    if x == u and y == -v:
        return None
    m = (3*x*x+A)/(2*y) if P == R else (v-y)/(u-x)
    z = m*m-x-u
    return z, m*(x-z)-y


def verify_base_certificate(row):
    cert = row['rational_base_certificate']
    t = V.poly([0, 1]); one = V.poly([1])
    if cert['type'] == 'rational':
        require(cert == {'type': 'rational', 't': '2*z/(1-z^2)',
                         'w': '(1+z^2)/(1-z^2)'}, 'rational map attachment')
        den, tn, wn = V.sub(one, V.power(t, 2)), V.scale(t, 2), V.add(one, V.power(t, 2))
        require(V.power(wn, 2) == V.add(V.power(tn, 2), V.power(den, 2)), 'conic parametrization')
        return 'rational parametrization'
    require(cert['type'] == 'elliptic' and cert['a_invariants'] == [0, 0, 0, -4, 1], 'pointed quartic model')
    require(cert['forward'] == ['X=2*(w+t^2)', 'Y=2*t*X+1'] and
            cert['inverse'] == ['t=(Y-1)/(2*X)', 'w=X/2-t^2'], 'quartic birational map attachment')
    # Work in Q[t,w]/(w^2-d). Store an element as (coefficient of 1, of w).
    d = V.poly(row['d'])
    def qa(a, b):
        return V.add(a[0], b[0]), V.add(a[1], b[1])
    def qs(a, b):
        return V.sub(a[0], b[0]), V.sub(a[1], b[1])
    def qm(a, b):
        return (V.add(V.mul(a[0], b[0]), V.mul(d, V.mul(a[1], b[1]))),
                V.add(V.mul(a[0], b[1]), V.mul(a[1], b[0])))
    X = V.scale(V.power(t, 2), 2), V.poly([2])
    Y = qa(qm((V.scale(t, 2), []), X), (one, []))
    require(qs(qm(Y, Y), qa(qs(qm(qm(X, X), X), qm((V.poly([4]), []), X)), (one, []))) == ([], []), 'quartic to elliptic identity')
    # Clearing the inverse denominators returns t and w identically.
    require(qs(Y, (one, [])) == qm((V.scale(t, 2), []), X), 'inverse t identity')
    require(qs((V.scale(X[0], Q(1, 2)), V.scale(X[1], Q(1, 2))), (V.power(t, 2), [])) == ([], one), 'inverse w identity')
    P = tuple(Q(v) for v in cert['point'])
    require(P == (Q(0), Q(1)) and P[1]**2 == P[0]**3-4*P[0]+1, 'base point equation')
    P2, P3 = elliptic_add(P, P, Q(-4)), None
    P3 = elliptic_add(P2, P, Q(-4))
    require(list(map(str, P2)) == cert['twice_point'] and
            list(map(str, P3)) == cert['three_times_point'], 'base point multiples')
    require(P3[0].denominator > 1, 'missing Lutz-Nagell nontorsion witness')
    return 'elliptic model with certified rational nontorsion point'


def verify_control(input_row, row):
    require(row['name'] == input_row['name'], 'control attachment')
    for k in ('D', 'square_factor', 'd', 's'):
        require(V.poly(row[k]) == V.poly(input_row[k]), 'control input projection')
    D, e, d, s, A, B = [V.poly(row[k]) for k in ('D', 'square_factor', 'd', 's', 'A', 'B')]
    one, p = V.poly([1]), row['smoothness_prime']
    require(p == 1009, 'control witness prime')
    require(D == V.mul(V.power(e, 2), d), 'literal squareclass identity')
    require(A == V.sub(V.mul(D, V.add(V.scale(s, 2), one)), one) and
            B == V.mul(D, V.power(s, 2)), 'reverse construction identity')
    delta = V.add(V.scale(V.power(A, 3), 4), V.scale(V.power(B, 2), 27))
    require([len(A), len(B), len(delta)] == [9, 13, 25], 'K3 and infinity degrees')
    for f in (delta, d, s, V.add(s, one)):
        smooth(f, p, 'nonsimple fibre, branch, or node')
    for f in (delta, s, V.add(s, one)):
        coprime(D, f, p, 'branch/node collision')
    expected = [([], V.mul(e, s)), (one, V.mul(e, V.add(s, one))),
                (V.sub(D, one), V.mul(e, V.sub(V.sub(one, D), s)))]
    require(len(row['points_and_sum']) == 3, 'two points and sum coverage')
    for actual, (x, r) in zip(row['points_and_sum'], expected):
        require((V.poly(actual['x']), V.poly(actual['y_over_w'])) == (x, r), 'point or sum attachment')
        require(V.add(V.add(V.power(x, 3), V.mul(A, x)), B) == V.mul(d, V.power(r, 2)), 'section equation')
        require(len(x)-1 <= 4 and 2*(len(r)-1)+len(d)-1 <= 12, 'section infinity pole')
    # The first two x coordinates differ by 1. Their chord slope is e*w.
    slope = V.sub(expected[1][1], expected[0][1])
    sum_x = V.sub(V.mul(d, V.power(slope, 2)), one)
    sum_y = V.sub(V.mul(slope, V.sub(expected[0][0], sum_x)), expected[0][1])
    require((sum_x, sum_y) == expected[2], 'sum group law')
    # Smooth branch and only I1 fibres; all three displayed sections miss O.
    heights = [8, 8, 8]
    cross = (heights[2]-heights[0]-heights[1])//2
    G = [[8, cross], [cross, 8]]
    require(row['height_matrix'] == G and row['height_determinant'] == 64-cross**2 == 48,
            'height determinant')
    g = (len(d)-3)//2
    require(row['normalization_genus'] == g and row['branch_degree'] == len(d)-1 and
            row['image_arithmetic_genus'] == 5 and row['nodes_per_image'] == 5-g, 'control genus')
    profile = [2]*4 + ([1]*4 if g == 1 else [2, 1, 1])
    require(genus_from_multiplicities(profile)['genus'] == g, 'node multiplicity profile')
    if g == 0:
        require(e == V.poly([0, 1]) and d == V.poly([1, 0, 1]), 'double-root control')
    else:
        require(e == one and d == V.poly([1, 1, 0, 0, 1]), 'elliptic control')
    # j is nonconstant. Evaluations are a certificate, not a parameter search.
    def ev(f, x):
        return sum(c*x**i for i, c in enumerate(f))
    require(ev(A, 0)**3*ev(B, 1)**2 != ev(A, 1)**3*ev(B, 0)**2, 'isotrivial parent')
    require(row['inherited_parent_rank'] == 'UNKNOWN' and not row['is_a_retained_mw17_parent'], 'parent rank scope')
    require(row['new_directions_modulo_full_inherited_group'] == 2, 'character independence conclusion')
    return {'name': row['name'], 'genus': g, 'nodes_per_image': 5-g,
            'height_matrix': G, 'new_independent_directions': 2,
            'infinite_rational_base': verify_base_certificate(row),
            'inherited_parent_rank': 'UNKNOWN'}


def verify(path):
    packet = json.loads((path / 'input.json').read_text())
    result = json.loads((path / 'result.json').read_text())
    require(packet['schema'] == 'common-quartic-singularities-input-v1' and
            result['schema'] == 'common-quartic-singularities-result-v1', 'certificate schema')
    require(result['input_sha256'] == sha256((path / 'input.json').read_bytes()).hexdigest(), 'input binding')
    require(packet['source'] == SOURCE and packet['source_sha256'] == sha256((ROOT / SOURCE).read_bytes()).hexdigest(), 'generic source binding')
    for relative, digest in packet['arithmetic_helpers'].items():
        require(sha256((ROOT / relative).read_bytes()).hexdigest() == digest, 'arithmetic helper binding')
    parents = json.loads((ROOT / SOURCE).read_text())['parents']
    require(packet['parents'] == NAMES == [p['name'] for p in parents] ==
            [p['name'] for p in result['parents']], 'complete parent coverage')
    require([r['name'] for r in packet['controls']] == [r['name'] for r in result['controls']] ==
            ['four-nodes-genus-one', 'five-nodes-genus-zero'], 'complete control coverage')
    require(result['fixed_parent_full_quartic_chart'] == 'UNKNOWN' and
            result['mw17_endpoint_complete'] is False, 'fixed-parent endpoint scope')
    gates = [verify_parent(a, b) for a, b in zip(parents, result['parents'])]
    controls = [verify_control(a, b) for a, b in zip(packet['controls'], result['controls'])]
    for row in result['parents'] + result['controls']:
        require(json.loads((path / (row['name'] + '-checkpoint.json')).read_text()) == row, 'checkpoint equality')
    return {'status': 'PASS_COMMON_QUARTIC_CONTROLS_AND_CONSTANT_X_EXCLUSIONS',
            'parents': gates, 'controls': controls,
            'fixed_parent_full_quartic_chart': 'UNKNOWN', 'mw17_endpoint_complete': False}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', type=Path, default=DEFAULT)
    p.add_argument('--output', type=Path)
    args = p.parse_args(); start = time.monotonic()
    answer = verify(args.input_dir)
    answer['elapsed_seconds'] = round(time.monotonic()-start, 6)
    answer['input_hashes'] = {name: sha256((args.input_dir / name).read_bytes()).hexdigest()
                            for name in ('input.json', 'result.json')}
    if args.output:
        with args.output.open('x') as stream:
            json.dump(answer, stream, indent=2, sort_keys=True); stream.write('\n')
    print(json.dumps(answer, sort_keys=True))


if __name__ == '__main__':
    main()

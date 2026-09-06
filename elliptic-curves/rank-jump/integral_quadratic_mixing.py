#!/usr/bin/env python3
"""Integral Galois mixing versus a rational/Sha switch: two small controls."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import local_collision as lc
import quadratic_norm_blocks as norms

PROTOCOL = Path(__file__).with_name('INTEGRAL_QUADRATIC_MIXING_PROTOCOL.json')
SMALL = r.OUT/'rank_jump_small_quotient_block_v1.json'
CONTRACTION = r.OUT/'rank_jump_generic_contraction_consequences_v1.json'
OUTPUT = r.OUT/'rank_jump_integral_quadratic_mixing_v1.json'
WORK = r.ROOT/'artifacts/local/rank-jump-integral-quadratic-mixing-v1'


def span(basis):
    out = {0}
    for v in basis:
        out |= {x ^ v for x in list(out)}
    return out


def module_counts(rplus, rminus, h):
    assert 0 <= h <= min(rplus, rminus)
    return {'trivial_summands': rplus-h, 'sign_summands': rminus-h,
            'regular_summands': h, 'integral_eigenspace_index': 2**h,
            'norm_index_plus': 2**(rplus-h), 'norm_index_minus': 2**(rminus-h)}


def compute():
    from sage.all import QQ, GF, PolynomialRing, NumberField, EllipticCurve, pari
    from sage.version import version
    pol = PolynomialRing(QQ, 'z'); z = pol.gen()
    F = NumberField(z*z+1, 'ii'); ii = F.gen()
    E = EllipticCurve(QQ, [0, 0, 0, 5, -3])
    Et = EllipticCurve(QQ, [0, 0, 0, 5, 3])
    EF = E.change_ring(F)
    T, Tc = EF(ii, 1+2*ii), EF(-ii, 1-2*ii)
    P, Q = E(4, -9), Et(1, -3)
    assert T+Tc == EF(P) and T-Tc == EF(-Q[0], ii*Q[1])
    assert 2*T == EF(P)+EF(-Q[0], ii*Q[1])
    cubic = z**3+5*z-3
    assert cubic.is_irreducible() and cubic.discriminant() == -743
    assert ((4-z)*(1+z)-(z*z+2)**2) % cubic == 0
    descents, reductions = [], []
    for curve, point in [(E, P), (Et, Q)]:
        out = pari.ellrank(pari.ellinit(curve.a_invariants()), 0, [list(point)[:2]])
        descents.append({'ainvs': list(map(str, curve.a_invariants())),
                         'raw': str(out), 'lower_bound': int(out[0]),
                         'upper_bound': int(out[1]), 'Sha2_mod_2Sha4_dimension': int(out[2])})
        red = curve.change_ring(GF(3)); pr = red(point)
        doubles = {2*q for q in red}
        assert pr not in doubles and red.order() == 4
        reductions.append({'prime': 3, 'point': list(map(str, point[:2])),
                           'group_order': 4, 'point_order': int(pr.order()),
                           'point_in_double_image': False})
    # A nonzero Kummer class is common; exact ranks bound the intersection above.
    exact_ranks = [x['lower_bound'] if x['lower_bound'] == x['upper_bound'] else None for x in descents]
    construct = {'curve': [0, 0, 0, 5, -3], 'twist': [0, 0, 0, 5, 3],
                 'prescribed_quadratic_point': ['i', '1+2*i'],
                 'trace': ['4', '-9'], 'rational_anti_trace': ['1', '-3'],
                 'shared_Kummer_identity': '(4-theta)*(1+theta)=(theta^2+2)^2',
                 'cubic_discriminant': -743, 'descents': descents,
                 'reduction_nondivisibility': reductions, 'exact_ranks': exact_ranks,
                 'common_rational_Kummer_dimension_at_least': 1,
                 'common_rational_Kummer_dimension_exact': 1 if min(x['upper_bound'] for x in descents) == 1 else 'UNKNOWN'}
    if exact_ranks == [1, 1]:
        construct['module'] = module_counts(1, 1, 1)
        construct['local_global_norm_defect_dimensions'] = [0, 0]

    old = r.read(SMALL)
    assert old['exact_MW_ranks'] == [1, 3] and old['full_2_Selmer_dimensions'] == [3, 3]
    # Coordinates (beta0,beta1,u0,u1), with theta=u0*u1.
    G0, G1 = span([4]), span([1, 2, 12])
    S0, S1 = span([1, 2, 4]), span([1, 2, 12])
    assert G0 & G1 == {0} and S0 & S1 == span([1, 2])
    assert G0 & S1 == {0} and G1 & S0 == span([1, 2])
    closed = {'coordinate_order': ['beta0', 'beta1', 'u0', 'u1'],
              'standard_scalar_labels': True, 'G0_basis': [4], 'G1_basis': [1, 2, 12],
              'S0_basis': [1, 2, 4], 'S1_basis': [1, 2, 12],
              'common_rational_Kummer_dimension': 0, 'common_Selmer_dimension': 2,
              'module': module_counts(1, 3, 0),
              'local_norm_quotient_dimensions': [1, 1],
              'local_global_norm_defect_dimensions': [0, 2],
              'boundary': 'Complete rational and Selmer spaces inherited from the small-control proof. No order-three genus-two label transport.'}
    cs = {x['case_index']: x for x in r.read(CONTRACTION)['rows']}
    production = []
    for row in r.read(norms.OUTPUT)['rows']:
        i = row['case_index']; d = cs[i]['exact_Selmer_dimension_drop']; k = row['CT_rank']
        production.append({'case_index': i, 'id': row['id'],
                           'generic_saturated_local_norm_codimension': d,
                           'certified_norm_defect_dimension_at_least': k,
                           'trivial_integral_summands_at_least': d+k,
                           'regular_summand_upper_bound': f'rank(E/Q)-{d+k}',
                           'with_rank_equal_to_known_witness': row['known_independent_rank']-d-k,
                           'last_column_condition': 'Conditional on no undiscovered free directions; not an exact-rank assertion.',
                           'regular_summands_exact': 'UNKNOWN'})
    paths = (Path(__file__), PROTOCOL, SMALL, CONTRACTION, norms.OUTPUT)
    return {'schema': 'rank-jump.integral-quadratic-mixing.v1',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths},
            'software': {'sage': version, 'pari': str(pari('version()'))},
            'status': 'PASS', 'constructive_control': construct, 'closed_control': closed,
            'production': production,
            'boundary': 'Integral 2-primary mixing, not free-rank prediction. The only new descents are on the two prescribed small curves.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['capture', 'worker', 'check']); args = p.parse_args()
    if args.mode == 'worker':
        r.write_new(WORK/'checkpoint.json', compute())
    elif args.mode == 'check':
        assert r.read(OUTPUT) == compute(); print('PASS integral quadratic mixing')
    else:
        WORK.mkdir(parents=True, exist_ok=True); path = WORK/'checkpoint.json'
        if not path.exists():
            with (WORK/'worker.log').open('x') as log:
                try:
                    proc = subprocess.run(['sage', '-python', str(Path(__file__).resolve()), 'worker'],
                                          cwd=r.ROOT, stdout=log, stderr=log, timeout=15)
                    reason = None if proc.returncode == 0 else 'worker failure'
                except subprocess.TimeoutExpired:
                    reason = '15-second timeout'
                if reason and not path.exists():
                    r.write_new(path, {'status': 'UNKNOWN', 'reason': reason})
        result = r.read(path); r.write_new(OUTPUT, result); print(result['status'])

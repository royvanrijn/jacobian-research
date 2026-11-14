"""Exact Puiseux-packet bounds and normal/root inclusion controls, PF2OMIT10."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import platform
import time
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []
KINDS = {'unramified': (1, 1), 'fold': (2, 1), 'cyclic3': (3, 1), 'cusp3': (3, 2)}
PACKETS = {'A0': ['unramified']*5,
           'A1': ['fold']+['unramified']*3,
           'A2': ['fold', 'fold', 'unramified'],
           'B': ['cusp3', 'fold'], 'C': ['cyclic3', 'fold'],
           'Wright_quartic': ['cusp3']}


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(kinds, poles):
    weights = sorted(q for kind, pole in zip(kinds, poles)
                     for q in [F(KINDS[kind][1]*pole, KINDS[kind][0])]*KINDS[kind][0])
    nu = sum(KINDS[kind][0]-1 for kind in kinds)
    delta = 2*sum(min(x, y) for i, x in enumerate(weights) for y in weights[i+1:])
    coeffs = [math.ceil(sum(weights[:i], F(0))) for i in range(1, len(weights)+1)]
    return weights, coeffs, nu, delta, math.ceil((delta-nu)/2)


def order(poly, epsilon):
    if s.expand(poly) == 0:
        return math.inf
    return min(m[0] for m, c in s.Poly(poly, epsilon).terms())


def control(label, kinds, poles, unit, epsilon, w, z):
    parts, rows, normal_disc = [], [], s.Integer(1)
    total = sum(KINDS[k][0] for k in kinds)
    for j, (kind, pole) in enumerate(zip(kinds, poles)):
        k, r = KINDS[kind]
        eq = w**k-epsilon
        coord = (j+1)*w**(r*pole)*(1+w if unit else 1)
        parts.append(s.resultant(eq, z-coord, w))
        normal_disc *= s.discriminant(eq, w)
        reduced = [s.Poly(s.rem(s.expand(coord**i), eq, w), w) for i in range(total)]
        rows.extend([[poly.nth(l) for poly in reduced] for l in range(k)])
    G = s.Poly(s.prod(parts), z).as_expr().expand()
    determinant = s.Matrix(rows).det(method='domain-ge').expand()
    raw_disc = s.discriminant(G, z)
    weights, coeffs, nu, delta, idx = bounds(kinds, poles)
    check(label+' separable full normal comparison', determinant != 0 and raw_disc != 0)
    check(label+' full trace discriminant/index-square identity',
          s.expand(raw_disc-normal_disc*determinant**2) == 0)
    check(label+' normalized different order', order(normal_disc, epsilon) == nu)
    check(label+' raw root-pair order bound', order(raw_disc, epsilon) >= delta)
    check(label+' normalization index bound', order(determinant, epsilon) >= idx)
    # R_(n+1)(T)=T^(n+1)*G(1/T); the retained root is T=0.
    G_poly = s.Poly(G, z)
    actual_coeffs = [order(G_poly.nth(total-i), epsilon) for i in range(1, total+1)]
    check(label+' all inverse coefficient bounds', all(a >= b for a, b in zip(actual_coeffs, coeffs)))
    check(label+' exact leading pole charge', actual_coeffs[-1] == sum(weights))
    check(label+' exact discriminant parity', (order(raw_disc, epsilon)-nu) % 2 == 0)
    if not unit and all(p == 1 for p in poles):
        check(label+' minimal-pole raw bound attained', order(raw_disc, epsilon) == delta)
        check(label+' minimal-pole index bound attained', order(determinant, epsilon) == idx)
    return dict(label=label, poles=poles, extra_unit_1_plus_w=unit,
                reciprocal_polynomial=str(G), normal_trace_discriminant=str(s.expand(normal_disc)),
                normalization_inclusion_determinant=str(determinant),
                raw_discriminant_order=order(raw_disc, epsilon),
                index_order=order(determinant, epsilon),
                inverse_coefficient_orders=[None if x == math.inf else x for x in actual_coeffs],
                zero_coefficient_order_encoding='null means coefficient identically zero, hence infinite order',
                is_global_Keller_realization=False)


def main():
    started = time.monotonic()
    epsilon, w, z = s.symbols('epsilon w z')
    prior = json.loads((HERE/'simple_radial_ledger.json').read_text())
    for row in prior['singleton_radial_ledger']:
        check(row['type']+' agrees with imported complete original packet',
              PACKETS[row['type']] == row['omitted_blocks'])
    expected = {'A0': ([1, 2, 3, 4, 5], 0, 20, 10),
                'A1': ([1, 1, 2, 3, 4], 1, 13, 6),
                'A2': ([1, 1, 2, 2, 3], 2, 10, 4),
                'B': ([1, 1, 2, 3, 3], 3, 11, 4),
                'C': ([1, 1, 1, 2, 2], 3, 7, 2),
                'Wright_quartic': ([1, 2, 2], 2, 4, 1)}
    rows, controls = [], []
    for name, kinds in PACKETS.items():
        weights, coeffs, nu, delta, idx = bounds(kinds, [1]*len(kinds))
        check(name+' exact minimum-pole table', (coeffs, nu, delta, idx) == expected[name])
        check(name+' complete singleton length', 1+sum(KINDS[k][0] for k in kinds) == (4 if name == 'Wright_quartic' else 6))
        check(name+' integral pair-valuation bound', delta.denominator == 1)
        rows.append(dict(packet=name, blocks=kinds,
                         minimum_pole_reciprocal_valuations=[str(q) for q in weights],
                         coefficient_order_lower_bounds=coeffs,
                         normalized_discriminant_order=nu,
                         raw_discriminant_order_lower_bound=int(delta), index_order_lower_bound=idx,
                         actual_original_prime_pole_orders=None,
                         actual_affine_images=None, actual_singleton_target_value=None))
        controls.append(control(name+'_p1', kinds, [1]*len(kinds), False, epsilon, w, z))
        for pole in [2, 3]:
            controls.append(control(name+'_p'+str(pole), kinds, [pole]*len(kinds), True, epsilon, w, z))
    for name, poles in [('B', [1, 2]), ('B', [2, 1]), ('C', [1, 2]), ('C', [2, 1])]:
        controls.append(control(name+'_mixed_'+''.join(map(str, poles)), PACKETS[name], poles, True, epsilon, w, z))
    # Check both weight-order regimes, without assuming pole orders are equal.
    for p, r in [(1, 1), (1, 2), (2, 1), (3, 2), (2, 3)]:
        b = bounds(PACKETS['B'], [p, r])
        c = bounds(PACKETS['C'], [p, r])
        check(f'cusp-fold exact general formula {p},{r}',
              sum(b[0]) == 2*p+r and b[3] == 4*p+r+12*min(F(2*p, 3), F(r, 2)))
        check(f'cyclic-fold exact general formula {p},{r}',
              sum(c[0]) == p+r and c[3] == 2*p+r+12*min(F(p, 3), F(r, 2)))
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT10',
        hypotheses='PF2OMIT7 complete residual r7 singleton packets with admissible primitive source T; marked quartic refinement imports PF2OMIT9',
        reciprocal_unit_times_boundary_power=True,
        reciprocal_block_valuations={k: dict(ambient_length=n, reciprocal_valuation=f'{r}*p/{n}',
                                            original_generic_e=prior['local_blocks'][k]['original_generic_e'],
                                            original_residue_f=1, omitted=True,
                                            ramifies=prior['local_blocks'][k]['original_generic_e'] > 1,
                                            actual_original_affine_image_equation=None)
                                     for k, (n, r) in KINDS.items()},
        packet_rows=rows, local_one_parameter_controls=controls,
        minimum_raw_and_index_bounds_attained_by_controls_only=True,
        global_prime_pole_identification_rule='Use one common p_D for all points of the same original omitted prime; do not assign them independent orders',
        actual_global_inverse=None, unrestricted_singleton_zero_properness_proved=False,
        full_goal_complete=False)
    target = HERE/'singleton_valuation_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/n for n in ['SINGLETON_INVERSE_VALUATIONS.md', 'ADMISSIBLE_INVERSE.md',
              'DISCRIMINANT_LEDGER.md', 'SIMPLE_RAMIFICATION_AND_RADIAL.md',
              'simple_radial_ledger.json', 'WRIGHT_IMAGE_RESTRICTION.md']]
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target), runtime_seconds=round(time.monotonic()-started, 6),
        python=platform.python_version(), sympy=s.__version__,
        control_count=len(controls),
        proof_boundary='Exact whole finite-algebra inclusion determinants and trace discriminants for local one-parameter controls. Universal pole, generic-line and root-product arguments are written proofs; controls are not global Keller realizations.',
        formal_verification=False, external_review=False, full_goal_complete=False)
    (HERE/'singleton_valuation_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'control_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()

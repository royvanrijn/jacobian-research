#!/usr/bin/env sage-python
"""Bounded, retrospective first-witness Jacobian Kummer diagnostic.

Freeze the generic inherited-line atlas before reading the first-witness
pencil. No point search, global class group, pilot input, or control selection.
"""
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, prime_range, power_mod

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves/det1092-rr-residual-jacobian-class-v1'
OUT = ART/'det1092_rr_residual_jacobian_class_v1.json'
ATLAS = ART/'det1092_rr_inherited_line_atlas_v1.json'
PROTOCOL = ART/'det1092_rr_residual_jacobian_protocol_v1.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(path, data):
    payload = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if path.exists():
        assert path.read_text() == payload
    else:
        path.write_text(payload)


def two_torsion_fixed_count(degrees):
    permutation = []
    offset = 0
    for degree in degrees:
        permutation.extend(offset+(i+1) % degree for i in range(degree))
        offset += degree
    assert offset == 6
    representatives = {min(mask, mask ^ 63) for mask in range(64) if mask.bit_count() % 2 == 0}
    fixed = 0
    for mask in representatives:
        moved = sum(((mask >> i) & 1) << permutation[i] for i in range(6))
        fixed += min(moved, moved ^ 63) == mask
    return fixed


def generic_atlas():
    paths = [ART/name for name in ['curve302_recovered_mw17_parent_v1.json',
             'det1092_first_centre_rr_net_v1.json', 'det1092_rr_net_reducible_locus_v2.json',
             'det1092_rr_inherited_point_v1.json']]
    parent, net, red, inherited = [json.loads(p.read_text()) for p in paths]
    R = PolynomialRing(QQ, 'T')
    A, B = [[R(row) for row in net[key]] for key in ['A', 'B']]
    source_parameters = [QQ(0), QQ(1), QQ(-1)]
    def at(row, t):
        return R(row['numerator'])(t)/R(row['denominator'])(t)
    curves = [EllipticCurve(QQ, [at(row, t) for row in parent['a_invariants']]) for t in source_parameters]
    assert all(E.discriminant() for E in curves)
    bases = [[E([at(row, t) for row in point]) for point in parent['basis_weierstrass_coordinates']]
             for E, t in zip(curves, source_parameters)]
    lines = []
    for index, pair in enumerate(red['section_pairs']):
        word = min(pair['words'], key=lambda w: (sum(abs(n) for n in w), tuple(w)))
        if not any(word):
            intercept, slope = map(QQ, inherited['point_formula']['unique_base_parameter'])
            values = [intercept+slope*t for t in source_parameters]
        else:
            values = []
            for E, basis, t in zip(curves, bases, source_parameters):
                P = sum((n*Q for n, Q in zip(word, basis)), E(0))
                assert P and P[2]
                av = A[0](t)+A[1](t)*P[0]+A[2](t)*P[1]
                bv = B[0](t)+B[1](t)*P[0]+B[2](t)*P[1]
                assert av != 0
                values.append(-bv/av)
            intercept, slope = values[0], values[1]-values[0]
            assert values[2] == intercept-slope
        lines.append({'pair_index': index, 'generic_word': word,
                      'intercept': str(intercept), 'slope': str(slope),
                      'source_evaluations': list(map(str, values))})
        retain(LOCAL/f'generic-line-{index:02d}.json', lines[-1])
    assert len(lines) == 22
    assert len({(r['intercept'], r['slope']) for r in lines}) == 22
    result = {
        'classification': 'verified application of the certified degree-one section intersection theorem',
        'status': 'PASS_GENERIC_INHERITED_RR_LINE_ATLAS',
        'line_equation': 'r=intercept+slope*T; intersect r=u+v*T',
        'source_parameters': list(map(str, source_parameters)),
        'lines': lines,
        'geometric_input': 'D.P_x=1: each known section maps to a line in the RR double plane; projection from the node is the original T-coordinate.',
        'selection_boundary': 'Only generic sections; no first-witness coordinate, exceptional point, control outcome or pilot input.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'script_sha256': sha(Path(__file__))
    }
    retain(ATLAS, result)
    return result


def audit(atlas, protocol):
    # The retrospective diagnostic begins only after the generic atlas is fixed.
    fp = ART/'det1092_first_witness_pencil_genus_gate_v1.json'
    first = json.loads(fp.read_text())
    R = PolynomialRing(QQ, 'T')
    q = R([QQ(row[0]) for row in first['q_coefficients_t_then_v']])
    c, u0 = QQ(first['branch_squareclass_scale']), QQ(first['u0'])
    assert q.degree() == 6 and q.gcd(q.derivative()) == 1
    assert c*q(0) != 0 and (c*q(0)).is_square()
    xvalues = []
    for line in atlas['lines']:
        intercept, slope = QQ(line['intercept']), QQ(line['slope'])
        assert slope != 0
        x = (u0-intercept)/slope
        assert c*q(x) != 0 and (c*q(x)).is_square()
        xvalues.append(x)
    base_index = next(i for i, line in enumerate(atlas['lines']) if not any(line['generic_word']))
    base_x = xvalues[base_index]
    inherited_indices = [i for i in range(22) if i != base_index]
    columns = [xvalues[i] for i in inherited_indices]+[QQ(0)]
    rows, trials = [], []
    torsion_free_prime = None
    for p in protocol['primes']:
        trial = {'p': p}
        if any(a.denominator() % p == 0 for a in list(q)+[base_x]+columns):
            trial['status'] = 'SKIP_DENOMINATOR'
        else:
            K = GF(p)
            P = PolynomialRing(K, 'X')
            X = P.gen()
            f = P([K(a) for a in q])
            if f.degree() != 6 or f.gcd(f.derivative()) != 1:
                trial['status'] = 'SKIP_BAD_SEXTIC'
            elif any(f(K(x)) == 0 for x in [base_x]+columns):
                trial['status'] = 'SKIP_BRANCH_REDUCTION'
            else:
                assert c.valuation(p) % 2 == 0
                factors = sorted([g.monic() for g, exponent in f.factor()
                                  if exponent == 1], key=lambda g: (g.degree(), list(map(int, g))))
                assert sum(g.degree() for g in factors) == 6
                bits = []
                for g in factors:
                    exponent = (ZZ(p)**g.degree()-1)//2
                    anchor = power_mod(P(K(base_x))-X, exponent, g)
                    assert anchor in [P(1), P(-1)]
                    row = []
                    for x in columns:
                        value = (power_mod(P(K(x))-X, exponent, g)*anchor) % g
                        assert value in [P(1), P(-1)]
                        row.append(int(value == P(-1)))
                    bits.append(row)
                raw = matrix(GF(2), bits)
                assert all(sum(raw[:, j].list(), GF(2)(0)) == 0 for j in range(len(columns)))
                parity = matrix(GF(2), 1, len(factors), [g.degree() % 2 for g in factors])
                quotient = parity.right_kernel().basis_matrix()
                block = quotient*raw
                rows.extend([list(map(int, row)) for row in block.rows()])
                degrees = [int(g.degree()) for g in factors]
                fixed = two_torsion_fixed_count(degrees)
                if fixed == 1 and torsion_free_prime is None:
                    torsion_free_prime = p
                trial.update({'status': 'PASS_LOCAL_KUMMER_BLOCK', 'degrees': degrees,
                              'factors': [list(map(int, g.list())) for g in factors],
                              'raw_character_rows': [list(map(int, row)) for row in raw.rows()],
                              'scalar_quotient_rows': [list(map(int, row)) for row in quotient.rows()],
                              'block_rows': [list(map(int, row)) for row in block.rows()],
                              'geometric_2_torsion_Frobenius_fixed_count': fixed})
        trials.append(trial)
        retain(LOCAL/f'prime-{p:04d}.json', trial)
    M = matrix(GF(2), rows)
    H = M[:, :-1]
    h_rank, total_rank = H.rank(), M.rank()
    separated = total_rank > h_rank
    separator = next((z for z in H.left_kernel().basis() if z*M[:, -1]), None)
    assert (separator is not None) == separated
    result = {
        'classification': 'retrospective verified application; not a blinded control panel',
        'status': 'PASS_RESIDUAL_JACOBIAN_CLASS_SEPARATED' if separated else 'INCONCLUSIVE_NO_SEPARATION_IN_FROZEN_LOCAL_PANEL',
        'curve': {'u': str(u0), 'v': '0', 'q': list(map(str, q.list())), 'scale': str(c)},
        'inherited_x_values': list(map(str, xvalues)), 'base_pair_index': base_index,
        'inherited_column_pair_indices': inherited_indices, 'diagnostic_x': '0',
        'Kummer_element': '(x-theta)/(base_x-theta) modulo algebra squares and rational scalars',
        'trials': trials, 'matrix_rows': rows,
        'inherited_character_rank': int(h_rank), 'total_character_rank': int(total_rank),
        'separating_row_combination': None if separator is None else list(map(int, separator)),
        'rational_Jacobian_2_torsion_zero_prime': torsion_free_prime,
        'independence_boundary': 'A separated class is outside H+2J(Q). A rational rank conclusion additionally needs torsion and inherited-rank control; absence of separation proves no global membership.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in [ATLAS, PROTOCOL, fp, Path(__file__)]}
    }
    retain(OUT, result)
    print(result['status'], 'ranks', h_rank, total_rank, 'J[2] zero prime', torsion_free_prime, flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    LOCAL.mkdir(parents=True, exist_ok=True)
    protocol = {
        'classification': 'frozen retrospective diagnostic protocol',
        'curve_choice': 'Previously fixed first-witness pencil u=u0, with v=0; not target-blind.',
        'generic_atlas': 'Exactly the22 previously certified section pairs, ordered as retained.',
        'primes': list(map(int, list(prime_range(17, 500))[:64])),
        'limits': {'wall_seconds': 25, 'RR_curves': 1, 'inherited_pairs': 22,
                   'primes': 64, 'point_searches': 0, 'global_Selmer_runs': 0,
                   'class_group_runs': 0, 'control_sweeps': 0, 'pilot_changes': 0},
        'failure_semantics': 'Retain every skip and block. No enlargement after a miss. A timeout is incomplete, not an obstruction.',
        'script_sha256': sha(Path(__file__))
    }
    retain(PROTOCOL, protocol)
    atlas = generic_atlas()
    audit(atlas, protocol)

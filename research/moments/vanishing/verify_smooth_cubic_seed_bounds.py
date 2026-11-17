"""Exact auxiliary degree checks for graded seed/cofactor prefix bounds.

Universal filtration and analytic no-constant arguments remain written
mathematics. Models check twisted divisions and every seed column through8.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
from math import comb
import hashlib
import json
import time
import sympy as s
from sympy.polys.rings import ring
from flint import fmpq, fmpq_poly, fmpq_mat

P = Path(__file__).resolve().parent
BASIS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
         (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
GRADES = [sum(m) for m in BASIS]


def model(label, diagonal, coupling, lower):
    R, x, y, z, u = ring('x,y,z,u', s.QQ)
    v = [x, y, z]
    aa = [s.QQ(coupling)/s.QQ(k) for k in diagonal]
    delta = 1-aa[0]*aa[1]*aa[2]
    assert delta != 0
    F = [v[i]**2-aa[i]*v[(i+1) % 3]*v[(i+2) % 3] for i in range(3)]
    one, zero = R.one, R.zero

    def mono(m):
        out = one
        for variable, exponent in zip(v, m):
            out *= variable**exponent
        return out

    def jac(m):
        if max(m) <= 1:
            return mono(m), [zero, zero, zero]
        support = [i for i, k in enumerate(m) if k]
        i = next(i for i, k in enumerate(m) if k >= 2)
        if len(support) == 1:
            h = [zero, zero, zero]
            h[i] = v[i]**(m[i]-2)
            tail = v[i]**(m[i]-2)*v[(i+1) % 3]*v[(i+2) % 3]*aa[i]
            if m[i] <= 3:
                return tail, h
            mm = list(m)
            mm[i] -= 2
            mm[(i+1) % 3] += 1
            mm[(i+2) % 3] += 1
            normal, more = jac(tuple(mm))
            return normal*aa[i], [h[k]+more[k]*aa[i] for k in range(3)]
        j = next(j for j in support if j != i)
        k = 3-i-j
        mm = list(m)
        mm[i] -= 2
        mm[j] -= 1
        factor = mono(mm).quo_ground(delta)
        h = [zero, zero, zero]
        h[i] = factor*v[j]
        h[j] = factor*aa[i]*v[k]
        h[k] = factor*aa[i]*aa[j]*v[i]
        return zero, h

    L = R.from_expr(lower)
    gi = [L.diff(variable) for variable in v]

    @lru_cache(None)
    def nf(m):
        normal, h = jac(m)
        assert mono(m)-normal == sum((F[i]*h[i] for i in range(3)), zero)
        if max(m) <= 1:
            return normal, tuple(h)
        tail = sum(((u*(v[i]*h[i]-h[i].diff(v[i]))-gi[i]*h[i]).quo_ground(3*diagonal[i])
                    for i in range(3)), zero)
        out = normal
        witness = [(u*h[i]).quo_ground(3*diagonal[i]) for i in range(3)]
        for exponent, coefficient in tail.items():
            mm = exponent[:3]
            assert sum(mm) < sum(m)
            nn, ww = nf(mm)
            factor = u**exponent[3]*coefficient
            out += factor*nn
            for i in range(3):
                witness[i] += factor*ww[i]
        return out, tuple(witness)

    def reduce(poly):
        out = zero
        witness = [zero, zero, zero]
        for exponent, coefficient in poly.items():
            nn, ww = nf(exponent[:3])
            factor = u**exponent[3]*coefficient
            out += factor*nn
            for i in range(3):
                witness[i] += factor*ww[i]
        assert u*(poly-out) == sum((u*witness[i].diff(v[i])
                                   +(3*diagonal[i]*F[i]+gi[i]-u*v[i])*witness[i]
                                   for i in range(3)), zero)
        return out

    poly = sum((diagonal[i]*v[i]**3 for i in range(3)), zero)-3*coupling*x*y*z+L

    def vector(normal):
        rows = [[] for _ in range(8)]
        for exponent, coefficient in normal.items():
            row = BASIS.index(exponent[:3])
            while len(rows[row]) <= exponent[3]:
                rows[row].append(fmpq(0))
            rows[row][exponent[3]] += fmpq(int(coefficient.numerator), int(coefficient.denominator))
        return [fmpq_poly(row) for row in rows]

    connection_columns = [reduce(poly*mono(m)) for m in BASIS]
    C = [vector(col) for col in connection_columns]
    seed = [[fmpq_poly([int(i == 0)]) for i in range(8)]]
    zeta = fmpq_poly([0, 1])
    for j in range(8):
        current = seed[-1]
        seed.append([-zeta*zeta*current[i].derivative()
                     +sum((C[k][i]*current[k] for k in range(8)), fmpq_poly()) for i in range(8)])
    for j in range(1, 9):
        for i, coefficient in enumerate(seed[j]):
            if coefficient:
                assert coefficient.degree() <= 3*j-GRADES[i]
                if i == 0:
                    assert coefficient.degree() <= 3*j-2
                if L == 0:
                    for exponent, scalar in enumerate(coefficient):
                        if scalar:
                            difference = 3*j-GRADES[i]
                            assert exponent <= difference <= 3*exponent
                            assert (difference-exponent) % 2 == 0
                            assert exponent >= (j if i == 0 else j+2-GRADES[i])
    assert seed[1] == vector(reduce(poly))
    assert seed[2] == vector(reduce(poly**2))
    assert seed[3] == vector(reduce(poly**3))
    point_rank = fmpq_mat([[seed[j][i](1) for j in range(8)] for i in range(8)]).rank()
    serialized = [[str(col) for col in connection_columns], str(L), str(coupling)]
    return {
        'label': label, 'diagonal': diagonal, 'coupling': coupling, 'lower_part': str(L),
        'basis': [list(m) for m in BASIS], 'basis_grades': GRADES,
        'connection_columns': [[[str(v) for v in coefficient] for coefficient in column] for column in C],
        'seed_columns_0_to_8': [[[str(v) for v in coefficient] for coefficient in column] for column in seed],
        'seed_row_max_degrees': [[int(coefficient.degree()) if coefficient else -1 for coefficient in col] for col in seed],
        'seed_rank_at_inverse_t_1': point_rank,
        'full_seed_rank_proved': point_rank == 8,
        'connection_divergence_columns_checked': 8,
        'independent_direct_seed_reductions_checked': [1, 2, 3],
        'graded_seed_columns_checked': 8,
        'homogeneous_exponent_support_checked': L == 0,
        'prior_matrix_hash_format': hashlib.sha256(json.dumps(serialized).encode()).hexdigest(),
        'normal_forms_cached': nf.cache_info().currsize}


def main():
    start = time.monotonic()
    x, y, z = s.symbols('x y z')
    dense = 2*x*x+3*x*y-5*x*z+7*y*y+11*y*z-13*z*z+17*x-19*y+23*z+29
    models = [model('diagonal_dense', [1, 1, 1], 0, dense),
              model('hesse_dense', [1, 1, 1], 2, dense),
              model('hesse_homogeneous', [1, 1, 1], 2, s.Integer(0)),
              model('asymmetric_homogeneous', [1, 2, 3], 2, s.Integer(0))]
    prior = json.loads((P/'tame_gaussian_verification.json').read_text())
    for current, old in zip(models[:3], prior['connection_models']):
        assert current['prior_matrix_hash_format'] == old['matrix_sha256']
    assert models[1]['full_seed_rank_proved'] and models[3]['full_seed_rank_proved']
    rows = []
    for count in range(8):
        for selected in combinations(range(1, 8), count):
            seed_rank = 1+count
            grade_sum = sum(GRADES[i] for i in selected)
            triangular = seed_rank*(seed_rank+1)//2
            A = 3*triangular-grade_sum-2
            homogeneous_min = triangular+2*(seed_rank-1)-grade_sum
            assert A-homogeneous_min == seed_rank*(seed_rank-1)
            assert A <= 94
            for j in range(1, seed_rank+1):
                degree = A+2-3*j
                assert A-degree == 3*j-2
                assert -j+(3*j-2) >= 0
            rows.append({'seed_rank': seed_rank, 'selected_rows': [0]+list(selected),
                         'grade_sum': grade_sum, 'zeroth_cofactor_degree_upper': A,
                         'homogeneous_zeroth_cofactor_valuation_lower': homogeneous_min,
                         'homogeneous_zeroth_polynomial_degree_upper': A-homogeneous_min})
    table = [max(row['zeroth_cofactor_degree_upper'] for row in rows if row['seed_rank'] == rank) for rank in range(1, 9)]
    assert table == [1, 6, 14, 25, 38, 54, 73, 94]
    general = []
    for dimension in range(1, 13):
        mu = 2**dimension
        sorted_grades = [degree for degree in range(dimension+1) for _ in range(comb(dimension, degree))]
        assert sum(sorted_grades) == dimension*mu//2 and len(sorted_grades) == mu
        assert all(degree <= rank for rank, degree in enumerate(sorted_grades[1:], 1))
        bounds = [3*rank*(rank+1)//2-sum(sorted_grades[1:rank])-2 for rank in range(1, mu+1)]
        expected = mu*(3*mu+3-dimension)//2-2
        assert max(bounds) == expected
        general.append({'variables': dimension, 'rank_upper': mu,
                        'smooth_cubic_gate': expected, 'homogeneous_smooth_cubic_gate': mu*(mu-1)})
    assert general[2]['smooth_cubic_gate'] == 94 and general[2]['homogeneous_smooth_cubic_gate'] == 56
    certificate = {'models': models, 'cofactor_grade_cases': rows, 'general_dimension_regressions': general,
                   'retained_receipts': {'tame_gaussian_verification.json': hashlib.sha256((P/'tame_gaussian_verification.json').read_bytes()).hexdigest()},
                   'scope': 'Exact models and degree arithmetic only; universal filtered reduction and retained analytic no-constant-period proof remain written.'}
    certificate_path = P/'smooth_cubic_seed_bounds_certificate.json'
    certificate_path.write_text(json.dumps(certificate, indent=2)+'\n')
    report = {
        'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'certificate_sha256': hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
        'sympy_version': s.__version__, 'runtime_seconds': round(time.monotonic()-start, 3),
        'basis_grades': GRADES, 'model_count': len(models), 'exact_connection_columns': 32,
        'exact_seed_columns_through_8': 32, 'direct_seed_reductions': 12,
        'all_grade_subsets_checked': len(rows), 'rank_specific_cubic_gates': table,
        'smooth_ternary_cubic_gate': 94, 'homogeneous_smooth_ternary_cubic_gate': 56,
        'general_dimension_regressions_through': 12,
        'full_seed_rank_8_models': [m['label'] for m in models if m['full_seed_rank_proved']],
        'optimality_proved': False, 'global1855_gate_improved': False,
        'universal_written_proof_certified_by_checker': False,
        'boundary': 'Exact twisted-divergence model columns, seed recurrences and cofactor arithmetic support the written universal filtration proof. The no-constant-period input is retained analytic mathematics; independent replay, external review and formal verification are not claimed.'}
    (P/'smooth_cubic_seed_bounds_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

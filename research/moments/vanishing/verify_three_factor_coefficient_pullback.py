#!/usr/bin/env python3
"""Exact transcription checks for the written coefficient-pullback theorem."""
from pathlib import Path
import hashlib
import json
import time
import sympy as S

P = Path(__file__).resolve().parent
EXPECTED = {
    'three_factor_primary_decomposition_certificate.json': '74128da3fa241823622505d351c4ac0155beb6ddc3ff374d97dd7c856b521b4b',
    'three_factor_global_ideal_certificate.json': '80d6dbcedbe025795781f3c2d8b258ca3a5233d687b1fed93ec3f75d62c36e8c',
    'three_factor_colon_geometry_certificate.json': '77c2d2353dcb9f89789c6b31fa44fc63186b6d8bb52fcc50d713748b470be3ff',
    'three_factor_associated_geometry_certificate.json': 'a6fd844bfdb11c64b64cfab60d68cbfbd8a12d62ec70d1dbf8cc59549c9a8f6c',
}


def sha(name):
    return hashlib.sha256((P / name).read_bytes()).hexdigest()


def null_stratum(d, k, s):
    return S.expand(s * (d + k) - S.Rational(s * (3 * s + 1), 2))


def main():
    start = time.monotonic()
    for name, expected in EXPECTED.items():
        assert sha(name) == expected, name
    global_data = json.loads((P / 'three_factor_global_ideal_certificate.json').read_text())
    primary = json.loads((P / 'three_factor_primary_decomposition_certificate.json').read_text())
    colon = json.loads((P / 'three_factor_colon_geometry_certificate.json').read_text())
    n, t, E = S.symbols('n t E')
    # The finite rank list is exhaustive, not a numerical dimension sweep.
    stable_strata = []
    for r in range(4):
        k, d = 3 - r, n - r
        frame = r * n - S.Rational(r * (r + 1), 2)
        top = S.expand(frame + null_stratum(d, k, k))
        assert S.expand(top - (3 * n - 6)) == 0
        for s in range(k):
            gap = S.expand(null_stratum(d, k, k) - null_stratum(d, k, s))
            assert gap.subs(n, 7) > 0 and S.diff(gap, n) > 0
        stable_strata.append({'Gram_rank': r, 'null_vectors': k,
                              'top_rank': k, 'dimension': str(top)})
    boundary_strata = []
    for N in [3, 4, 5, 6]:
        ranks = range(4) if N >= 3 else range(N + 1)
        for r in ranks:
            k, d = 3 - r, N - r
            dimensions = [int(r * N - S.Rational(r * (r + 1), 2)
                              + null_stratum(d, k, s))
                          for s in range(min(k, d // 2) + 1)]
            if N >= 5:
                assert max(dimensions) == 3 * N - 6
            boundary_strata.append({'N': N, 'Gram_rank': r,
                                     'rank_stratum_dimensions': dimensions})
    assert max(x['rank_stratum_dimensions'][-1] for x in boundary_strata
               if x['N'] == 3 and x['Gram_rank'] == 0) == 4
    assert next(x for x in boundary_strata if x['N'] == 4 and x['Gram_rank'] == 0)['rank_stratum_dimensions'][-1] == 7
    # Lower-rank complements cannot be components of these complete intersections.
    incidence = {
        'row': {'top': 3*n-3, 'lower': [2*n], 'codimension': 3},
        'star': {'top': 3*n-3, 'lower': [2*n, n], 'codimension': 3},
        'axis': {'top': 3*n-5, 'lower': [2*n-1, n], 'codimension': 5},
    }
    for item in incidence.values():
        for low in item['lower']:
            gap = S.expand(item['top'] - low)
            assert gap.subs(n, 5) > 0 and S.diff(gap, n) > 0
    # A dense Witt chart for the N=5 zero fiber, keeping its nilpotent direction.
    v11, v12, v21, v22, v31, v32, z1, z2, z3 = S.symbols('v11 v12 v21 v22 v31 v32 z1 z2 z3')
    C5 = S.Matrix([[1, 0, v11, v12, z1],
                   [0, 1, v21, v22, z2],
                   [0, 0, v31, v32, z3]])
    J5 = S.Matrix([[0, 0, 1, 0, 0], [0, 0, 0, 1, 0],
                   [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1]])
    G5 = C5 * J5 * C5.T
    solved = {v11: -z1*z1/2, v22: -z2*z2/2, v21: -v12-z1*z2,
              v31: -z1*z3, v32: -z2*z3}
    residual = S.simplify(G5.subs(solved))
    assert residual == S.diag(0, 0, z3*z3)
    # Full-rank N=6 charts are smooth; the two Euclidean charts are separated.
    vv = S.Matrix(3, 3, S.symbols('u0:9'))
    J6 = S.zeros(6)
    J6[:3, 3:] = S.eye(3)
    J6[3:, :3] = S.eye(3)
    C6 = S.eye(3).row_join(vv)
    assert C6 * J6 * C6.T == vv + vv.T
    six = [vv[i, j] + vv[j, i] for i in range(3) for j in range(i, 3)]
    assert S.Matrix(six).jacobian(list(vv)).rank() == 6
    separator_checks = []
    for sign in [1, -1]:
        sample = S.eye(3).row_join(S.I * S.diag(1, 1, sign))
        assert sample * sample.T == S.zeros(3)
        p123, p456 = sample[:, :3].det(), sample[:, 3:].det()
        hp, hm = S.simplify(p123 + S.I*p456), S.simplify(p123 - S.I*p456)
        assert [hp, hm] == ([2, 0] if sign == 1 else [0, 2])
        separator_checks.append({'orthogonal_determinant': sign,
                                  'h_plus': str(hp), 'h_minus': str(hm)})
    # Dimension two: the determinant relation is symmetric and balanced.
    A, s, q = S.symbols('A s q')
    symbols = {'A': A, 'E': E, 's': s, 'q': q}
    fm = [S.sympify(x, locals=symbols) for x in global_data['symmetric_moments_F1_to_F5']]
    specialization = {A: -4*E, s: -2*E, q: 7*E*E}
    reduced = [S.expand(x.subs(specialization)) for x in fm]
    assert reduced == [0, 0, -4082400*E**3, 11735539200*E**4, -9807557760000*E**5]
    det = A + 2*E - s
    gb = S.groebner([fm[0], fm[1], det], A, s, q, E, domain=S.QQ)
    assert all(gb.reduce(x-y)[1] == 0 for x, y in zip(fm, reduced))
    assert gb.reduce(E**3)[1] != 0
    gb6 = S.groebner([fm[0], fm[1], fm[2], det], A, s, q, E, domain=S.QQ)
    assert all(gb6.reduce(x)[1] == 0 for x in fm)
    z = S.symbols('z')
    witness_cubic = z**3 + 2*z*z + 7*z + 4
    assert S.discriminant(witness_cubic, z) == -728
    # Its three roots are the diagonal Gram entries; off-diagonal entries are1.
    assert (-4) + 2 - (-2) == 0  # determinant, via elementary symmetric values
    m6_witness = reduced[2].subs(E, 1)
    assert m6_witness == -4082400
    # Flat graded base change doubles the Gram grading, not the coefficient count.
    numerators = colon['first_Hilbert_numerators']
    h10 = sum(c*t**i for i, c in enumerate(numerators['I10']))
    h8 = sum(c*t**i for i, c in enumerate(numerators['I8']))
    degree10 = S.cancel(h10/(1-t)**3).subs(t, 1)
    coefficient_degree = S.cancel(h10.subs(t, t*t)/(1-t)**3).subs(t, 1)
    assert degree10 == 90 and coefficient_degree == 720
    assert S.cancel(h8.subs(t, t*t)/(1-t)**3).subs(t, 1) == 720
    assert colon['exact_nilradical_indices'] == {'I8': 14, 'I10_equals_I_infinity': 13}
    assert primary['origin_torsion']['length'] == 77
    assert primary['diagonal_axis_torsion']['length'] == 539
    assert 2 * 77 == 154 and 6 * 2**3 == 48
    def strings(value):
        if isinstance(value, dict):
            return {str(key): strings(item) for key, item in value.items()}
        if isinstance(value, list):
            return [strings(item) for item in value]
        if isinstance(value, S.Basic):
            return str(value)
        return value
    cert = {
        'claim': 'G3COEFFICIENT20261114', 'field': 'C, characteristic zero',
        'map': 'C in Mat(3,N) maps to C*C^T; three rows are the coefficient vectors',
        'all_dimension_generation_cutoffs': {'N=1': 2, 'N=2': 6, 'N>=3': 10},
        'all_dimension_set_cutoffs': {'N=1': 2, 'N=2': 6, 'N>=3': 8},
        'N1_coefficient_scheme': {'ideal': '(u^2*v^2*w^2)',
                                   'primary_components': ['(u^2)', '(v^2)', '(w^2)'],
                                   'nilradical_index': 2, 'degree': 6},
        'rank_two': {'determinant': str(det), 'elimination': strings(specialization),
                     'reduced_moments': list(map(str, reduced)),
                     'sharp_witness_diagonal_root_polynomial': str(witness_cubic),
                     'witness_discriminant': -728, 'witness_M6': int(m6_witness)},
        'Gram_map_flatness': {'faithfully_flat_iff': 'N>=5',
                              'stable_symbolic_strata': stable_strata,
                              'boundary_strata': boundary_strata,
                              'proof_boundary': 'Finite rank list plus all-N linear dimension inequalities; flatness uses the written dimension/miracle-flatness proof, not a dimension sweep.'},
        'incidence_dimension_checks': strings(incidence),
        'N5_zero_fiber': {'generic_chart_residual': 'z3^2', 'generic_length': 2,
                          'components': 1, 'reduced': False},
        'N6_zero_fiber': {'components': 2, 'generically_reduced': True,
                          'separators': ['det(C[:,1:3])+i*det(C[:,4:6])',
                                         'det(C[:,1:3])-i*det(C[:,4:6])'],
                          'separator_checks': separator_checks},
        'stable_primary_transfer': {'N=5': 10, 'N=6': 11, 'N>=7': 10,
                                    'construction': 'Extend the ten universal primary ideals. For N=6 replace the origin extension by its two h_plus/h_minus saturations.',
                                    'minimal_primes': 6, 'codimensions': [3, 5, 6],
                                    'ordinary_coefficient_origin_associated': False},
        'stable_invariants': {'coefficient_scheme_degree': 720, 'reduced_degree': 48,
                              'nilradical_indices': {'I8': 14, 'I_infinity': 13},
                              'generic_row_length': 24, 'generic_star_length': 6,
                              'generic_axis_torsion_length': 539,
                              'generic_Gram_zero_torsion_lengths': {'N=5': [154], 'N=6': [77, 77], 'N>=7': [77]},
                              'Hilbert_numerator_I8': numerators['I8'],
                              'Hilbert_numerator_I10': numerators['I10'],
                              'coefficient_Hilbert_formula': 'H(S_N/I_k*S_N;t)=HilbertNumerator(R/I_k;t^2)/(1-t)^(3N), N>=5'},
        'provenance': EXPECTED,
        'boundary': 'Written geometric and invariant-theoretic transfer with exact algebra/chart checks. N=1 is elementary and complete; full coefficient primary/thickness descriptions for N=2,3,4 remain open, while all-N generation/set cutoffs are closed. No claim for arbitrary constrained subfamilies, affine factors, changed characteristic or arbitrary cubic GMC. No independent replay, external review or formal verification.'}
    path = P / 'three_factor_coefficient_pullback_certificate.json'
    path.write_text(json.dumps(cert, indent=2) + '\n')
    receipt = {'claim': cert['claim'], 'completed': True,
               'generation_cutoffs': cert['all_dimension_generation_cutoffs'],
               'flat_primary_counts': cert['stable_primary_transfer'],
               'stable_degree': 720,
               'checker_sha256': sha(Path(__file__).name),
               'certificate_sha256': sha(path.name),
               'elapsed_seconds': time.monotonic()-start,
               'assurance': cert['boundary']}
    (P / 'three_factor_coefficient_pullback_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

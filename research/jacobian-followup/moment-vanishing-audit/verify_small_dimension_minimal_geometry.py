#!/usr/bin/env python3
"""Exact checks for scoped N=3,4 geometry; geometric arguments remain written."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import re
import time
import sympy as S
from verify_nodal_product_reduction import moments

P = Path(__file__).resolve().parent
EXPECTED = {
    'verify_nodal_product_reduction.py': '7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6',
    'three_factor_primary_decomposition_certificate.json': '74128da3fa241823622505d351c4ac0155beb6ddc3ff374d97dd7c856b521b4b',
    'three_factor_associated_geometry_certificate.json': 'a6fd844bfdb11c64b64cfab60d68cbfbd8a12d62ec70d1dbf8cc59549c9a8f6c',
    'three_factor_colon_geometry_certificate.json': '77c2d2353dcb9f89789c6b31fa44fc63186b6d8bb52fcc50d713748b470be3ff',
}


def sha(name):
    return hashlib.sha256((P / name).read_bytes()).hexdigest()


def clean(name):
    return [x.strip() for x in (P / name).read_text().splitlines()
            if x.strip() and not x.lstrip().startswith('//')]


def run_check(stem, accepted=True):
    rec = json.loads((P / (stem + '_run.json')).read_text())
    assert rec['terminal'] and rec['completed'] == accepted
    assert rec['input_sha256'] == sha(stem + '.sing')
    assert rec['log_sha256'] == sha(stem + '_stdout.txt')
    if accepted:
        assert not rec['timed_out'] and rec['exit_code'] == 0
        assert '?' not in (P / (stem + '_stdout.txt')).read_text()
    return rec


def main():
    start = time.monotonic()
    for name, expected in EXPECTED.items():
        assert sha(name) == expected, name
    accepted = ['ternary_minimal_adic_geometry', 'origin_torsion_linear_map']
    failed = ['ternary_minimal_slices', 'ternary_quaternary_null_section',
              'rank_one_linear_colon', 'quaternary_hilbert_section',
              'unmixed_depth_audit']
    for stem in accepted:
        run_check(stem)
    for stem in failed:
        rec = run_check(stem, False)
        assert rec['timed_out'] or stem == 'ternary_minimal_slices'
    run_check('ternary_minimal_slices_v2')
    run_check('ternary_isotropic_pair_local')

    x, y, z, t = S.symbols('x y z t')
    J3 = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    charts = {
        'ZERO_FACTOR': (S.Matrix([[x, y, z], [1, 0, 0], [0, 1, 0]]), S.eye(3)),
        'ISOTROPIC_ROW': (S.Matrix([[1, 0, x], [0, 1, y], [1, 1, z]]), J3),
        'ISOTROPIC_PAIR': (S.Matrix([[S.Rational(1, 2), 0, 1], [1, 0, x], [1, z, y]]), J3),
    }
    source = (P / 'ternary_minimal_adic_geometry.sing').read_text()
    log = clean('ternary_minimal_adic_geometry_stdout.txt')
    expected_hf = {'ZERO_FACTOR': [1, 3, 5, 7, 8, 8, 7, 5, 3, 1],
                   'ISOTROPIC_ROW': [1, 2, 3, 4, 4, 4, 3, 2, 1],
                   'ISOTROPIC_PAIR': [1, 2, 2, 2, 2, 2, 1]}
    slices = {}
    for label, (C, J) in charts.items():
        gram = C * J * C.T
        block = source.split('ring r' + label + '=')[1]
        actual = block.split('ideal I=')[1].split(';')[0].split(',')
        expected = moments(gram.tolist(), 10)[1::2]
        assert len(actual) == len(expected) == 5
        for f, g in zip(actual, expected):
            f = S.sympify(f.replace('^', '**'), locals={'x': x, 'y': y, 'z': z})
            assert not f.has(S.Float)
            assert S.Poly(f, x, y, z).monic() == S.Poly(g, x, y, z).monic()
        length = int(log[log.index(label + '_LENGTH') + 1])
        steps = [[int(a) for a in log[i+1:i+4]] for i, marker in enumerate(log)
                 if marker == label + '_ADIC']
        hf = [steps[0][1]] + [b[1]-a[1] for a, b in zip(steps, steps[1:])]
        assert hf == expected_hf[label] and sum(hf) == length
        assert steps[-1] == [len(hf), length, 0]
        assert all(a[2] > 0 for a in steps[:-1])
        if label == 'ISOTROPIC_ROW':
            assert 'G=std(sat(G,ideal(2*z + 1)))' in block
        elif label == 'ISOTROPIC_PAIR':
            assert 'std(sat(I,ideal(y+2)))' in block
        else:
            # Finite homogeneous CI at this parameter proves the CI condition
            # on a nonempty open set, beyond the sampled all-five ideal.
            uu, vv, ww = S.symbols('uu vv ww')
            squared = []
            for f in expected[:3]:
                squared.append(sum(coef*uu**(mon[0]//2)*vv**(mon[1]//2)*ww**(mon[2]//2)
                                   for mon, coef in S.Poly(f, x, y, z).terms()))
                assert all(all(k % 2 == 0 for k in mon) for mon in S.Poly(f, x, y, z).monoms())
            assert S.groebner(squared, uu, vv, ww).is_zero_dimensional
        # A full-rank group-and-slice Jacobian makes the orbit argument explicit.
        orbit_rank = None
        if label != 'ZERO_FACTOR':
            C0 = C.subs({x: 0, y: 0, z: 0})
            vectors = []
            for i, j in [(0, 1), (0, 2), (1, 2)]:
                K = S.zeros(3); K[i, j] = 1; K[j, i] = -1
                omega = K * J.inv()
                assert omega*J + J*omega.T == S.zeros(3)
                vectors.append(S.Matrix(list(C0 * omega)))
            for i in range(3):
                scale = S.zeros(3); scale[i, i] = 1
                vectors.append(S.Matrix(list(scale * C0)))
            for v in [x, y, z]:
                vectors.append(S.Matrix(list(C.diff(v))))
            jac = S.Matrix.hstack(*vectors)
            orbit_rank = int(jac.rank())
            assert orbit_rank == 9 and S.Matrix.hstack(*vectors[:6]).rank() == 6
        slices[label] = {'Gram': [[str(q) for q in row] for row in gram.tolist()],
                         'length': length, 'maximal_ideal_Hilbert_function': hf,
                         'maximal_ideal_nilpotence_index': len(hf),
                         'group_and_slice_Jacobian_rank': orbit_rank}
    # Rank-stratum enumeration is exhaustive for these two fixed dimensions.
    flat = {}
    for N in [3, 4]:
        dims = {}
        for r in range(4):
            k, d = 3-r, N-r
            vals = [r*N-r*(r+1)//2 + s*(d+k)-s*(3*s+1)//2
                    for s in range(min(k, d//2)+1)]
            dims[r] = max(vals)
        assert dims == ({0: 4, 1: 4, 2: 3, 3: 3} if N == 3 else
                        {0: 7, 1: 6, 2: 6, 3: 6})
        flat[N] = dims
    # Dense N=4 Witt chart; six Gram entries have exactly this one relation.
    s, u = S.symbols('s u')
    v = S.Matrix(3, 2, S.symbols('v0:6'))
    U = S.Matrix([[1, 0], [0, 1], [s, u]])
    C4 = U.row_join(v)
    J4 = S.zeros(4); J4[:2, 2:] = S.eye(2); J4[2:, :2] = S.eye(2)
    g = C4*J4*C4.T
    relation = g[2, 2]+s*s*g[0, 0]+u*u*g[1, 1]+2*s*u*g[0, 1]-2*s*g[0, 2]-2*u*g[1, 2]
    assert S.expand(relation) == 0
    assert S.Matrix([g[0, 0], g[1, 1], g[0, 1], g[0, 2], g[1, 2]]).jacobian(list(v)).rank() == 5
    # A null pair in a nondegenerate three-space has generic fiber length2.
    q = S.symbols('q')
    pair = S.Matrix([[1, 0, x], [q, z, y]])
    gp = pair*J3*pair.T
    gbpair = S.groebner([gp[0, 0], gp[0, 1], gp[1, 1]], x, y, z, domain=S.QQ.frac_field(q))
    assert [f.as_expr() for f in gbpair.polys] == [x, y, z*z]
    # Check that the retained standard basis and the accepted new input agree.
    matrix_source = (P / 'origin_torsion_linear_map.sing').read_text()
    original_prefix = (P / 'three_factor_colon_geometry_v3.sing').read_text().split('ideal G8=std(I);')[0]
    assert matrix_source.startswith(original_prefix)
    sat_basis = (P / 'three_factor_origin_coordinate_stdout.txt').read_text().split('ORIGIN_SAT_BASIS\n')[1].split('FULL_HILBERT')[0].strip()
    assert matrix_source.split('ideal GS=')[1].split(';')[0].strip() == sat_basis
    lines = clean('origin_torsion_linear_map_stdout.txt')
    assert int(lines[lines.index('TORSION_BASIS_SIZE')+1]) == 77
    i = lines.index('TORSION_BASIS_DEGREES')+1
    degrees = list(map(int, lines[i:i+77])); M = S.zeros(77)
    for i, marker in enumerate(lines):
        if marker == 'MATRIX_ENTRY':
            row, col = [int(q)-1 for q in lines[i+1:i+3]]
            assert M[row, col] == 0 and degrees[row] == degrees[col]+1
            M[row, col] = S.Rational(lines[i+3])
    primary = json.loads((P / 'three_factor_primary_decomposition_certificate.json').read_text())
    ht = primary['origin_torsion']['Hilbert_function']
    assert Counter(degrees) == {i: a for i, a in enumerate(ht) if a}
    ranks, kernel = {}, {}
    for d in sorted(set(degrees)):
        cols = [i for i, a in enumerate(degrees) if a == d]
        rows = [i for i, a in enumerate(degrees) if a == d+1]
        rank = int(M.extract(rows, cols).rank()); ranks[d] = rank
        if len(cols) > rank:
            kernel[d] = len(cols)-rank
    assert kernel == {20: 10, 21: 5, 30: 3}
    Hk = sum(a*t**d for d, a in kernel.items())
    Ht = sum(a*t**d for d, a in enumerate(ht))
    assert S.expand((1-t)*Ht+t*Hk) == 9*t**18+6*t**19+3*t**20
    colon = json.loads((P / 'three_factor_colon_geometry_certificate.json').read_text())
    num = colon['first_Hilbert_numerators']['I10']
    Nh = S.Poly(S.expand(sum(a*t**d for d, a in enumerate(num))+t*Hk*(1-t)**5), t)
    section_num = [int(Nh.nth(i)) for i in range(Nh.degree()+1)]
    assert section_num == json.loads((P / 'quaternary_hilbert_section_run.json').read_text())['Hilbert_numerator']
    axis_length = primary['diagonal_axis_torsion']['length']
    lower_bound = 3*axis_length + sum(kernel.values())
    assert lower_bound == 1635
    assert 3*(48+24*6+12*4) == 720 and 3*(1+6+4) == 33
    assert 3*24*8+6*6*4 == 720 and 3*8+6*4 == 48
    provenance_names = list(EXPECTED) + ['three_factor_origin_coordinate_stdout.txt',
        'three_factor_colon_geometry_v3.sing', 'three_factor_minimal_primary_stdout.txt',
        'origin_linear_map_analysis.json']
    for stem in accepted + failed + ['ternary_minimal_slices_v2', 'ternary_isotropic_pair_local']:
        provenance_names += [stem+'.sing', stem+'_stdout.txt', stem+'_run.json', 'run_'+stem+'.py']
    cert = {
        'claim': 'G34MINIMAL20261114', 'field': 'C, characteristic zero',
        'N3': {'minimal_components': 9, 'components_each_type': 3,
               'types': ['zero factor', 'nonzero isotropic row closure', 'isotropic parallel pair'],
               'component_degrees': [1, 6, 4], 'generic_lengths': [48, 24, 12],
               'generic_maximal_ideal_indices': [10, 9, 7],
               'degree': 720, 'reduced_degree': 33, 'dimension': 6},
        'N4': {'minimal_components': 9, 'row_components': 3, 'chiral_star_components': 6,
               'component_degrees': [8, 4], 'generic_lengths': [24, 6],
               'degree': 720, 'reduced_degree': 48, 'dimension': 9,
               'associated_primes_on_Gram_nonzero_open': 12,
               'horizontal_axis_components': 3,
               'horizontal_axis_intrinsic_generic_torsion': 1078,
               'rank_two_Gram_zero_components_are_embedded_associated': 2,
               'full_associated_prime_count_lower_bound': 14},
        'fiber_dimensions_by_Gram_rank': flat,
        'flat_loci': {'N3': 'Gram rank >=2', 'N4': 'Gram nonzero'},
        'ternary_slices': slices,
        'quaternary_generic_null_section': {
            'ideal': 'I10 specialized by c=-a-b-2d-2e-2f',
            'associated_primes': 'six restrictions of the minimal Gram primes and the origin, exactly seven',
            'kernel_Hilbert_function': kernel, 'kernel_length': 18,
            'torsion_quotient_Hilbert': {18: 9, 19: 6, 20: 3},
            'origin_torsion_length_lower_bound': lower_bound,
            'first_Hilbert_numerator_over_5_variables': section_num,
            'matrix_degree_ranks': ranks,
            'full_origin_torsion_length': 'UNKNOWN;1635 is only a proved lower bound'},
        'provenance': {name: sha(name) for name in provenance_names},
        'boundary': 'Written geometric/orbit/flat-locus proofs and exact symbolic/CAS checks. Full N=3,4 associated-prime lists, primary decompositions and global nilradical indices remain open. The generic null-section has exactly seven associated primes but its full intrinsic origin torsion remains unknown. No result is inferred from a timed-out basis/colon run. No independent replay, external review or formal verification.'}
    path = P / 'small_dimension_minimal_geometry_certificate.json'
    path.write_text(json.dumps(cert, indent=2)+'\n')
    receipt = {'claim': cert['claim'], 'completed': True, 'checker_sha256': sha(Path(__file__).name),
               'certificate_sha256': sha(path.name), 'minimal_component_counts': [9, 9],
               'degrees': [720, 720], 'reduced_degrees': [33, 48],
               'section_kernel_length': 18, 'section_origin_torsion_lower_bound': lower_bound,
               'elapsed_seconds': time.monotonic()-start, 'assurance': cert['boundary']}
    (P / 'small_dimension_minimal_geometry_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()

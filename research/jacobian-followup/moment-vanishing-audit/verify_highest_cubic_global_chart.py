"""Replay complete highest-cubic coverage and a global finite-module bound.

The finite-field Groebner step is an explicit Singular dependency. A
degreewise integer-matrix rank argument, not an affine point lift, gives
the characteristic-zero upper bound. No purity classification or new
all-marker safety theorem is asserted here.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import subprocess
import time

import sympy as S

from conformal_euler_cubic import harmonic_basis, harmonic_parts, laplace, q, sphere_average, z, w, t
from highest_cubic_chart import (
    LW, LT, LZ, BM2, BM1, B0, B1, B2, CM3, CM2, CM1, C0, C1,
    VARIABLES13, WEIGHTS13, VARIABLES, WEIGHTS, BASE, DEPENDENT,
    DEPENDENT_WEIGHTS, phase13, eliminated_cminus3, phase, central_phase,
)
from highest_cubic_fibre import exact_moment

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def lower(f, parameter):
    return S.expand(f.subs({z: z-parameter*t-parameter**2*w/2,
                           t: t+parameter*w, w: w}, simultaneous=True))


def coeff(f, monomial):
    return S.Poly(f, z, w, t).coeff_monomial(monomial)


def source_polynomial(row, variables):
    output = []
    for powers, coefficient in row['terms']:
        monomial = '*'.join(v+(f'^{k}' if k != 1 else '')
                            for v, k in zip(variables, powers) if k)
        output.append(str(coefficient)+('*'+monomial if monomial else ''))
    return '+'.join(output).replace('+-', '-') or '0'


def singular_source(data, prime):
    return ('ring R='+str(prime)+',('+','.join(data['variables'])+'),wp('
            +','.join(map(str, data['weights']))+');\noption(redSB);\nideal I='
            +','.join(source_polynomial(row, data['variables']) for row in data['moments'])
            +';\nideal G=std(I);print("DIM");print(dim(G));print("VDIM");print(vdim(G));'
            +'print("LEADING_BEGIN");lead(G);print("LEADING_END");'
            +'if(dim(G)==0){ideal K=kbase(G);print("FOOTPRINT_BEGIN");K;'
            +'print("FOOTPRINT_END");};quit;\n')


def monomial_powers(expression, variables):
    result = [0]*len(variables)
    if expression == '1':
        return tuple(result)
    for factor in expression.split('*'):
        name, *power = factor.split('^')
        assert name in variables and len(power) <= 1
        result[variables.index(name)] += int(power[0]) if power else 1
    return tuple(result)


def inspect_log(text, data, expected_length, expected_top):
    assert re.search(r'(?m)^DIM\n0\nVDIM\n'+str(expected_length)+r'\n', text)
    assert '?' not in text and 'error' not in text.lower()
    footprint = []
    for line in text.split('FOOTPRINT_BEGIN\n')[1].split('FOOTPRINT_END')[0].splitlines():
        if line.startswith('K['):
            footprint.append(monomial_powers(line.split('=', 1)[1], data['variables']))
    assert len(footprint) == len(set(footprint)) == expected_length
    hilbert = Counter(sum(e*v for e, v in zip(pp, data['weights'])) for pp in footprint)
    assert max(hilbert) == expected_top
    leading = []
    for line in text.split('LEADING_BEGIN\n')[1].split('LEADING_END')[0].splitlines():
        if line.startswith('_['):
            leading.append(monomial_powers(line.split('=', 1)[1], data['variables']))
    # Pure powers in the initial ideal are a transparent check of dimension zero.
    pure_bounds = []
    for j in range(len(data['variables'])):
        pure_bounds.append(min(pp[j] for pp in leading
                               if pp[j] and sum(bool(e) for e in pp) == 1))
    # Downward closure catches a malformed/truncated footprint independently.
    points = set(footprint)
    for pp in footprint:
        for j, e in enumerate(pp):
            if e:
                down = list(pp)
                down[j] -= 1
                assert tuple(down) in points
    return {'length': expected_length, 'maximum_weight': expected_top,
            'leading_monomial_count': len(leading), 'pure_initial_power_bounds': pure_bounds,
            'weighted_hilbert_function': [[d, h] for d, h in sorted(hilbert.items())],
            'generating_monomials': [list(pp) for pp in footprint]}


def run(replay_cas=True):
    started = time.monotonic()
    s, u, rho = S.symbols('s u rho', nonzero=True)
    hb = harmonic_basis()
    assert all(laplace(f) == 0 for level in hb.values() for f in level.values())
    assert lower(q, s) == q
    jacobian = S.Matrix([lower(v, s) for v in (z, w, t)]).jacobian((z, w, t))
    assert jacobian.det() == 1
    swap = {z: w, w: z, t: -t}
    assert S.expand(q.subs(swap, simultaneous=True)) == q
    assert S.Matrix([w, z, -t]).jacobian((z, w, t)).det() == 1

    # Seven FIXED rational rotations cover H3 != 0. No root or support choice.
    low_polynomials = [coeff(lower(hb[3][j], s), w**3) for j in range(-3, 4)]
    assert all(f != 0 and S.Poly(f, s).degree() == j+3
               for j, f in zip(range(-3, 4), low_polynomials))
    evaluations = S.Matrix([[f.subs(s, k) for f in low_polynomials] for k in range(7)])
    cover_determinant = evaluations.det()
    assert cover_determinant != 0
    gamma = S.symbols('gamma_minus3 gamma_minus2 gamma_minus1 gamma0 gamma1 gamma2 gamma3')
    full_cubic = sum(a*hb[3][j] for a, j in zip(gamma, range(-3, 4)))
    moved = lower(full_cubic, u)
    assert coeff(moved, z**3) == gamma[6]
    assert coeff(moved, z*z*t) == gamma[5]-3*u*gamma[6]
    # In each gamma3 != 0 open, amplitude and this u are invertible adapters.

    F13 = phase13()
    mu2 = sphere_average(F13**2)
    assert S.expand(mu2-S.Rational(4, 35)*(CM3-eliminated_cminus3())) == 0
    assert sphere_average(phase()) == 0
    # Every term has coefficient weight + circular spatial charge exactly3.
    for pp, scalar in S.Poly(phase(), *VARIABLES, z, w, t).terms():
        assert sum(k*d for k, d in zip(pp[:len(VARIABLES)], WEIGHTS))+pp[-3]-pp[-2] == 3
    assert max(WEIGHTS) == 5 and min(WEIGHTS) == 1

    # Explicit finite nullcone incidence and its five retained base parameters.
    null_image = {
        LW: -LZ*s*s/2, LT: -LZ*s,
        BM2: -s**3*B1/2-3*s**4*B2/4,
        BM1: -3*s*s*B1/2-2*s**3*B2,
        B0: -s*B1-s*s*B2,
        CM2: -S.Rational(5, 2)*C1*s**3-S.Rational(9, 2)*s**5,
        CM1: -3*C1*s*s-S.Rational(9, 2)*s**4,
        C0: -s*C1-s**3,
    }
    positive_seed = (LZ*z+(B1+2*s*B2)*z*t+B2*z*z
                     +(C1+3*s*s)*hb[3][1]+3*s*z*z*t+z**3)
    assert S.expand(phase().subs(null_image, simultaneous=True)-lower(positive_seed, s)) == 0
    assert S.expand(eliminated_cminus3().subs(null_image, simultaneous=True)
                    +S.Rational(5, 8)*C1*s**4+S.Rational(5, 4)*s**6) == 0
    null_variables = (LZ, B1, B2, C1, s)
    null_coordinates = [S.expand(v.subs(null_image, simultaneous=True)) for v in VARIABLES]
    null_jacobian = S.Matrix(null_coordinates).jacobian(null_variables)
    assert null_jacobian.subs({LZ: 1, B1: 1, B2: 1, C1: 0, s: 0}).rank() == 5

    # The complete normalized Rodrigues graph is a second safe family.
    v = S.symbols('v')
    rod_image = {LW: -2*u*u, LT: 2*u*v, LZ: v*v,
                 BM2: 0, BM1: 0, B0: 0, B1: -2*u, B2: -2*v,
                 CM2: 0, CM1: 0, C0: 0, C1: 0}
    rod = S.expand(phase().subs(rod_image, simultaneous=True))
    assert rod == S.expand(z**3-2*u*z*t-2*v*z*z-2*u*u*w+2*u*v*t+v*v*z)
    h = (z-v)/u
    assert S.expand(rod+2*u*u*(w+h*t-z*h*h/2)) == 0

    # The imported odd Gaussian classification also gives exact odd-sphere
    # nullcone coverage, not merely a cutoff. Check its homogeneous adapter.
    aa, bb, cc, dd = S.symbols('odd_a odd_b odd_c odd_d')
    oddG = z*(aa*z*z+bb*z*t+cc*t*t+dd*z*w)
    odd_parts = harmonic_parts(oddG)
    assert S.expand(odd_parts[1]-(cc+2*dd)*z/5) == 0
    assert S.expand(odd_parts[3]-aa*z**3-bb*z*z*t
                    -(4*cc-2*dd)*hb[3][1]/5) == 0
    assert S.Matrix([[(S.diff(f, x)) for x in (cc, dd)]
                     for f in ((cc+2*dd)/5, (4*cc-2*dd)/5)]).det() == -S.Rational(2, 5)

    # Dense formula, low independent spatial checks, then all retained inputs.
    Fc = central_phase()
    assert sphere_average(Fc) == 0 and sphere_average(Fc**2) == 0
    test_point = dict(zip(DEPENDENT, range(1, 8)))
    Fnumeric = Fc.subs(test_point)
    for n in range(1, 7):
        computed = sum(S.Rational(c.numerator, c.denominator)
                       *math.prod(test_point[x]**e for x, e in zip(DEPENDENT, pp))
                       for pp, c in exact_moment(n).items())
        assert computed == sphere_average(Fnumeric**n)
    results = {}
    for prefix, expected_length, expected_top in [(14, 15638, 60), (18, 15528, 56)]:
        data = json.loads((HERE/f'highest_cubic_fibre_prefix{prefix}.json').read_text())
        assert data['variables'] == list(map(str, DEPENDENT))
        assert data['weights'] == list(DEPENDENT_WEIGHTS)
        assert data['base_variables'] == list(map(str, BASE))
        assert [row['n'] for row in data['moments']] == list(range(3, prefix+1))
        for row in data['moments']:
            scale = Fraction(row['raw_equals_primitive_times'])
            decoded = {tuple(pp): scale*c for pp, c in row['terms']}
            assert len(decoded) == len(row['terms']) and decoded == exact_moment(row['n'])
            assert all(sum(k*d for k, d in zip(pp, DEPENDENT_WEIGHTS)) == 3*row['n']
                       for pp in decoded)
        source_name = f'highest_cubic_fibre_prefix{prefix}_char65521.sing'
        log_name = f'highest_cubic_fibre_prefix{prefix}_char65521.log'
        assert (HERE/source_name).read_text() == singular_source(data, 65521)
        if replay_cas:
            result = subprocess.run(['Singular', '-q', str(HERE/source_name)],
                                    capture_output=True, text=True, timeout=45)
            assert result.returncode == 0 and not result.stderr
            text = result.stdout
            # The exact output, including the footprint, is deterministic here.
            assert text == (HERE/log_name).read_text()
        else:
            text = (HERE/log_name).read_text()
        record = inspect_log(text, data, expected_length, expected_top)
        record['input_sha256'] = digest(source_name)
        record['log_sha256'] = digest(log_name)
        record['source_moment_json_sha256'] = digest(f'highest_cubic_fibre_prefix{prefix}.json')
        results[str(prefix)] = record
    # No moments above18 can add a relation in the characteristic-zero fibre:
    # every homogeneous piece of weight >56 is already zero by the modular rank bound.
    assert 3*19 > results['18']['maximum_weight']

    receipt = {
        'claim': 'SIC2B33GLOBALGRADED20261117', 'status': 'PASS',
        'checker_sha256': digest(Path(__file__).name),
        'chart_model_sha256': digest('highest_cubic_chart.py'),
        'fibre_model_sha256': digest('highest_cubic_fibre.py'),
        'fixed_rational_rotation_parameters': list(range(7)),
        'rotation_cover_evaluation_determinant': str(cover_determinant),
        'normalized_chart_variables_before_mu2': list(map(str, VARIABLES13)),
        'normalized_chart_variables_after_mu2': list(map(str, VARIABLES)),
        'positive_weights': list(WEIGHTS), 'retained_base': list(map(str, BASE)),
        'base_finiteness': 'For N14 rank/generator upper bound15638, for N18 upper bound15528. No equality or freeness asserted.',
        'first14_global_zero_locus_dimension': 7,
        'all_order_global_pure_locus_dimension': 7,
        'odd_sphere_pure_classification': 'Exactly the one-sided nullcone; dimension5, using retained G3ODDPURE20261115 plus the checked homogeneous/polar adapter',
        'nilpotent_or_semisimple_quadratic_and_D_boundary': 'ALL INCLUDED',
        'central_fibre_complete_moment_ideal': 'I_infty=I18 over Q; finite support z^3; rational length at most15528',
        'central_fibre_ordinary_nilpotence_upper_bound': 19,
        'generic_component_safety': 'OPEN', 'full_Image33': 'OPEN',
        'Singular_dependency_for_finiteness': True, 'msolve_dependency': False,
        'CAS_replayed_this_run': replay_cas,
        'proof_boundary': 'A global dimension/finite-module theorem, not a global component decomposition, safe cover, exact characteristic-zero fibre length, local rigidity or all-marker safety. Whole chart remains positive dimensional.',
        'fibres': results, 'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'highest_cubic_global_chart_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'fibres'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--retained-cas-output', action='store_true',
                        help='check retained exact CAS evidence without rerunning Singular')
    args = parser.parse_args()
    run(not args.retained_cas_output)

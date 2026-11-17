"""Exact eigenvector-stratum exclusion under the cubic contraction condition.

This verifies the entire rank<=1 cyclic-vector stratum for semisimple
axial B. The mixed-eigenspace chart is deliberately NOT claimed closed.
All accepted exclusion witnesses are multiplied directly over QQ; no
standard-basis oracle is needed to check the proved result.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import time

import flint
import sympy as S

from conformal_euler_cubic import harmonic_basis, laplace, sphere_average, z, w, t
from axial_contraction_zero import (
    A, B, C, D, E, F, G, U, VARIABLES, exact_moment, raw_expression,
    amplitude_moment, linear_normalized_expression, central_linear_expression,
)
from nilpotent_kernel_unit_check import parse_polynomial

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    hb = harmonic_basis()
    gamma = S.symbols('gm3 gm2 gm1 g0 g1 g2 g3')
    full_cubic = sum(c*hb[3][j] for c, j in zip(gamma, range(-3, 4)))
    axial = t*t-z*w
    assert laplace(laplace(axial*full_cubic)) == (
        72*gamma[3]*t+24*gamma[4]*z+24*gamma[2]*w)
    qmatrix = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    operator = qmatrix*S.hessian(axial, (z, t, w))/2
    assert operator == S.diag(-S.Rational(1, 2), 1, -S.Rational(1, 2))
    vector = S.Matrix([C, B, A])
    wedge = vector.cross(operator*vector)
    assert wedge == S.Matrix([-3*A*B/2, 0, 3*B*C/2])
    # Over a field, wedge0 is exactly b0 OR(a=c0), including H1=0.

    data = json.loads((HERE/'axial_contraction_zero_prefix12.json').read_text())
    assert data['variables'] == list(map(str, VARIABLES))
    assert [r['n'] for r in data['moments']] == list(range(2, 13))
    for row in data['moments']:
        scale = Fraction(row['raw_equals_primitive_times'])
        decoded = {tuple(pp): scale*c for pp, c in row['terms']}
        assert len(decoded) == len(row['terms']) and decoded == exact_moment(row['n'])
    phase = A*z+B*t+C*w+axial+D*w**3+E*w*w*t+F*z*z*t+G*z**3
    values = dict(zip(VARIABLES, range(1, 8)))
    for n in range(1, 7):
        computed = sum(S.Rational(c.numerator, c.denominator)
                       *math.prod(values[v]**k for v, k in zip(VARIABLES, pp))
                       for pp, c in exact_moment(n).items())
        assert computed == sphere_average(phase.subs(values)**n)

    assert S.expand(raw_expression(2)-(2*A*C+B*B)/3-S.Rational(1, 5)
                    -4*E*F/105-4*D*G/35) == 0
    assert S.expand(raw_expression(3)-2*(B*B-A*C)/5-S.Rational(2, 35)
                    +4*D*G/35) == 0

    # Central linear eigenvector: products, then a direct univariate unit.
    central = json.loads((HERE/'axial_contraction_zero_central_unit.json').read_text())
    central_polynomials = [central_linear_expression(n) for n in (4, 5)]
    assert central_polynomials == [S.sympify(f, locals={'b': B}) for f in central['moments']]
    assert S.expand(sum(S.sympify(c, locals={'b': B})*f
                        for c, f in zip(central['multipliers'], central_polynomials))) == 1
    central2 = S.expand(raw_expression(2).subs({A: 0, C: 0}))
    central3 = S.expand(raw_expression(3).subs({A: 0, C: 0}))
    assert S.expand(central2.subs({D: ((1+7*B*B)/2)/G,
                    E: (-(27+77*B*B)/4)/F}, simultaneous=True)) == 0
    assert S.expand(central3.subs(D, ((1+7*B*B)/2)/G)) == 0
    # The displayed divisions are algebraic checks of polynomial product
    # substitutions; the proof uses products directly and does NOT assume
    # any individual d,e,f,g nonzero on this central branch.

    # Verify restoring the quadratic amplitude and dividing odd moments.
    lam = S.Symbol('lam', nonzero=True)
    test = {A: 2, D: 3, E: 4, F: 5, G: 6, U: 9}
    for central_value in (0, 1):
        current = 2*z+central_value*t+w+3*axial+3*w**3+4*w*w*t+5*z*z*t+6*z**3
        for n in range(1, 6):
            assert amplitude_moment(n, central_value).subs(test)*3**(n % 2) == sphere_average(current**n)
    assert S.expand(amplitude_moment(2, 1)-(2*A+1)/3-U/5-4*D*G/35-4*E*F/105) == 0
    assert S.expand(amplitude_moment(3, 1)-2*(1-A)/5-2*U/35+4*D*G/35) == 0

    normalized = json.loads((HERE/'axial_contraction_zero_linear_normalized_charts.json').read_text())
    for name, chart in normalized['charts'].items():
        for row in chart['moments']:
            if row['n'] <= 8:
                assert linear_normalized_expression(name, row['n']) == S.sympify(
                    row['primitive'], locals={str(v): v for v in VARIABLES})

    # Nonisotropic transverse eigenvector. No cubic coefficient is
    # required nonzero; a=c1 and the nonzero quadratic amplitude remain.
    relation = 217+54*D*G+4*E*F
    assert S.expand(amplitude_moment(3, 0).subs({A: 1, U: 7+2*D*G})) == 0
    assert S.expand(amplitude_moment(2, 0).subs({A: 1, U: 7+2*D*G})-relation/105) == 0
    chart = normalized['charts']['transverse_nonisotropic']
    context = flint.fmpq_mpoly_ctx.get(tuple(chart['variables']))
    source_strings = chart['relations']+[row['primitive'] for row in chart['moments'][:5]]
    equations = [parse_polynomial(text, chart['variables'], context) for text in source_strings]
    log_name = 'axial_contraction_zero_transverse_nonisotropic_unit_witness.log'
    lines = (HERE/log_name).read_text().split('WITNESS_BEGIN\n')[1].split('WITNESS_END')[0].strip().splitlines()
    multipliers = [parse_polynomial(text, chart['variables'], context) for text in lines]
    assert len(multipliers) == len(equations) == 6
    assert sum((m*f for m, f in zip(multipliers, equations)), context.constant(0)) == context.constant(1)
    counts = [len(m) for m in multipliers]

    # Isotropic transverse eigenvector: dg1/2 and ef-27/4 force d,e
    # to be units. Exact polynomial rows give e^2=d^2=0, a contradiction.
    assert S.expand(raw_expression(2).subs({A: 0, B: 0, C: 1,
                    F: -S.Rational(27, 4)/E, G: S.Rational(1, 2)/D}, simultaneous=True)) == 0
    assert S.expand(raw_expression(3).subs({A: 0, B: 0, C: 1,
                    F: -S.Rational(27, 4)/E, G: S.Rational(1, 2)/D}, simultaneous=True)) == 0
    chart = normalized['charts']['transverse_isotropic']
    context = flint.fmpq_mpoly_ctx.get(tuple(chart['variables']))
    equations = [parse_polynomial(row['primitive'], chart['variables'], context)
                 for row in chart['moments'][:3]]
    totals = [context.constant(0), context.constant(0)]
    for line in (HERE/'axial_contraction_zero_transverse_isotropic_witness.log').read_text().split(
            'WITNESS_BEGIN\n')[1].split('WITNESS_END')[0].strip().splitlines():
        i, j, value = line.split('|')
        totals[int(j)-1] += parse_polynomial(value, chart['variables'], context)*equations[int(i)-1]
    assert totals == [context.gen(1)**2, context.gen(0)**2]

    # Export the closed safe cone by importing the already certified
    # seven equations of the harmonic-square quadratic cone, dropping
    # the old H1=0 restrictions, and adding only the declared conditions.
    cone = json.loads((HERE/'boundary_axial_safe_container.json').read_text())
    assert cone['base_generators'][:3] == ['h1_z', 'h1_t', 'h1_w']
    quadratic_names = S.symbols('h2_w2 h2_wt h2_t2minuszw h2_zt h2_z2')
    linear_names = S.symbols('h1_z h1_t h1_w')
    quadratic = sum(v*hb[2][j] for v, j in zip(quadratic_names, range(-2, 3)))
    matrix = qmatrix*S.hessian(quadratic, (z, t, w))/2
    vec = S.Matrix([linear_names[2], linear_names[1], linear_names[0]])
    wedge_generators = [str(S.expand(f)) for f in vec.cross(matrix*vec)]
    contracted = laplace(laplace(quadratic*full_cubic))
    contracted_poly = S.Poly(contracted, z, t, w)
    contraction_generators = [str(contracted_poly.coeff_monomial(v)) for v in (z, t, w)]
    container = {
        'claim': 'SIC2B33AXIALCOVEIG20261117',
        'status': 'CLOSED_PURE_SAFE_CONTAINER, NOT_ENTIRELY_PURE',
        'variables': list(map(str, (*linear_names, *quadratic_names, *gamma))),
        'constant_harmonic': '0 by the first moment',
        'quadratic_cone_source': 'boundary_axial_safe_container.json',
        'quadratic_cone_source_sha256': digest('boundary_axial_safe_container.json'),
        'quadratic_cone_generators': cone['base_generators'][3:],
        'linear_eigenvector_generators': wedge_generators,
        'cubic_contraction_generators': contraction_generators,
        'linear_harmonic': 'arbitrary subject to v wedge Bv=0',
        'cubic_harmonic': 'arbitrary subject to Delta^2(H2*H3)=0',
        'pure_field_classification': 'Every pure point is one-sided; nonzero semisimple axial part has no first8 point; nilpotent rank1 and B0 import retained theorems.',
        'full_Image_cutoff': 'm>coordinate_degree Q',
        'entire_cone_pure': False, 'unrestricted_mixed_eigenspace_axial_case': 'OPEN',
        'unrestricted_Image33': 'OPEN',
    }
    assert len(container['quadratic_cone_generators']) == 7
    (HERE/'axial_contraction_zero_eigenvector_safe_container.json').write_text(json.dumps(container, indent=2)+'\n')
    receipt = {
        'claim': container['claim'], 'status': 'PASS',
        'checker_sha256': digest(Path(__file__).name),
        'model_sha256': digest('axial_contraction_zero.py'),
        'raw_input_sha256': digest('axial_contraction_zero_prefix12.json'),
        'normalized_input_sha256': digest('axial_contraction_zero_linear_normalized_charts.json'),
        'raw_moments_regenerated': list(range(2, 13)),
        'independent_spatial_checks': list(range(1, 7)),
        'quadratic_amplitude_checks': list(range(1, 6)),
        'normalized_polynomials_regenerated': 'all three charts through8, including the unproved mixed chart as INPUT ONLY',
        'central_eigenvector_prefix': 5,
        'transverse_isotropic_eigenvector_prefix': 6,
        'transverse_nonisotropic_eigenvector_prefix': 8,
        'nonisotropic_exact_unit_multiplier_term_counts': counts,
        'nonisotropic_unit_sha256': digest(log_name),
        'isotropic_square_memberships_sha256': digest('axial_contraction_zero_transverse_isotropic_witness.log'),
        'closed_safe_container_sha256': digest('axial_contraction_zero_eigenvector_safe_container.json'),
        'Singular_dependency_for_new_exclusions': False,
        'msolve_dependency': False,
        'Singular_role': 'Produced candidate coefficients; accepted proof uses independent exact QQ polynomial multiplication.',
        'imported_boundary_scope': 'Harmonic-square cone description, rank1 contraction-zero safety, and all-order odd safety retain their source assurance. No inherited largest axial or rank2 QQ unit is used for the new exclusions.',
        'mixed_eigenspace_chart': 'OPEN; modular unit is not a QQ exclusion, and incomplete homogeneous bases are not GB certificates.',
        'global_three_memberships': 'ALL OPEN', 'full_Image33': 'OPEN',
        'boundary': 'Complete cyclic-vector rank<=1 exclusion on the semisimple axial contraction-zero stratum, plus a closed conditional safe cone. No unrestricted axial purity exclusion, full Image, coefficient-ring safety or independent specialist/external/formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'axial_contraction_zero_eigenvector_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()

#!/usr/bin/env python3
"""Exact exhaustive H1=0/rank-one-H2 certificate, without a CAS unit oracle.

All exclusion steps are exported polynomial identities after invertible
normalizations and monic eliminations. Does not solve the full (3,3) problem.
"""
import hashlib
import json
from pathlib import Path

from boundary_rank1_quadratic import (
    S, LOWEST, COEFFICIENTS, exact_moment, raw_expression, next_coefficient_chart,
)
from conformal_euler_cubic import z, t, w, q, sphere_average


def zero(expression):
    assert S.expand(expression) == 0, expression


def main():
    root = Path(__file__).resolve().parent
    cm3 = LOWEST
    cm2, cm1, y, x, c2, c3 = COEFFICIENTS
    all_variables = (cm3, *COEFFICIENTS)
    basis = (w**3, w*w*t, w*(t*t-z*w/2), t**3-3*z*t*w,
             z*(t*t-z*w/2), z*z*t, z**3)
    phase = w*w+sum(coefficient*marker for coefficient, marker in zip(all_variables, basis))
    for marker in basis:
        zero(2*S.diff(marker, z, w)+S.diff(marker, t, 2))
    control = dict(zip(all_variables, [2, 4, 6, 8, 10, 12, 14]))
    for n in range(1, 7):
        assert raw_expression(n).subs(control) == sphere_average(phase.subs(control)**n)
    s, rho = S.symbols('s rho', nonzero=True)
    rotation = {z: z-s*t-s*s*w/2, t: t+s*w, w: w}
    zero(q.subs(rotation, simultaneous=True)-q)
    transformed = S.Poly(S.expand(phase.subs(rotation, simultaneous=True)), z, t, w)
    assert transformed.coeff_monomial(z**3) == c3
    assert transformed.coeff_monomial(z*z*t) == c2-3*s*c3
    # When c3=0 the next charged harmonic has coefficient x-2s*c2.
    assert S.expand(transformed.coeff_monomial(z*t*t).subs(c3, 0)) == x-2*s*c2
    for charge, marker in zip(range(-3, 4), basis):
        zero(rho**2*marker.subs({z: rho*z, w: w/rho}, simultaneous=True)
             -rho**(charge+2)*marker)
    zero(rho**2*(w*w).subs(w, w/rho)-w*w)

    # No positive charge: the middle coefficient is forced zero by mu2.
    assert raw_expression(2).subs({x: 0, c2: 0, c3: 0}) == y*y/7
    # Highest positive charge1: mu3 is a nonzero scalar times its square.
    assert raw_expression(3).subs({c2: 0, c3: 0}) == 2*x*x/105

    # Highest charge2, normalized c2=1, and x=0 by a unipotent rotation.
    case2 = {n: next_coefficient_chart(n, 2).as_expr() for n in range(2, 6)}
    assert case2[2] == 15*y*y+4*cm2 and case2[3] == -y
    f4 = S.expand(case2[4].subs({y: 0, cm2: 0}))
    f5 = S.expand(case2[5].subs({y: 0, cm2: 0}))
    assert f4 == 16*cm1*cm3+208 and f5 == cm1*cm1
    h4, h5 = (13-cm1*cm3)/S.Integer(2704), cm3**2/S.Integer(169)
    zero(h4*f4+h5*f5-1)

    # Highest charge3, normalized c3=1 and c2=0. No denominator cases remain:
    # the first two equations solve cm1 and cm3 with constant pivots2 and12.
    moments = {n: next_coefficient_chart(n, 3).as_expr() for n in range(2, 11)}
    assert moments[2] == 15*y*y+5*x*cm1+12*cm3
    assert moments[3] == x*x-2*cm1
    substitutions = {cm1: x*x/2, cm3: -5*y*y/4-5*x**3/24}
    reduced, reduction_scales = {}, {}
    for n in range(4, 11):
        raw = S.Poly(moments[n].subs(substitutions, simultaneous=True), cm2, y, x)
        primitive = raw.clear_denoms()[1].primitive()[1]
        factor = raw.LC()/primitive.LC()
        assert factor != 0 and raw == factor*primitive
        reduced[n], reduction_scales[n] = primitive.as_expr(), factor
    u = S.Symbol('u')
    f = {n: S.expand(value.subs(cm2, u-y*x/2)) for n, value in reduced.items()}
    v = 27*y*y+4*x**3
    A = S.Matrix([u*u+5, v*v-2112*x, 287*x*v-2358*y*u,
                  x*x*v*v-36*y*x*u*v-1728*y*y-16*x**3])
    Q6 = 3888*y*y*x-1728*y*u+576*x**4
    transform = S.Matrix([[0, S.Rational(1, 4), 0, 0],
                          [1, -12*x, 0, 0],
                          [-v/12, (12*x*v-Q6/4)/12, S.Rational(1, 12), 0],
                          [0, -81*y*y, 0, 4]])
    assert transform.det() != 0 and not transform.det().free_symbols
    assert (transform*S.Matrix([f[n] for n in [4, 5, 6, 7]])-A).applyfunc(S.expand) == S.zeros(4, 1)

    # Explicit ideal memberships, rather than trusting the final Groebner basis.
    e0, e1, e2, e3 = [S.eye(4)[:, i] for i in range(4)]
    B_rep = -(2358*e3+7974*x*x*e1-36*x*v*e2)/288
    C_rep = (287*x*v+2358*y*u)*e2-287**2*x*x*e1+2358**2*y*y*e0
    B = 14148*y*y+58607*x**3
    C = 5*2358**2*y*y+287**2*2112*x**3
    zero((B_rep.T*A)[0]-B)
    zero((C_rep.T*A)[0]-C)
    linear = S.Matrix([[14148, 58607], [5*2358**2, 287**2*2112]])
    assert linear.det() == 831910506804
    inverse = linear.inv()
    y2_rep = inverse[0, 0]*B_rep+inverse[0, 1]*C_rep
    x3_rep = inverse[1, 0]*B_rep+inverse[1, 1]*C_rep
    v_rep = 27*y2_rep+4*x3_rep
    x_rep = (v*v_rep-e1)/2112
    yu_rep = (287*v*x_rep-e2)/2358
    y_rep = (y*e0-u*yu_rep)/5
    zero((x_rep.T*A)[0]-x)
    zero((y_rep.T*A)[0]-y)
    for generator in A:
        assert S.reduced(generator, [u*u+5, y, x], u, y, x)[1] == 0
    for n in [8, 9]:
        assert S.reduced(f[n], [u*u+5, y, x], u, y, x)[1] == 0
    quotients, remainder = S.reduced(f[10], [u*u+5, y, x], u, y, x)
    assert remainder == -771876864
    unit_rep = -(quotients[0]*e0+quotients[1]*y_rep+quotients[2]*x_rep)/remainder
    zero((unit_rep.T*A)[0]+f[10]/remainder-1)

    # The surviving first9 pair is an exact exclusion control, not pure.
    alpha = S.Symbol('alpha')
    candidate = w*w+z**3+alpha*w*w*t
    candidate_moments = []
    for n in range(1, 11):
        value = sphere_average(candidate**n)
        remainder_n = S.rem(S.Poly(value, alpha), S.Poly(alpha*alpha+5, alpha)).as_expr()
        candidate_moments.append(remainder_n)
    assert candidate_moments[:9] == [0]*9 and candidate_moments[9] != 0
    sharp_phase = w*w+w*(t*t-z*w/2)
    assert sphere_average(sharp_phase**2*z*z) == S.Rational(2, 315)

    certificate = {
        'field': 'QQ', 'case2_moments': {str(n): str(value) for n, value in case2.items()},
        'case2_reduced_sources': [str(f4), str(f5)],
        'case2_unit_multipliers': [str(h4), str(h5)],
        'case3_primitive_moments': {str(n): str(value) for n, value in moments.items()},
        'case3_monic_eliminations': {str(key): str(value) for key, value in substitutions.items()},
        'case3_reduction_scales': {str(n): str(value) for n, value in reduction_scales.items()},
        'case3_coordinate_change': 'u=cminus2+c0*c1/2',
        'case3_reduced_sources': {str(n): str(value) for n, value in f.items()},
        'case3_equivalent_generators': list(map(str, A)),
        'case3_generator_transform': [[str(value) for value in row] for row in transform.tolist()],
        'case3_x_membership': list(map(str, x_rep.applyfunc(S.expand))),
        'case3_y_membership': list(map(str, y_rep.applyfunc(S.expand))),
        'case3_unit_multipliers_on_A': list(map(str, unit_rep.applyfunc(S.expand))),
        'case3_unit_multiplier_on_moment10': str(1/remainder),
        'first9_candidate_condition': 'alpha^2+5=0',
        'first9_candidate_mu10': str(candidate_moments[9]),
    }
    output = root/'boundary_rank1_quadratic_certificate.json'
    output.write_text(json.dumps(certificate, indent=2)+'\n')
    hz, ht, hw = S.symbols('h1_z h1_t h1_w')
    b0, b1, b2, b3, b4 = S.symbols('h2_w2 h2_wt h2_t2minuszw h2_zt h2_z2')
    operator = S.Matrix([[-b2/2, b1/2, b0], [b3/2, b2, b1/2], [b4, b3/2, -b2/2]])
    minors = [S.expand(operator.extract(rows, columns).det())
              for rows in [(0, 1), (0, 2), (1, 2)]
              for columns in [(0, 1), (0, 2), (1, 2)]]
    container = {'coefficient_field': 'QQ', 'linear_harmonic': 'h1_z*z+h1_t*t+h1_w*w',
                 'quadratic_harmonic': 'h2_w2*w^2+h2_wt*w*t+h2_t2minuszw*(t^2-z*w)'
                                      '+h2_zt*z*t+h2_z2*z^2',
                 'cubic_harmonic': 'unrestricted',
                 'defining_equations': list(map(str, [hz, ht, hw, *minors])),
                 'scope': 'Every all-order pure FIELD point is safe. Rank1 nonzero '
                          'quadratic case proved by this certificate; zero quadratic '
                          'case uses the retained odd-harmonic theorem. This container '
                          'is not asserted to consist entirely of pure points and '
                          'does not cover the full pure locus.',
                 'full_Image_cutoff': 'm>coordinate_degree Q'}
    (root/'boundary_rank1_safe_container.json').write_text(json.dumps(container, indent=2)+'\n')
    report = {
        'claim': 'SIC2B33RANK1QUAD20261117', 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_sha256': hashlib.sha256((root/'boundary_rank1_quadratic.py').read_bytes()).hexdigest(),
        'certificate_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'complete_stratum': 'H1=0 and nonzero trace-free quadratic operator of rank1',
        'normalized_phase': 'w^2+arbitrary harmonic ternary cubic',
        'pure_classification': 'Exactly the cubics of torus charges-3,-2,-1 in an '
            'oriented linear frame with H2=w^2. No semistable pure point in this stratum.',
        'set_theoretic_purity_prefix': 10, 'least_consecutive_set_theoretic_prefix': 10,
        'case2_exclusion_prefix': 5, 'case3_exclusion_prefix': 10,
        'case3_prefix7_to9_reduced_algebra': 'QQ[u]/(u^2+5), x=y=0',
        'case3_moment10_unit_remainder': str(remainder),
        'first9_false_candidate_mu10': str(candidate_moments[9]),
        'full_Image_cutoff': 'm>coordinate_degree Q',
        'strict_cutoff_boundary_control': 'h((w^2+w*(t^2-zw/2))^2*z^2)=2/315',
        'msolve_dependency': False, 'global_unrestricted_Image33': 'OPEN',
        'boundary': 'Exact polynomial identities cover every highest-charge chart '
            'after a complete invertible normalization. Classification of this entire '
            'rank1 quadratic stratum is complete, not the whole D=0 boundary or '
            'full(3,3) Image theorem. No ordinary ideal stabilization is asserted. The '
            'zero-quadratic case is separately covered by the retained odd-slice '
            'theorem. No independent specialist or formal verification.'}
    (root/'boundary_rank1_quadratic_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

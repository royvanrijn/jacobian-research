"""Algebraic controls for the written top-two-boundary-jet proof.

The proper-graph/holomorphic-period argument is an internal written proof;
this script verifies its local formulas and the exact axial pole units.
It does not machine-certify the analytic or global geometric argument.
"""
from pathlib import Path
import hashlib
import json
import time

import sympy as S

from conformal_euler_cubic import harmonic_basis, laplace, z, t, w
from axial_contraction_zero import (
    A, B, C, D, E, F, G, U, raw_expression, amplitude_moment,
    linear_normalized_expression,
)

HERE = Path(__file__).resolve().parent
CLAIM = 'SIC2B33TOPTWOGATE20261117'


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    s, u, tau = S.symbols('s u tau')
    Z = (s*s-u*u)/u
    T = s/u
    W = -1/(2*u)
    assert S.cancel(2*Z*W+T*T) == 1
    jacobian = S.Matrix([[S.diff(T, v) for v in (u, s)],
                         [S.diff(W, v) for v in (u, s)]]).det()
    assert S.cancel(jacobian/(2*W)) == 1/(2*u*u)

    hb = harmonic_basis()
    cubic_coefficients = S.symbols('gm3 gm2 gm1 g0 g1 g2 g3')
    quadratic_coefficients = S.symbols('bm2 bm1 b0 b1 b2')
    cubic = sum(c*hb[3][j] for c, j in zip(cubic_coefficients, range(-3, 4)))
    quadratic = sum(c*hb[2][j] for c, j in zip(quadratic_coefficients, range(-2, 3)))
    linear = A*z+B*t+C*w
    substitutions = {z: s*s-u*u, t: s, w: -S.Rational(1, 2)}
    homogenized = S.expand((cubic+u*quadratic+u*u*linear).subs(substitutions))
    boundary_cubic = S.expand(cubic.subs(substitutions).subs(u, 0))
    boundary_quadratic = S.expand(quadratic.subs(substitutions).subs(u, 0))
    assert S.diff(homogenized, u).subs(u, 0) == boundary_quadratic
    assert S.expand(S.diff(homogenized, s).subs(u, 0)-S.diff(boundary_cubic, s)) == 0
    # K_tau=-u^3 gives the factor u^(3-2) in each of the two
    # implicit-function volume formulas. Unit Jacobians, smoothness,
    # properness and the global period theorem are proved in the note.
    K = homogenized-tau*u**3
    assert S.diff(K, tau) == -u**3

    axial = t*t-z*w
    contraction = laplace(laplace(axial*cubic))
    assert contraction == (24*cubic_coefficients[2]*w
                           +72*cubic_coefficients[3]*t
                           +24*cubic_coefficients[4]*z)
    restricted_cubic = D*w**3+E*w*w*t+F*z*z*t+G*z**3
    binary = S.expand(8*restricted_cubic.subs({z:s*s, t:s, w:-S.Rational(1, 2)}))
    assert binary == -D+2*E*s+8*F*s**5+8*G*s**6
    assert S.expand(axial.subs({z:s*s, t:s, w:-S.Rational(1, 2)})) == 3*s*s/2
    assert binary.subs(s, 0) == -D
    assert S.diff(binary, s).subs(s, 0) == 2*E
    reversed_binary = S.Poly(binary, s).homogenize(u).as_expr().subs({s:1, u:s}, simultaneous=True)
    assert reversed_binary.subs(s, 0) == 8*G
    assert S.diff(reversed_binary, s).subs(s, 0) == 8*F
    # The two projective poles therefore give d=e0 OR f=g0.

    mu2 = raw_expression(2).subs({D:0, E:0})
    mu3 = raw_expression(3).subs({D:0, E:0})
    ac, bsq = S.symbols('ac bsq')
    equations = [S.expand(v).subs({A*C:ac, B**2:bsq}) for v in (mu2, mu3)]
    solution = S.solve(equations, (ac, bsq))
    assert solution == {ac:-S.Rational(16,105), bsq:-S.Rational(31,105)}
    assert solution[ac]/solution[bsq] == S.Rational(16,31)
    assert 1/solution[bsq] == -S.Rational(105,31)
    for zero in ({D:0, E:0}, {F:0, G:0}):
        for n in (2, 3):
            assert amplitude_moment(n, 1).subs(
                {**zero, A:S.Rational(16,31), U:-S.Rational(105,31)}) == 0

    data = json.loads((HERE/'axial_contraction_zero_boundary_pole_units.json').read_text())
    counts = {}
    for name, zero, variables in [('minus_pole', {D:0,E:0}, (F,G)),
                                   ('plus_pole', {F:0,G:0}, (D,E))]:
        row = data['charts'][name]
        assert row['variables'] == list(map(str, variables))
        equations = [S.expand(linear_normalized_expression('mixed_linear', n).subs(zero))
                     for n in (4, 5, 6)]
        decoded = [S.sympify(f, locals={str(v):v for v in variables}) for f in row['equations']]
        assert decoded == equations
        root = dict(zip(variables, map(S.Rational, row['linear_root'])))
        assert [f.subs(root) for f in equations[:2]] == [0, 0]
        constant = equations[2].subs(root)
        assert constant == S.Rational(198023471767658880,121)
        assert constant == S.Rational(row['sixth_moment_at_root'])
        multipliers = [S.sympify(f, locals={str(v):v for v in variables}) for f in row['multipliers']]
        assert S.expand(sum(f*a for f, a in zip(equations, multipliers))) == 1
        counts[name] = [len(S.Poly(a,*variables).terms()) for a in multipliers]
        assert counts[name] == [2, 2, 1]

    old = json.loads((HERE/'axial_contraction_zero_eigenvector_safe_container.json').read_text())
    exported = {k:v for k,v in old.items() if k not in
                ('claim', 'linear_eigenvector_generators', 'linear_harmonic',
                 'pure_field_classification', 'unrestricted_mixed_eigenspace_axial_case')}
    exported.update({
        'claim':CLAIM,
        'linear_harmonic':'arbitrary; NO linear eigenvector condition',
        'pure_field_classification':'All pure points are one-sided; every nonzero semisimple axial C0 point is excluded by the written top-two-jet gate and direct pole units. Rank1 nilpotent/B0 import retained classifications.',
        'unrestricted_mixed_eigenspace_axial_C0_pure_case':'EXCLUDED ALL-ORDER, by the geometric proof; unconditional first8 ordinary ideal still uncomputed',
        'assurance':'Internal written proper-graph/period proof with exact local algebra controls; no external/formal verification. Direct rational pole units; no new CAS oracle.',
    })
    assert len(exported['quadratic_cone_generators']) == 7
    assert len(exported['cubic_contraction_generators']) == 3
    assert 'linear_eigenvector_generators' not in exported
    target = 'axial_contraction_zero_full_safe_container.json'
    (HERE/target).write_text(json.dumps(exported, indent=2)+'\n')
    receipt = {
        'claim':CLAIM, 'status':'PASS_ALGEBRAIC_CONTROLS',
        'checker_sha256':digest(Path(__file__).name),
        'model_sha256':digest('axial_contraction_zero.py'),
        'proof_sha256':digest('CUBIC_TOP_TWO_JET_GATE_AND_FULL_AXIAL_CONTRACTION.md'),
        'pole_units_sha256':digest('axial_contraction_zero_boundary_pole_units.json'),
        'closed_safe_container_sha256':digest(target),
        'boundary_charts':'Explicit conic chart, quadric identity, volume Jacobian, and top-two jets; reciprocal polynomial checks infinity.',
        'pole_unit_multiplier_term_counts':counts,
        'first3_solve':{'a*c':'-16/105','b^2':'-31/105'},
        'new_Singular_or_msolve_oracle':False,
        'analytic_dependency':'GPROPERHOLPERIOD20261115: internal written holomorphic-period theorem; not machine certified here.',
        'assurance':'The checker verifies algebra only. Smooth proper graph, finite critical set, volume extension and the period application are internal written geometric/analytic arguments. No independent specialist/external/formal verification.',
        'finite_prefix':'Six suffices on the two pole patterns. No unconditional full axial C0 consecutive-prefix bound asserted.',
        'mixed_axial_C0_all_order_pure_case':'EXCLUDED',
        'mixed_axial_C0_unconditional_first8_ordinary_unit':'UNCOMPUTED; failed large computations retained',
        'full_Image33':'OPEN', 'three_global_radical_memberships':'ALL OPEN',
        'elapsed_seconds':round(time.monotonic()-started,3),
    }
    (HERE/'cubic_top_two_jet_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()

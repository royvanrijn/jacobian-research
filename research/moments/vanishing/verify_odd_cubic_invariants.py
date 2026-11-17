"""Exact auxiliary algebra for the odd-cubic classification and prefix lower bounds.

This does not produce a 19-zero-moment coefficient witness or check the
retained analytic classification, Reynolds theorem, or Nullstellensatz.
"""
from pathlib import Path
import hashlib
import itertools
import json
import time
import sympy as s


def main():
    started = time.monotonic()
    path = Path(__file__).resolve()
    x = s.symbols('x y z')
    h = s.symbols('h0:7')
    l = s.Matrix(s.symbols('l0:3'))
    c = s.Matrix(s.symbols('c0:3'))
    H = (h[0]*x[0]**3+h[1]*x[0]**2*x[1]+h[2]*x[0]*x[1]**2
         +h[3]*x[1]**3+h[4]*x[0]**2*x[2]+h[5]*s.prod(x)
         +h[6]*x[1]**2*x[2]-(3*h[0]+h[2])*x[0]*x[2]**2
         -(h[1]+3*h[3])*x[1]*x[2]**2-(h[4]+h[6])*x[2]**3/3)
    lap = lambda f: s.expand(sum(s.diff(f, v, 2) for v in x))
    assert lap(H) == 0
    F = H+sum(v*v for v in x)*sum(c[i]*x[i] for i in range(3))
    assert s.expand(lap(F)-10*sum(c[i]*x[i] for i in range(3))) == 0
    tensor = [[[s.diff(H, x[i], x[j], x[k])/6 for k in range(3)]
               for j in range(3)] for i in range(3)]

    def contract(a, b, d):
        return s.expand(sum(tensor[i][j][k]*a[i]*b[j]*d[k]
                            for i, j, k in itertools.product(range(3), repeat=3)))

    w = l.cross(c)
    unsquared = [contract(l, l, w), contract(l, c, w), contract(c, c, w)]
    inv = s.Matrix([l.dot(l), l.dot(c), c.dot(c), contract(l, l, l),
                    contract(l, l, c), contract(l, c, c), contract(c, c, c)]
                   + [f*f for f in unsquared])
    variables = list(l)+list(c)+list(h)
    frame = dict(zip(list(l)+list(c), [1, 0, 0, 0, 1, 0]))
    point = frame | dict(zip(h, range(1, 8)))
    columns = [0, 1, 4, 6, 7, 8, 9, 10, 11, 12]
    J = inv.jacobian(variables)
    minor = s.factor(J[:, columns].subs(frame).det())
    assert minor == 8*h[4]*h[5]*h[6]/6561
    value = minor.subs(point)
    assert value == s.Rational(560, 2187)
    # Universal infinitesimal SO(3) invariance, plus one reflection, proves
    # O(3) invariance. The action on coefficients is H(x) -> H(R^T x).
    monomials = [x[0]**3, x[0]**2*x[1], x[0]*x[1]**2, x[1]**3,
                 x[0]**2*x[2], s.prod(x), x[1]**2*x[2]]
    infinitesimal_checks = 0
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        A = s.zeros(3); A[i, j] = 1; A[j, i] = -1
        dH = s.Poly(-sum((A*s.Matrix(x))[k]*s.diff(H, x[k]) for k in range(3)), x)
        tangent = list(A*l)+list(A*c)+[dH.coeff_monomial(m) for m in monomials]
        result = J*s.Matrix(tangent)
        assert all(s.expand(f) == 0 for f in result)
        infinitesimal_checks += len(result)
    R = s.diag(-1, 1, 1)
    reflected = s.Poly(H.subs({x[0]: -x[0]}), x)
    substitutions = dict(zip(list(l)+list(c), list(R*l)+list(R*c)))
    substitutions.update(dict(zip(h, [reflected.coeff_monomial(m) for m in monomials])))
    assert all(s.expand(f.xreplace(substitutions)-f) == 0 for f in inv)

    # Homogeneous restriction: the second vector is the covariant H(c,c,-).
    covector = s.Matrix([contract(c, c, s.eye(3)[:, i]) for i in range(3)])
    substitution = dict(zip(l, covector))
    assert s.expand(inv[5].subs(substitution, simultaneous=True)
                    -inv[0].subs(substitution, simultaneous=True)) == 0
    assert s.expand(inv[6].subs(substitution, simultaneous=True)
                    -inv[1].subs(substitution, simultaneous=True)) == 0
    assert s.expand(unsquared[2].subs(substitution, simultaneous=True)) == 0
    # Chain rule at an exact point avoids expanding degree-22 polynomials.
    hom_variables = list(c)+list(h)
    hom_point = dict(zip(list(c)+list(h), [0, 1, 0]+list(range(1, 8))))
    full_point = hom_point | {l[i]: covector[i].subs(hom_point) for i in range(3)}
    selected = [0, 1, 2, 3, 4, 7, 8]
    map_jacobian = s.Matrix(list(covector)+list(c)+list(h)).jacobian(hom_variables)
    hom_J = J[selected, :].subs(full_point)*map_jacobian.subs(hom_point)
    pivots = list(hom_J.rref()[1])
    assert len(pivots) == 7
    hom_minor = s.factor(hom_J[:, pivots].det())
    assert hom_minor != 0

    Z, W, T, a, b, cc, d, e, lam = s.symbols('Z W T a b cc d e lam')
    Zp = Z+lam*T-lam**2*W/2
    Tp = T-lam*W
    assert s.expand(2*Zp*W+Tp**2-(2*Z*W+T**2)) == 0
    P = Zp*(a*Zp**2+b*Zp*Tp+cc*Tp**2+d*Zp*W+e)
    polys = s.Poly(P, Z, W, T)
    coef = s.Matrix(polys.coeffs())
    CJ = coef.jacobian([a, b, cc, d, e, lam]).subs({lam: 0, a: 1, b: 2, cc: 3, d: 4, e: 1})
    row_pivots = list(CJ.T.rref()[1])
    assert len(row_pivots) == 6
    normal_minor = CJ[row_pivots, :].det()
    assert normal_minor != 0
    base = s.Poly(P.subs(lam, 0), Z, W, T)
    weights = sorted({powers[0]-powers[1] for powers, _ in base.terms()})
    assert weights == [1, 2, 3]

    def wick(f):
        out = 0
        for (i, j, k), coefficient in s.Poly(s.expand(f), Z, W, T).terms():
            if i == j and k % 2 == 0:
                out += coefficient*s.factorial(i)*s.factorial2(k-1)
        return s.expand(out)

    sharpness = []
    for degree in range(1, 5):
        actual = wick(W**degree*(Z*T*T)**degree)
        assert actual == s.factorial(degree)*s.factorial2(2*degree-1)
        assert wick(W**degree*(Z*T*T)**(degree+1)) == 0
        sharpness.append({'D': degree, 'M_D': str(actual), 'M_D_plus_1': '0'})
    sample = P.subs({lam: s.Rational(2, 3), a: 1, b: -2, cc: 3, d: 1, e: 2})
    pure_sample = [wick(sample**k) for k in range(1, 7)]
    assert pure_sample == [0]*6
    nonnull = [wick((W+Z**2*T)**k) for k in range(1, 7)]
    assert nonnull == [0]*5+[360]
    report = {
        'completed': True, 'checker_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'sympy_version': s.__version__, 'runtime_seconds': round(time.monotonic()-started, 3),
        'harmonic_trace_factor': 10,
        'odd_independent_invariants': 10,
        'odd_invariant_degrees': [2, 2, 2, 4, 4, 4, 4, 10, 10, 10],
        'odd_Jacobian_columns': columns, 'odd_Jacobian_minor': str(minor),
        'odd_Jacobian_minor_at_point': str(value),
        'universal_Lie_derivative_checks': infinitesimal_checks,
        'universal_reflection_checks': 10,
        'homogeneous_independent_invariants': 7,
        'homogeneous_invariant_indices': selected,
        'homogeneous_Jacobian_columns': pivots,
        'homogeneous_Jacobian_minor_at_point': str(hom_minor),
        'homogeneous_covariant_at_point': [str(v.subs(hom_point)) for v in covector],
        'normal_form_Jacobian_rank': 6, 'normal_form_Jacobian_minor': str(normal_minor),
        'normal_form_positive_circular_weights': weights,
        'sharp_mixed_examples': sharpness,
        'normal_form_sample_pure_M1_to_M6': [str(v) for v in pure_sample],
        'nonnull_sample_pure_M1_to_M6': [str(v) for v in nonnull],
        'proposed_global_lower_bound': 20,
        'proposed_homogeneous_lower_bound': 14,
        'explicit_19_zero_moment_witness_computed': False,
        'explicit_13_zero_moment_homogeneous_witness_computed': False,
        'M20_nonzero_proved': False, 'prefix20_sufficiency_proved': False,
        'classification_and_invariant_dimension_argument_checked_by_machine': False,
        'boundary': 'Auxiliary exact algebra only. The odd classification uses retained internal analytic theorems; Reynolds, finite generation, the finite-module argument and algebraic-point existence remain written mathematics. No independent replay, external review or formal verification.'
    }
    output = path.with_name('odd_cubic_invariants_verification.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

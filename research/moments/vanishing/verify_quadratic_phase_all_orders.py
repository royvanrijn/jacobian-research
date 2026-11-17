"""Exact two-jet classifier plus small controls for the written all-order proof."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def classify_two_jet(a0, b, a, phase, variables):
    """Exact rational/polynomial-parameter equations; parameter strata are not solved."""
    n = len(variables)
    phase = s.expand(phase)
    assert phase == 0 or s.Poly(phase, *variables).total_degree() <= 2
    b, a = s.Matrix(b), s.Matrix(a)
    assert a.shape == (n, n) and a == a.T and b.shape == (n, 1)
    assert not any(v in s.sympify(c).free_symbols for v in variables for c in [a0, *b, *a])
    parameters = set().union(*(s.sympify(c).free_symbols for c in [a0, *b, *a, phase]))-set(variables)
    generators = tuple(variables)+tuple(sorted(parameters, key=str))
    for coefficient in [a0, *b, *a, phase]:
        s.Poly(coefficient, *generators, domain=s.QQ)
    h = s.hessian(phase, variables)
    ell = s.Matrix([s.diff(phase, v).subs(dict.fromkeys(variables, 0)) for v in variables])
    equations = [s.expand(c) for c in [a0, *(h*b), (b.T*ell)[0], *((a*h)**n)] if c != 0]
    equations = list(dict.fromkeys(c for c in equations if c != 0))
    pure = phase == 0 or not equations
    return {'pure_at_given_exact_coefficients': pure,
            'zero_phase_branch': phase == 0,
            'nonzero_phase_branch_equations': [str(c) for c in equations],
            'pure_prefix_sufficient': n+1,
            'safe_if_pure': True,
            'common_mixed_slope': 3*(n+1),
            'higher_operator_coefficients_unrestricted': True}


def classify_symbol(symbol, phase, variables, derivatives):
    assert len(variables) == len(derivatives)
    assert set(variables).isdisjoint(s.sympify(symbol).free_symbols), 'Symbol coefficients must be constant in x'
    poly = s.Poly(symbol, *derivatives)
    n = len(variables)
    b, a = s.zeros(n, 1), s.zeros(n)
    for i in range(n):
        b[i] = poly.coeff_monomial(derivatives[i])
        a[i, i] = poly.coeff_monomial(derivatives[i]**2)
        for j in range(i):
            a[i, j] = a[j, i] = poly.coeff_monomial(derivatives[i]*derivatives[j])/2
    return classify_two_jet(poly.coeff_monomial(1), b, a, phase, variables)


def apply_symbol(symbol, polynomial, variables, derivatives):
    answer = 0
    for mon, coeff in s.Poly(symbol, *derivatives).terms():
        term = polynomial
        for v, order in zip(variables, mon):
            if order:
                term = s.diff(term, v, order)
        answer += coeff*term
    return s.expand(answer)


def run():
    x, y, z, w = s.symbols('x y z w')
    dx, dy, dz, dw = s.symbols('dx dy dz dw')
    t, alpha, beta = s.symbols('t alpha beta')
    zero = lambda f: s.expand(f) == 0
    def iterate(symbol, phase, variables, derivatives, m):
        for _ in range(m):
            phase = apply_symbol(symbol, phase, variables, derivatives)
        return phase
    cases = [
        ('higher_order_only', dx**3+dy**5, x*y+x+1, (x, y), (dx, dy), True),
        ('drift_quadratic', dy+dx**3, x**2, (x, y), (dx, dy), True),
        ('free_higher_coefficients', dy+alpha*dx**3+beta*dx**7, x**2, (x, y), (dx, dy), True),
        ('linear_phase_higher_jets', dy+dx**2+dx**5, x+1, (x, y), (dx, dy), True),
        ('first_moment_cancellation', dy+dx**2+dx**3+dx*dy, x**2-2*y, (x, y), (dx, dy), False),
        ('constant_symbol', 1+dx**3, x**2, (x, y), (dx, dy), False),
        ('zero_phase_unit', 1+dx**3, s.Integer(0), (x, y), (dx, dy), True),
        ('four_variable_nilpotent', 2*dx*dz+2*dy*dw+dx**3+dz**3+dy**2*dw**2,
         x*w+y**2/2, (x, y, z, w), (dx, dy, dz, dw), True),
    ]
    rows = []
    for label, symbol, phase, variables, derivatives, expected in cases:
        result = classify_symbol(symbol, phase, variables, derivatives)
        assert result['pure_at_given_exact_coefficients'] == expected
        tests = [iterate(symbol, phase**m, variables, derivatives, m) == 0
                 for m in range(1, len(variables)+2)]
        assert all(tests) == expected
        result.update(case=label, first_prefix_controls=tests)
        rows.append(result)

    # Pure truncation equality needs the first equation,not merely degree<=2.
    bad_phase = x*y
    symbol = dy+dx**2*dy
    difference = iterate(symbol, bad_phase**2, (x, y), (dx, dy), 2)-iterate(dy, bad_phase**2, (x, y), (dx, dy), 2)
    assert difference == 8
    assert apply_symbol(symbol, bad_phase, (x, y), (dx, dy)) == x

    sharp = []
    for d in range(1, 4):
        value = iterate(dy+dx**3, y**d*x**(6*d), (x, y), (dx, dy), 3*d)
        expected = s.factorial(3*d)*s.factorial(6*d)/s.factorial(2*d)
        assert value == expected
        assert iterate(dy+dx**3, y**d*x**(2*(3*d+1)), (x, y), (dx, dy), 3*d+1) == 0
        assert iterate(dy, y**d*x**(6*d), (x, y), (dx, dy), 3*d) == 0
        sharp.append({'d': d, 'boundary_m': 3*d, 'marked_value': str(value), 'next_marked_zero': True})

    # Explicit scalar potential for a nontrivial nilpotent quadratic phase.
    variables = (x, y, z, w)
    a = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    p = x*w+y**2/2+x+z+1
    h = s.hessian(p, variables)
    ell = s.Matrix([1, 0, 1, 0])
    j = a*h
    assert j**4 == s.zeros(4)
    resolvent = sum((t**k*j**k for k in range(4)), s.zeros(4))
    xv = s.Matrix(variables)
    u = s.expand((xv.T*h*resolvent*xv)[0]/2+(ell.T*resolvent*xv)[0]+1
                 +t*(ell.T*resolvent*a*ell)[0]/2)
    gu = s.Matrix([s.diff(u, v) for v in variables])
    assert zero(s.diff(u, t)-(gu.T*a*gu)[0]/2)
    assert s.degree(u, t) <= 4
    assert u.subs(t, 0) == p
    mixed_prefix = []
    symbol = 2*dx*dz+2*dy*dw+dx**3+dz**3
    for m in (1, 2, 3, 4, 5):
        value = iterate(symbol, z*p**m, variables, (dx, dy, dz, dw), m)
        mixed_prefix.append({'m': m, 'zero': value == 0})

    # Arithmetic control of the bounded number of higher-weight factors.
    factor_controls = 0
    for ntime in range(5):
        for d in range(1, 5):
            m = 3*(ntime+1)*d+1
            for r in range(2*d+1):
                assert m-r-d > ntime*(r+d)
                factor_controls += 1
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__, 'classifier_cases': rows,
            'first_equation_truncation_mutation_difference': str(difference),
            'sharp_higher_jet_cutoff_controls': sharp,
            'nilpotent_quadratic_scalar_formula_HJ_checked': True,
            'nilpotent_quadratic_scalar_time_degree': int(s.degree(u, t)),
            'nontrivial_mixed_prefix_controls': mixed_prefix,
            'higher_factor_inequality_controls': factor_controls,
            'finite_symbol_two_jet_classifier_implemented': True,
            'symbolic_parameter_equations_returned_without_solving': True,
            'universal_proof_inferred_from_finite_controls': False,
            'general_higher_degree_GVC_solver_implemented': False,
            'general_Gaussian_uniformity_proved': False,
            'ordinary_nonreduced_ideal_answer': False,
            'formal_nonpolynomial_phases_covered': False,
            'higher_symbol_terms_dropped_from_mixed_contractions': False,
            'Dvorsky_replayed_again': False,
            'protected_ternary_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('quadratic_phase_all_orders_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

"""Small exact controls of coefficient jets and marked infinity levels."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def run():
    x, t, a, u, v = s.symbols('x t a u v')
    c0, c1, c2, c3 = s.symbols('c0 c1 c2 c3')
    phase = c0+c1*x+c2*x**2+c3*x**3
    def gauss(f):
        value = 0
        for (k,), coeff in s.Poly(s.expand(f), x).terms():
            if k % 2 == 0:
                value += coeff*(s.factorial2(k-1) if k else 1)
        return s.expand(value)
    jet_controls = 0
    for k in range(1, 4):
        for m in range(k, 5):
            left = s.diff(gauss(phase**m)/s.factorial(m), c1, k)
            right = gauss(x**k*phase**(m-k))/s.factorial(m-k)
            assert s.expand(left-right) == 0
            jet_controls += 1

    # Genuine retained pure quartic family;these are only finite transcription controls.
    z, w, q, e = s.symbols('Z W T e')
    pe = w*(1+e*z)-e*(1+3*e*z/2+e**2*z**2/2)*q**2
    def wick(f):
        answer = 0
        for (az, aw, aq), coeff in s.Poly(s.expand(f), z, w, q).terms():
            if az == aw and aq % 2 == 0:
                answer += coeff*s.factorial(az)*(s.factorial2(aq-1) if aq else 1)
        return s.expand(answer)
    quartic_rows = []
    for m in range(1, 5):
        pure = wick(pe**m)
        marked = wick(z*pe**m)
        assert pure == s.diff(pure, e) == 0
        assert s.expand(marked-s.factorial(m)*e**(m-1)) == 0
        ambient = wick(s.diff((pe+v*z)**(m+1), v).subs(v, 0))/s.factorial(m+1)
        assert s.expand(ambient-e**(m-1)) == 0
        quartic_rows.append({'moment': m, 'pure': '0', 'marked_Z': str(marked),
                             'ambient_linear_jet_time_index': m+1})

    # Source Algorithm9's integral root can vary across central parameters.
    root_rows = []
    for c in (2, 1, 0, -1, -3):
        nonnegative = [root for root in (0, 1-c) if root >= 0]
        root_rows.append({'central_c': c, 'max_nonnegative_root_of_s_times_s_plus_c_minus1': max(nonnegative)})

    # Same abstract cyclic O module,unboundedly shifted marked generators.
    level_rows = []
    for n in (0, 1, 3, 5):
        assert s.diff(t**n, t, n) == s.factorial(n)
        germ = u**(-n)
        assert s.expand(u*s.diff(germ, u)+germ-(1-n)*germ) == 0
        assert s.expand(u**n*germ) == 1
        if n:
            assert s.expand(u**(n-1)*germ) == 1/u
        level_rows.append({'generator_degree': n, 'V_index_for_sharp_h': 1-n})

    # A regular-at-infinity chart with a nonpolynomial normal derivative.
    regular_rows = []
    for n in (0, 2, 4):
        exponent = t/(1-t)
        transformed = s.cancel(exponent.subs(t, 1/u))
        assert transformed == 1/(u-1)
        matrix_entry = s.diff(transformed, u)
        assert s.limit(matrix_entry, u, 0) == -1
        period = t**n+a*s.exp(exponent)
        assert period.subs(a, 0) == t**n
        assert s.diff(period, a) == s.exp(exponent)
        assert s.cancel(s.diff(s.exp(exponent), t)/s.exp(exponent)-1/(1-t)**2) == 0
        assert [s.cancel(u**n*a), s.cancel(u**n*u**(-n))] == [a*u**n, 1]
        regular_rows.append({'N': n, 'infinity_connection_ordinary': True,
                             'polynomial_fiber_degree': n, 'normal_derivative_nonpolynomial_by_double_log_pole': True})

    # Tail generator evaluation does not annihilate its parameter-D-module closure.
    tail_rows = []
    for n in (2, 4):
        f = a*s.exp(t)+t**n
        tail = s.diff(f, t, n+1)
        assert tail.subs(a, 0) == 0
        assert s.diff(tail, a) == s.exp(t)
        tail_rows.append({'N': n, 'tail_restricts_zero': True,
                          'normal_derivative_of_tail_restricts_exp_t': True})

    # Non-characteristic for M alone does not establish the localization hypothesis.
    fiber_normal = s.Matrix([0, 1])
    divisor_normal = s.Matrix([1, -1])
    infinity_normal = s.Matrix([1, 0])
    assert s.Matrix.hstack(fiber_normal, divisor_normal).det() != 0
    assert infinity_normal-divisor_normal == fiber_normal
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__, 'ambient_linear_coefficient_jet_controls': jet_controls,
            'quartic_pure_line_and_normal_jet_controls': quartic_rows,
            'central_parameter_integral_root_controls': root_rows,
            'marked_generator_V_level_controls': level_rows,
            'regular_infinity_chart_controls': regular_rows,
            'finite_length_tail_mutation_controls': tail_rows,
            'localized_noncharacteristic_condition_distinguished': True,
            'universal_ambient_jet_proof_inferred_from_samples': False,
            'quartic_all_order_identity_reproved': False,
            'comprehensive_restriction_producer_executed': False,
            'general_V_filtration_or_Hom_producer_executed': False,
            'Gaussian_global_uniformity_proved': False,
            'irregular_family_restriction_theorem_imported': False,
            'finite_isomorphism_classes_used_as_generator_degree_bound': False,
            'ordinary_nonreduced_uniformity_proved': False,
            'protected_ternary_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('coefficient_restriction_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

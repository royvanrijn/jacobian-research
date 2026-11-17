"""Finite exact controls; no general inverse producer or new witness search."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def run():
    x, y, z, t, h = s.symbols('x y z t h')
    zero = lambda f: s.expand(f) == 0

    def grad(f, xs):
        return s.Matrix([s.diff(f, v) for v in xs])

    def lap(f, a, xs):
        return s.expand(sum(a[i, j]*s.diff(f, u, v)
                            for i, u in enumerate(xs)
                            for j, v in enumerate(xs)))

    def drift(f, b, xs):
        return s.expand(sum(v*s.diff(f, u) for v, u in zip(b, xs)))

    def op(f, a, b, xs, c=0):
        return s.expand(lap(f, a, xs)+drift(f, b, xs)+c*f)

    def iterate(f, a, b, xs, m, c=0):
        for _ in range(m):
            f = op(f, a, b, xs, c)
        return f

    sharp = []
    a2 = s.diag(1, 0)
    b2 = s.Matrix([0, 1])
    xs2 = (x, y)
    for d in range(1, 5):
        value = iterate(y**d*x**(2*d), a2, b2, xs2, 2*d)
        expected = s.factorial(2*d)**2/s.factorial(d)
        assert value == expected
        assert iterate(y**d*x**(2*d+1), a2, b2, xs2, 2*d+1) == 0
        assert iterate(y**d*x**(2*d), a2, s.zeros(2, 1), xs2, 2*d) == 0
        sharp.append({'d': d, 'drift_boundary_m': 2*d,
                      'marked_value': str(value), 'next_marked_zero': True,
                      'homogeneous_quadratic_value_at_same_index': '0'})

    # A first-moment cancellation is not separate quadratic and drift nullity.
    bad = x**2-2*y
    assert op(bad, a2, b2, xs2) == 0
    assert drift(bad, b2, xs2) == -2
    assert iterate(bad**2, a2, b2, xs2, 2) == 16
    ubad = x**2/(1-2*t)-2*y
    fbad = s.Matrix([x*(1-2*t), y])
    jbad = a2*s.hessian(bad, xs2)
    left = op(ubad, a2, b2, xs2).subs(dict(zip(xs2, fbad)), simultaneous=True)
    right = s.trace(jbad*(s.eye(2)-t*jbad).inv())+drift(bad, b2, xs2)
    assert s.cancel(left-right) == 0

    # Drift quotient with all mixed quadratic blocks present.
    xs = (x, y, z)
    a = s.Matrix([[0, 1, 2], [1, 0, 3], [2, 3, 5]])
    b = s.Matrix([0, 0, 1])
    p = x**3+x+y+1
    u = (x+t)**3+x+y+1+t
    g = s.Matrix(xs)+t*a*grad(u, xs)
    forward = s.Matrix(xs)-t*a*grad(p, xs)
    def compose(fs, gs):
        return s.Matrix([s.expand(f.subs(dict(zip(xs, gs)), simultaneous=True)) for f in fs])
    assert compose(forward, g) == s.Matrix(xs)
    assert compose(g, forward) == s.Matrix(xs)
    assert drift(p, b, xs) == drift(u, b, xs) == 0
    assert op(u, a, b, xs) == 0
    assert zero(s.diff(u, t)-(grad(u, xs).T*a*grad(u, xs))[0]/2)
    jmat = a*s.hessian(p, xs)
    assert jmat**3 == s.zeros(3)
    assert compose([op(u, a, b, xs)], forward)[0] == 0
    assert all(s.trace(jmat**k) == 0 for k in range(1, 5))
    quotient_controls = 0
    a_reduced = s.Matrix([[0, 1], [1, 0]])
    for k in (1, 2):
        for m in range(5):
            full = iterate(p**(m+k), a, b, xs, m)
            reduced = iterate(p**(m+k), a_reduced, s.zeros(2, 1), xs2, m)
            assert zero(full-reduced)
            coefficient = s.expand(u**k).coeff(t, m)
            expected = s.factorial(k)*full/(2**m*s.factorial(m)*s.factorial(m+k))
            assert zero(coefficient-expected)
            quotient_controls += 1
    for m in range(1, 5):
        assert iterate(p**m, a, b, xs, m) == 0

    # The heat-log degree lemma survives drift, including a nonpure cancellation.
    heat_rows = []
    for label, mat, vec, phase, pure in [
        ('linear_phase_sharp_drift', a2, b2, x, True),
        ('first_moment_cancellation', a2, b2, bad, False),
    ]:
        gp = grad(phase, xs2)
        def conjugated(f):
            return s.expand(op(f, mat, vec, xs2)+2*h*(gp.T*mat*grad(f, xs2))[0]
                            +(h*op(phase, mat, vec, xs2)+h**2*(gp.T*mat*gp)[0])*f)
        terms = [s.Integer(1)]
        raw = s.Integer(1)
        for m in range(1, 4):
            raw = conjugated(raw)
            terms.append(s.expand(raw/(2**m*h**m*s.factorial(m))))
        logs = [h*phase, terms[1], s.expand(terms[2]-terms[1]**2/2),
                s.expand(terms[3]-terms[1]*terms[2]+terms[1]**3/3)]
        powers = [int(term.as_powers_dict().get(h, 0))
                  for poly in logs for term in s.expand(poly).as_ordered_terms() if term != 0]
        assert all(power <= 1 for power in powers)
        assert all(power == 1 for power in powers) == pure
        heat_rows.append({'case': label, 's_degree_at_most_one': True,
                          'nonpositive_s_powers_present': not pure})

    # Characteristic-zero trace lemma: trace2 alone cannot establish nilpotence.
    deceptive = s.diag(1, s.I)
    assert s.trace(deceptive**2) == 0
    assert s.trace(deceptive**3) != 0
    trace_rows = []
    for n in range(1, 5):
        nil = s.zeros(n)
        for i in range(n-1):
            nil[i, i+1] = 1
        assert nil**n == s.zeros(n)
        assert all(s.trace(nil**k) == 0 for k in range(2, n+2))
        trace_rows.append({'dimension': n, 'trace_powers_checked': list(range(2, n+2))})

    # A nonzero constant symbol coefficient cannot have a nonzero polynomial kernel.
    test = x**3+2*y**2+1
    value = op(test, a2, b2, xs2, c=1)
    assert s.Poly(value, x, y).coeff_monomial(x**3) == 1
    assert value != 0
    assert iterate(s.Integer(0), a2, b2, xs2, 1, c=1) == 0
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__, 'sharp_drift_cutoff_controls': sharp,
            'first_moment_cancellation_rejected_at_second_moment': True,
            'rational_trace_identity_checked': True,
            'quotient_with_cross_terms_both_inverse_compositions': True,
            'quotient_and_shifted_power_controls': quotient_controls,
            'Laurent_heat_log_controls': heat_rows,
            'trace_lemma_controls': trace_rows,
            'trace_square_only_mutation_rejected': True,
            'constant_symbol_branch_checked': True,
            'all_order_trace_lemma_inferred_from_samples': False,
            'universal_inverse_or_closed_locus_producer_executed': False,
            'general_Gaussian_uniformity_proved': False,
            'ordinary_nonreduced_ideal_answer': False,
            'operators_of_order_at_least_three_classified': False,
            'higher_degree_quadraticization_excluded': False,
            'external_quartic_witness_replayed': False,
            'protected_ternary_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('order_two_drift_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))

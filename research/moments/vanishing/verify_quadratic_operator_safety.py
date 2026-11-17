"""Small exact controls for the written quadratic-operator theorem."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def run():
    x, y, z, t, h, e = s.symbols('x y z t h e')
    xs = (x, y)
    expand = s.expand
    zero = lambda f: expand(f) == 0

    def grad(f, variables=xs):
        return s.Matrix([s.diff(f, v) for v in variables])

    def lap(f, a, variables=xs):
        return expand(sum(a[i, j]*s.diff(f, u, v)
                          for i, u in enumerate(variables)
                          for j, v in enumerate(variables)))

    def it(f, a, m):
        for _ in range(m):
            f = lap(f, a)
        return f

    def comp(fs, gs):
        return s.Matrix([expand(f.subs(dict(zip(xs, gs)), simultaneous=True)) for f in fs])

    # A polynomial family crosses rank2 to rank0 at e0; all equations are in Q[e].
    a = s.Matrix([[0, e], [e, 0]])
    p0 = x**3
    c = s.Matrix([1, 1])
    p = p0+x+y+1
    u = (x+t*e)**3+x+y+1+t*e
    g = s.Matrix(xs)+t*a*grad(u)
    forward = s.Matrix(xs)-t*a*grad(p)
    assert all(zero(f-v) for f, v in zip(comp(forward, g), xs))
    assert all(zero(f-v) for f, v in zip(comp(g, forward), xs))
    assert zero(s.diff(u, t)-(grad(u).T*a*grad(u))[0]/2)
    gp = comp(grad(p), g)
    assert zero(u-comp([p], g)[0]+t*(gp.T*a*gp)[0]/2)
    assert zero(lap(u, a))
    assert (a*s.hessian(p, xs))**2 == s.zeros(2)
    # The scalar term te cannot be discarded even though the projected gradient agrees.
    wrong = u-t*e
    assert a*grad(wrong) == a*grad(u)
    assert not zero(s.diff(wrong, t)-(grad(wrong).T*a*grad(wrong))[0]/2)
    assert g.subs(e, 0) == s.Matrix(xs)
    assert u.subs(e, 0) == p
    controls = 0
    for k in (1, 2):
        for j in range(0, 3*k+2):
            rhs = s.factorial(k)*it(p**(j+k), a, j)/(2**j*s.factorial(j)*s.factorial(j+k))
            assert zero(expand(u**k).coeff(t, j)-rhs)
            controls += 1
    for j in range(1, 5):
        assert zero(it(p**j, a, j))
        controls += 1

    # Degenerate rank1 operator: a time-dependent kernel potential is essential.
    singular = s.diag(1, 0)
    ps = x*y+z*y+z**2+1
    us = ps+t*y**2/2
    gs = s.Matrix([x+t*y, y])
    fs = s.Matrix(xs)-t*singular*grad(ps)
    assert comp(fs, gs) == s.Matrix(xs)
    assert comp(gs, fs) == s.Matrix(xs)
    assert singular*grad(us) == singular*grad(ps)
    assert not zero(s.diff(ps, t)-(grad(ps).T*singular*grad(ps))[0]/2)
    assert zero(s.diff(us, t)-(grad(us).T*singular*grad(us))[0]/2)
    for j in range(1, 5):
        assert zero(it(ps**j, singular, j))
        controls += 1

    # Laurent heat logarithm via the conjugated heat operator, not guessed cumulants.
    log_rows = []
    for label, mat, phase, pure in [
        ('rank_drop_affine', a, p, True),
        ('singular_kernel', singular, ps, True),
        ('harmonic_non_HN', s.eye(2), x**3-3*x*y**2+2*x+1, False),
        ('nonharmonic', s.eye(2), x**2+y+1, False),
    ]:
        gp = grad(phase)
        def conjugated(f):
            return expand(lap(f, mat)+2*h*(gp.T*mat*grad(f))[0]
                          +(h*lap(phase, mat)+h**2*(gp.T*mat*gp)[0])*f)
        terms = [s.Integer(1)]
        raw = s.Integer(1)
        for j in range(1, 4):
            raw = conjugated(raw)
            terms.append(expand(raw/(2**j*h**j*s.factorial(j))))
        logs = [h*phase, terms[1],
                expand(terms[2]-terms[1]**2/2),
                expand(terms[3]-terms[1]*terms[2]+terms[1]**3/3)]
        powers = []
        for term in logs:
            exponents = [int(mon.as_powers_dict().get(h, 0)) for mon in expand(term).as_ordered_terms() if mon != 0]
            assert all(power <= 1 for power in exponents)
            powers += exponents
        if pure:
            assert all(power == 1 for power in powers)
        else:
            assert any(power <= 0 for power in powers)
        jmat = mat*s.hessian(phase, xs)
        tests = [zero(it(phase**j, mat, j)) for j in (1, 2)]
        assert all(tests) == pure
        assert (jmat**2 == s.zeros(2)) == pure
        # The s^1 coefficient solves Hamilton--Jacobi even in the non-HN controls.
        potentials = [expand(term).coeff(h, 1) for term in logs]
        for j in range(3):
            rhs = sum((grad(potentials[k]).T*mat*grad(potentials[j-k]))[0]
                      for k in range(j+1))/2
            assert zero((j+1)*potentials[j+1]-rhs)
        log_rows.append({'case': label, 'pure_first_two': tests,
                         's_degree_at_most_one': True,
                         'nonpositive_s_power_present': not pure,
                         'HJ_coefficient_controls': 3})

    # The common cutoff increases with active rank; these only check arithmetic.
    cutoff_rows = []
    for d in (2, 3, 4, 5):
        previous = 0
        for n in range(1, 6):
            b = (d-1)**(n-1)
            bound = 2*n*b-n+1
            assert bound > previous
            cutoff_rows.append({'n': n, 'D': d, 'C': bound})
            previous = bound
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__,
            'affine_rank_drop_both_inverse_compositions': True,
            'full_scalar_affine_term_mutation_rejected': True,
            'singular_kernel_potential_mutation_rejected': True,
            'pure_and_shifted_contraction_controls': controls,
            'Laurent_heat_log_controls': log_rows,
            'cutoff_arithmetic_controls': cutoff_rows,
            'universal_Laurent_log_proof_inferred_from_samples': False,
            'general_inverse_or_closed_locus_producer_executed': False,
            'general_Gaussian_uniformity_proved': False,
            'ordinary_nonreduced_ideal_answer': False,
            'operators_with_lower_order_terms_covered': False,
            'higher_degree_quadraticization_excluded': False,
            'external_quartic_witness_replayed': False,
            'protected_ternary_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('quadratic_operator_safety_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
